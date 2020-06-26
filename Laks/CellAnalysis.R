library(dplyr)
library(GenomicRanges)
library(matrixStats)
library(ScRNAClone)
library(TailRank)

SNV_PATH <- "/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs_exon.csv"
SC_READS_PATH <- "/Users/seonghwanjun/data/cell-line/smartseq3/ReadsExon/"
annotation_file <- "/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs.csv"

# Read in single cell data.
sc <- CombineSingleCellReads(SC_READS_PATH)
length(unique(sc$Cell))

# Read in Laks data.
exon_snv <- read.table(SNV_PATH, header=T)
laks_snv <- read.table(annotation_file, header=T, sep=",")
laks_snv$locs <- paste(laks_snv$chrom, laks_snv$coord, sep=":")
locs <- unique(laks_snv$locs)
clone_names <- rep("NA", length(locs))
for (i in 1:length(locs)) {
    muts <- laks_snv[laks_snv$locs %in% locs[i],]
    clone_names[i] <- paste(muts[muts$is_present == 1,"clone_id"], collapse = "_")
}

library(plyr)
ret <- ldply(strsplit(locs, ":"), function(x) {
    cbind("CHR"=x[1], "POS"=x[2])
})
ret$CHR <- as.character(ret$CHR)
ret$POS <- as.numeric(as.character(ret$POS))
ret$clone_name <- clone_names
snv2clone <- ret
snv2clone.gr <- ConstructGranges(snv2clone$CHR, snv2clone$POS, width = 0)
head(snv2clone)
detach("package:plyr", unload=TRUE)

ret <- sc %>% group_by(Cell) %>% summarise(muts = ID[which(d-a > 0)])

# Get mutation profile for each cell.
cell_ids <- unique(sc$Cell)
n_cells <- length(cell_ids)
clones <- unique(clone_names)
clones <- clones[order(nchar(clones), clones)]
n_clones <- length(clones)
mut_profile <- matrix(0, nrow = n_cells, ncol = n_clones)
clones[clones == ""] <- "NA"
mut_names <- rep("", n_cells)
for (i in 1:n_cells) {
    cell_id <- cell_ids[i]
    mut_ids <- subset(ret, Cell == cell_id)$muts
    if (length(mut_ids) == 0)
        next
    snvs <- exon_snv[exon_snv$ID %in% mut_ids,]
    snvs.gr <- ConstructGranges(snvs$CHR, snvs$POS, width = 0)
    overlap <- findOverlaps(snvs.gr, snv2clone.gr)
    mut_status <- unique(snv2clone[overlap@to,"clone_name"])
    mut_profile[i, clones %in% mut_status] <- 1
    mut_names[i] <- paste(mut_status, collapse = ",")
}
mut_profile.df <- as.data.frame(mut_profile)
names(mut_profile.df) <- clones
colSums(mut_profile.df)
strsplit(mut_names, ",")

valid_clusters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                    "A_B", "C_D", "E_F", "G_H",
                    "G_H_I", "A_B_C_D", "E_F_G_H_I",
                    "A_B_C_D_E_F_G_H_I")
valid.mut.df <- mut_profile.df[,names(mut_profile.df) %in% valid_clusters]
colSums(valid.mut.df)
valid.mut.df[which(valid.mut.df$A_B == 1),]
valid.mut.df[which(valid.mut.df$C_D == 1),]
valid.mut.df[rowSums(valid.mut.df) == 2,]

invalid.mut.df <- mut_profile.df[,!(names(mut_profile.df) %in% valid_clusters)]
colSums(invalid.mut.df)

# Make a plot similar to Figure 3D in Laks but from the exons.
laks_snv.gr <- ConstructGranges(laks_snv$chrom, laks_snv$coord, width = 0)
exon_snv.gr <- ConstructGranges(exon_snv$CHR, exon_snv$POS, width = 0)
overlaps <- findOverlaps(laks_snv.gr, exon_snv.gr)
laks_snv.exon <- laks_snv[overlaps@from,]

#snv2clone2 <- snv2clone[order(snv2clone$clone_name, nchar(snv2clone$clone_name), decreasing = TRUE),]
snv2clone2 <- snv2clone[order(nchar(snv2clone$clone_name), snv2clone$clone_name, decreasing = F),]
snv2clone2$locs <- paste(snv2clone2$CHR, snv2clone2$POS, sep=":")
laks_snv.exon$locs <- factor(laks_snv.exon$locs, levels = snv2clone2$locs)
library(ggplot2)
p <- ggplot(laks_snv.exon) + geom_tile(aes(x = clone_id, y = locs, fill = as.factor(is_present)))
p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p
length(unique(laks_snv.exon$locs))
