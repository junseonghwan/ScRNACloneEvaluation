rm(list=ls())

library(dplyr)
library(GenomicRanges)
library(ScRNAClone)
library(sabre)

snv_file <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/snv.txt"
annotation_file <- "/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs.csv"
gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/ground_truth.txt"

snv <- read.table(snv_file, header=T)
dim(snv)
snv.gr <- ConstructGranges(snv$CHR, snv$POS, width = 0)

laks_snv <- read.table(annotation_file, header=T, sep=",")
laks_snv.gr <- ConstructGranges(laks_snv$chrom, laks_snv$coord, width = 0)

overlaps <- findOverlaps(laks_snv.gr, snv.gr)
gt <- laks_snv[overlaps@from,]
head(gt)

# For each SNV, get the names of clones.
gt$loc <- paste(gt$chrom, gt$coord, sep=":")
locs <- unique(gt$loc)
n_snvs <- length(locs)
snv$clone_name <- rep("", n_snvs)
for (i in 1:n_snvs) {
    temp <- subset(gt, loc == locs[i])
    snv$clone_name[i] <- paste(temp$clone_id[temp$is_present==1], collapse = "_")
}
dim(snv)
unique(snv$clone_name)
snv[which(snv$clone_name == ""),]

# The tree is as follows:
# (((A,B), (C,D)), ((E,F), ((G,H),I))).
# We will remove any annotation that is not consistent with the tree.
valid_clusters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                    "A_B", "C_D", "E_F", "G_H",
                    "G_H_I", "A_B_C_D", "E_F_G_H_I",
                    "A_B_C_D_E_F_G_H_I")
snv_with_valid_clones <- snv[snv$clone_name %in% valid_clusters,]
dim(snv_with_valid_clones)
unique(snv_with_valid_clones$clone_name)
snv$clone_name

write.table(x = snv_with_valid_clones, file = gt_file, row.names = F, col.names = T, quote = F)

