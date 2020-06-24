rm(list=ls())

library(dplyr)
library(GenomicRanges)
library(ScRNAClone)
library(sabre)
all_snv <- read.table("/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs.csv", header=T, sep=",")
exon_snv <- read.table("/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs_exon.csv", header=T)

all_snv.gr <- ConstructGranges(all_snv$chrom, all_snv$coord, width = 0)
exon_snv.gr <- ConstructGranges(exon_snv$CHR, exon_snv$POS, width = 0)
overlaps <- findOverlaps(all_snv.gr, exon_snv.gr)
gt <- all_snv[overlaps@from,]
head(gt)

# For each SNV, get the names of clones that it belongs to.
gt$loc <- paste(gt$chrom, gt$coord, sep=":")
locs <- unique(gt$loc)
n_snvs <- length(locs)
exon_snv$clone_name <- rep("", n_snvs)
for (i in 1:n_snvs) {
    temp <- subset(gt, loc == locs[i])
    exon_snv$clone_name[i] <- paste(temp$clone_id[temp$is_present==1], collapse = "_")
}

exon_snv <- exon_snv[,-which(names(exon_snv) == "loc")]

# The tree is as follows:
# (((A,B), (C,D)), ((E,F), ((G,H),I))).
# We will remove any annotation that is not consistent with the tree.
library(ape)
tree <- ape::read.tree(text = "(((A,B), (C,D)), ((E,F), ((G,H),I)));")
plot(tree)
valid_clusters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                    "A_B", "C_D", "E_F", "G_H",
                    "G_H_I", "A_B_C_D", "E_F_G_H_I",
                    "A_B_C_D_E_F_G_H_I")
exon_snv_valid_clones <- exon_snv[exon_snv$clone_name %in% valid_clusters,]

write.table(x = exon_snv, file = "/Users/seonghwanjun/data/cell-line/phylo/ov2295_exon_ground_truth.txt", row.names = F, col.names = T, quote = F)

