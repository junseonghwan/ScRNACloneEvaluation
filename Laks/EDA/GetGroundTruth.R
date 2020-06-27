rm(list=ls())

library(dplyr)
library(GenomicRanges)
library(ScRNAClone)
library(sabre)

# There are 3 ground truth files to generate:
# 1. All exon SNVs identified.
# 2. Trimmed SNVs.
# 3. SSM's in PhyloWGS comparison study.

annotation_file <- "/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs.csv"

snv_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/ssm.txt"
trimmed_snv_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm.txt"
pwgs_snv_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/PWGS/snv.txt"

gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/ssm_gt.txt"
trimmed_gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm_gt.txt"
pwgs_gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/PWGS/gt.txt"

get_gt <- function(snv, laks_snv, outfile) {
    snv.gr <- ConstructGranges(snv$CHR, snv$POS, width = 0)
    laks_snv.gr <- ConstructGranges(laks_snv$chrom, laks_snv$coord, width = 0)

    overlaps <- findOverlaps(laks_snv.gr, snv.gr)
    gt <- laks_snv[overlaps@from,]

    # For each SNV, get the names of clones.
    gt$loc <- paste(gt$chrom, gt$coord, sep=":")
    locs <- unique(gt$loc)
    n_snvs <- length(locs)
    snv$clone_name <- rep("", n_snvs)
    for (i in 1:n_snvs) {
        temp <- subset(gt, loc == locs[i])
        snv$clone_name[i] <- paste(temp$clone_id[temp$is_present==1], collapse = "_")
    }

    # The tree is as follows:
    # (((A,B), (C,D)), (E,F)), ((G,H),I))).
    # We will remove any annotation that is not consistent with the tree.
    valid_clusters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                        "A_B", "C_D", "E_F", "G_H",
                        "G_H_I", "A_B_C_D", "E_F_G_H_I",
                        "A_B_C_D_E_F_G_H_I")
    snv_with_valid_clones <- snv[(snv$clone_name %in% valid_clusters),]
    write.table(x = snv_with_valid_clones, file = outfile, row.names = F, col.names = T, quote = F)
    return(snv_with_valid_clones)
}

snv <- read.table(snv_file, header=T)
trimmed_snv <- read.table(trimmed_snv_file, header=T)
pwgs_snv <- read.table(pwgs_snv_file, header=T)
laks_snv <- read.table(annotation_file, header=T, sep=",")

ret1 <- get_gt(snv, laks_snv, gt_file)
ret1 <- get_gt(trimmed_snv, laks_snv, trimmed_gt_file)
ret1 <- get_gt(pwgs_snv, laks_snv, pwgs_gt_file)

