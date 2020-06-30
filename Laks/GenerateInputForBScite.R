args = commandArgs(trailingOnly=TRUE)

print(args)

SNV_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm.txt"
SC_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_sc.txt"
OUTPUT_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/B-SCITE"
SC_VAR_READ_THRESHOLD <- 3

library(dplyr)
library(matrixStats)
library(reshape2)
library(ScRNAClone)
library(TailRank)

# B-SCITE expects the bulk data
bulk <- read.table(SNV_PATH, header=T)
sc <- read.table(SC_PATH, header=T)

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}

# Bulk data preparation.
bulk$a <- bulk$d - bulk$b
bulk_bscite <- bulk[,c("ID", "CHR", "POS", "b", "a")]
names(bulk_bscite) <- c("ID", "Chromosome", "Position", "MutantCount", "ReferenceCount")
bulk_bscite$INFO <- "NA"
bulk_bscite_file <- paste(OUTPUT_PATH, "/trimmed_bscite.bulk", sep="")
write.table(bulk_bscite, file=bulk_bscite_file, sep="\t", quote=F, row.names = F, col.names = T)

# Single cell data preparation.
if (dim(sc)[1] > 0) {
    sc$b <- sc$d - sc$a
    sc$ID <- factor(sc$ID, levels = bulk$ID)

    num_cells <- length(unique(sc$Cell))
    X.df <- dcast(sc, ID ~ Cell, value.var = "b")
    X <- as.matrix(X.df[,-1])
    mean(X.df[!is.na(X.df[,2]),2] == subset(sc, Cell == "c0")$b)
    X[is.na(X)] <- 0
    N.df <- dcast(sc, ID ~ Cell, value.var = "d")
    N <- as.matrix(N.df[,-1])
    N[is.na(N)] <- 0

    # Default value is set to 3=NA due to nature of scRNA-seq.
    Z <- matrix(3, nrow = dim(X)[1], ncol = dim(X)[2])
    # Generate another SC data with default value 0=Absent.
    Z0 <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])

    Z[X > SC_VAR_READ_THRESHOLD] <- 1
    # Uncomment if the default value is set to 0.
    Z0[(N == 0) & (Z != 1)] <- 3
    bscite_output_name <- paste(OUTPUT_PATH, "/trimmed_bscite.SC", sep="")
    bscite_output_name0 <- paste(OUTPUT_PATH, "/trimmed_bscite0.SC", sep="")
    write.table(Z, bscite_output_name, row.names = F, col.names = F, quote = F)
    write.table(Z0, bscite_output_name0, row.names = F, col.names = F, quote = F)
} else {
    Z <- matrix(3, nrow = dim(bulk)[1], ncol = 1)
    bscite_output_name <- paste(OUTPUT_PATH, "/trimmed_bscite.SC", sep="")
    write.table(Z, bscite_output_name, row.names = F, col.names = F, quote = F)
}
