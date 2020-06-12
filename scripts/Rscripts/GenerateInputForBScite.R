args = commandArgs(trailingOnly=TRUE)

print(args)
rep_path <- args[1]
#rep_path <- "/Users/seonghwanjun/data/simulation/binary/case0/sim0/rep0"

library(dplyr)
library(matrixStats)
library(reshape2)
library(ScRNAClone)
library(TailRank)

SC_VAR_READ_THRESHOLD <- 3

# B-SCITE expects the bulk data
bulk <- read.table(paste(rep_path, "/genotype_ssm.txt", sep=""), header=T)
sc <- read.table(paste(rep_path, "/simul_sc.txt", sep=""), header=T)

# Bulk data preparation.
bulk$a <- bulk$d - bulk$b
bulk_bscite <- bulk[,c("ID", "CHR", "POS", "b", "a")]
names(bulk_bscite) <- c("ID", "Chromosome", "Position", "MutantCount", "ReferenceCount")
bulk_bscite$INFO <- "NA"
bulk_bscite_file <- paste(rep_path, "/simul_bscite.bulk", sep="")
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

    Z <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])

    Z[X > SC_VAR_READ_THRESHOLD] <- 1
    Z[(N == 0) & (Z != 1)] <- 3
    bscite_output_name <- paste(rep_path, "/simul_bscite.SC", sep="")
    write.table(Z, bscite_output_name, row.names = F, col.names = F, quote = F)
} else {
    Z <- matrix(3, nrow = dim(bulk)[1], ncol = 1)
    bscite_output_name <- paste(rep_path, "/simul_bscite.SC", sep="")
    write.table(Z, bscite_output_name, row.names = F, col.names = F, quote = F)
}
