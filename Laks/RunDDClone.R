args = commandArgs(trailingOnly=TRUE)
print(args)
SNV_PATH <- args[1]
SC_PATH <- args[2]
MCMC_ITER <- as.numeric(args[3])
SEED <- as.numeric(args[4])
SC_READ_THRESHOLD <- as.numeric(args[5])

SNV_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm.txt"
SC_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_sc.txt"
OUTPUT_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/ddClone/"
MCMC_ITER <- 2000
SEED <- 157
SC_READ_THRESHOLD <- 3

library(ddclone)
library(dplyr)
library(reshape2)

# Load the bulk data.
bulk <- read.table(SNV_PATH, header=T)
sc <- read.table(SC_PATH, header=T)

bulkDat <- data.frame("mutation_id" = bulk$ID,
                      "ref_counts" = bulk$d - bulk$b,
                      "var_counts" = bulk$b,
                      "normal_cn" = 2,
                      "minor_cn" = bulk$MajorCN,
                      "major_cn" = bulk$MinorCN)

if (dim(sc)[1] > 0) {
    sc$ID <- factor(sc$ID, levels = bulk$ID)
    sc$b <- sc$d - sc$a
    sc$z <- as.numeric(sc$b >= SC_READ_THRESHOLD)
    sum(sc$z)
    genDat <- dcast(sc, formula = Cell ~ ID, value.var = "z")
    genDat <- genDat[,-1]
    genDat[is.na(genDat)] <- 0
} else {
    n_snvs <- dim(bulk)[1]
    genDat <- as.data.frame(matrix(0, nrow = 1, ncol = n_snvs))
    names(genDat) <- bulk$ID
}

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = genDat, outputPath = OUTPUT_PATH, nameTag = '')
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
                      numOfIterations = mcmc_iter, thinning = 10, burnIn = 0,
                      seed = seed)

# Output the results.
df <- ddCloneRes$df
output_file <- paste(OUTPUT_PATH, "results.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)
