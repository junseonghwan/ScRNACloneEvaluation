args = commandArgs(trailingOnly=TRUE)
print(args)
SNV_PATH <- args[1]
SC_MUT_MATRIX_PATH <- args[2]
MCMC_ITER <- as.numeric(args[3])
SEED <- as.numeric(args[4])

SNV_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_wo_cn.txt"
SC_MUT_MATRIX_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/mut_matrix.txt"
OUTPUT_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ddClone/"
MCMC_ITER <- 2000
SEED <- 157

library(ddclone)
library(dplyr)
library(reshape2)

# Load the bulk data.
bulk <- read.table(SNV_PATH, header=T)

bulkDat <- data.frame("mutation_id" = bulk$ID,
                      "ref_counts" = bulk$d - bulk$b,
                      "var_counts" = bulk$b,
                      "normal_cn" = 2,
                      "minor_cn" = bulk$MajorCN,
                      "major_cn" = bulk$MinorCN)

sc_mut_matrix <- read.table(SC_MUT_MATRIX_PATH, header=T)
head(sc_mut_matrix)
# Check that the mutations are correctly ordered.
mean(colnames(sc_mut_matrix) == as.character(bulkDat$mutation_id))

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = sc_mut_matrix, outputPath = OUTPUT_PATH, nameTag = '')
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
                      numOfIterations = mcmc_iter, thinning = 10, burnIn = 0,
                      seed = seed)

# Output the results.
df <- ddCloneRes$df
output_file <- paste(OUTPUT_PATH, "results.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)

# We may remove cells that have mutation at exactly one loci as these would more than likely to confuse ddClone.
sc_mut_matrix2 <- sc_mut_matrix[rowSums(sc_mut_matrix) > 1,]
dim(sc_mut_matrix2)
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = sc_mut_matrix2, outputPath = OUTPUT_PATH, nameTag = '')
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
                      numOfIterations = mcmc_iter, thinning = 10, burnIn = 0,
                      seed = seed)

# Output the results.
df <- ddCloneRes$df
output_file <- paste(OUTPUT_PATH, "results2.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)
