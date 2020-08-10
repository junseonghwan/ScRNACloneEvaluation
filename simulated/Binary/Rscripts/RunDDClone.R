args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
REP_PATH <- args[2]
MCMC_ITER <- as.numeric(args[3])
THRESHOLD <- as.numeric(args[4])
#SEED <- 157
#REP_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/"
#MCMC_ITER <- 5
#THRESHOLD <- 3

OUTPUT_PATH <- paste(REP_PATH, "ddClone", sep="/")
#OUTPUT_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/ddClone"

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}
setwd(OUTPUT_PATH)

library(ddclone)
library(dplyr)
library(reshape2)

# Load the bulk data.
SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
SC_PATH <- paste(REP_PATH, "simul_sc.txt", sep="/")

bulk <- read.table(SNV_PATH, header=T)

bulkDat <- data.frame("mutation_id" = bulk$ID,
                      "ref_counts" = bulk$d - bulk$b,
                      "var_counts" = bulk$b,
                      "normal_cn" = 2,
                      "minor_cn" = bulk$major_cn,
                      "major_cn" = bulk$minor_cn)

sc <- read.table(SC_PATH, header=T)
sc$b <- sc$d - sc$a
sc$ID <- factor(sc$ID, levels = bulk$ID)
cells <- unique(sc$Cell)
sc$Cell <- factor(sc$Cell, levels = cells)

sc_var <- dcast(sc, Cell ~ ID, value.var = "b")
sc_depth <- dcast(sc, Cell ~ ID, value.var = "d")
sc_var[is.na(sc_var)] <- 0
sc_depth[is.na(sc_depth)] <- 0

sc_var <- as.matrix(sc_var[,-1])
sc_depth <- as.matrix(sc_depth[,-1])

# Call variant using threshold.
sc_mut_matrix <- matrix(0, ncol = length(bulk$ID), nrow = length(cells))
sc_mut_matrix[sc_var >= THRESHOLD] <- 1
rownames(sc_mut_matrix) <- as.character(cells)
colnames(sc_mut_matrix) <- as.character(bulk$ID)

ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = sc_mut_matrix, outputPath = OUTPUT_PATH, nameTag = '')

start <- proc.time()
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
                      numOfIterations = MCMC_ITER, thinning = 10, burnIn = 0,
                      seed = SEED)
end <- proc.time()
diff <- end - start
elapsed_time_seconds <- diff[3]

# Output the results.
df <- ddCloneRes$df
output_file <- paste(OUTPUT_PATH, "results.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)

timing_file <- paste(OUTPUT_PATH, "timing.txt", sep="/")
write.table(elapsed_time_seconds, file = timing_file, row.names = F, col.names = F, quote=F)
