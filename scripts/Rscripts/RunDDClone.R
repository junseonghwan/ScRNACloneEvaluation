args = commandArgs(trailingOnly=TRUE)
print(args)
data_path <- args[1]
mcmc_iter <- as.numeric(args[2])
print(mcmc_iter)
#data_path <- "/Users/seonghwanjun/data/simulation/binary/case4/sim0/rep0"

library(ddclone)
library(dplyr)
library(reshape2)

SC_READ_THRESHOLD <- 5

# Load the bulk data.
bulk <- read.table(paste(data_path, "genotype_ssm.txt", sep="/"), header=T)
sc <- read.table(paste(data_path, "simul_sc.txt", sep="/"), header=T)

bulkDat <- data.frame("mutation_id" = bulk$ID,
                      "ref_counts" = bulk$d - bulk$b,
                      "var_counts" = bulk$b,
                      "normal_cn" = 2,
                      "minor_cn" = bulk$minor_cn,
                      "major_cn" = bulk$major_cn)

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

output_path <- paste(data_path, "ddClone", sep="/")
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = genDat, outputPath = output_path, nameTag = '')
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = output_path, tumourContent = 1.0,
                      numOfIterations = mcmc_iter, thinning = 10, burnIn = 0,
                      seed = 1)

# Output the results.
df <- ddCloneRes$df
output_file <- paste(output_path, "results.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)
