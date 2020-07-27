args = commandArgs(trailingOnly=TRUE)
print(args)
data_path <- args[1]
#data_path <- "/Users/seonghwanjun/data/simulation/binary/case3/sim0/rep0"

library(dplyr)
library(matrixStats)
library(ScRNAClone)
library(TailRank)

# Process simulation data to get sc hyperparameters.
ssm <- read.table(paste(data_path, "simul_ssm.txt", sep="/"), header=T, sep="\t")
sc <- read.table(paste(data_path, "simul_sc.txt", sep="/"), header=T)
sc_hp_outfile <- paste(data_path, "simul_sc_hp.txt", sep="/")

# Now, we estimate the hyperparameters for the single cell reads.
n_snvs <- dim(ssm)[1]
hyper_params.df <- data.frame(ID = ssm$ID, alpha = rep(1, n_snvs), beta = rep(1, n_snvs), delta0 = rep(0.5, n_snvs))
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
for (i in 1:n_snvs) {
  id <- as.character(ssm$ID[i])
  temp2 <- subset(temp, ID == id)
  if (dim(temp2)[1] > 0) {
    hyper_params.df[i,2:4] <- EstimateHyperparameters(temp2$b, temp2$d)
  }
}
write.table(hyper_params.df, sc_hp_outfile, row.names = F, col.names = T, quote = F, sep="\t")
