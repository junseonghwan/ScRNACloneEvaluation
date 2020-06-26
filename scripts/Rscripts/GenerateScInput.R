args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space with following fields.")
    print("OUTPUT_PATH: Path containing bulk data for each of the four workflows.")
    print("sc_reads_path: Path to processed single cell reads.")
    stop("Terminating...")
}

library(dplyr)
library(matrixStats)
library(ScRNAClone)
library(TailRank)

#rm(list=ls())
#config_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"

config_file <- args[1]
config <- read.table(config_file, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
EXON_SNV_PATH <- values[keys == "EXON_SNV_PATH"]
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])

sc_final_outfile <- paste(OUTPUT_PATH, "sc.txt", sep="/")
sc_hp_final_outfile <- paste(OUTPUT_PATH, "sc_hp.txt", sep="/")

# Get all mutations
snv <- read.table(EXON_SNV_PATH, header = T)
# Extract read counts from single cells. This function will write the data to file.
sc <- CombineSingleCellReads(SC_READS_PATH)
n_cells <- length(unique(sc$Cell))
n_cells
length(unique(sc$ID))
length(snvs$ID)
dim(sc)

sc <- subset(sc, ID %in% snvs$ID)
sc$ID <- factor(sc$ID, snvs$ID)

# Cell coverage by site?
#sc_final <- subset(sc, ID %in% snv_final$ID)
temp <- sc %>% group_by(ID) %>% dplyr::summarise(n_d = sum(d > 0), n_b = sum(d-a > 0))
temp2 <- temp[temp$n_b >= MIN_CELLS,]
dim(temp2)
summary(temp2$n_b/n_cells)
temp2[which.max(temp2$n_b/n_cells),]
temp2[(temp2$n_b > 0) & (temp2$n_d - temp2$n_b > 0),]
temp2

# Now, we estimate the hyperparameters for the single cell reads.
n_snvs <- dim(snv)[1]
hyper_params.df <- data.frame(ID = snv$ID, alpha = rep(1, n_snvs), beta = rep(1, n_snvs), delta0 = rep(0.5, n_snvs))
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
for (i in 1:n_snvs) {
    id <- snv_final$ID[i]
    temp2 <- subset(temp, ID == id)
    if (dim(temp2)[1] > 0) {
        hyper_params.df[i,2:4] <- EstimateHyperparameters(temp2$b, temp2$d)
    }
}

sc <- sc[,-which(names(sc) == "b")]
write.table(sc, sc_final_outfile, row.names = F, col.names = T, quote = F, sep="\t")
write.table(hyper_params.df, sc_hp_final_outfile, row.names = F, col.names = T, quote = F, sep="\t")
