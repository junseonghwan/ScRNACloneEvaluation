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

#config_file <- "/Users/seonghwanjun/bulk-sc/scripts/Rscripts/ConfigLocal.txt"

config_file <- args[1]
config <- read.table(config_file, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])

ssm_file <- paste(OUTPUT_PATH, "ov2295_clone_snvs_exon.csv", sep="/")
ssm_final_outfile <- paste(OUTPUT_PATH, "final_ssm.txt", sep="/")
sc_final_outfile <- paste(OUTPUT_PATH, "final_sc.txt", sep="/")
sc_hp_final_outfile <- paste(OUTPUT_PATH, "final_sc_hp.txt", sep="/")

# Get all mutations
snvs <- read.table(ssm_file, header = T)
# Extract read counts from single cells. This function will write the data to file.
sc <- CombineSingleCellReads(SC_READS_PATH)
n_cells <- length(unique(sc$Cell))
length(unique(sc$ID))
length(snvs$ID)

# Ensure at least MIN_CELLS are found with variant.
temp <- sc %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d-a > 0))
temp2 <- temp[temp$n_b >= 1,]
dim(temp2)

snv_final <- subset(snvs, ID %in% temp2$ID)
dim(snv_final)
write.table(snv_final, ssm_final_outfile, row.names = F, col.names = T, quote = F, sep= "\t")

# Cell coverage by site?
summary(temp2$n_b/n_cells)
temp2[which.max(temp2$n_b/n_cells),]
temp2[(temp2$n_b > 0) & (temp2$n_d - temp2$n_b > 0),]

# How many SNVs are covered on average by each cell?
temp <- subset(sc, ID %in% temp2$ID) %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
n_snvs <- dim(snv_final)[1]
summary(temp$n_d/n_snvs)
summary(temp$n_b/n_snvs)
temp[which.max(temp$n_b/n_snvs),]

# Take a subset of the sc data for the final SNVs.
sc_final <- subset(sc, ID %in% snv_final$ID)
length(unique(sc_final$ID))
write.table(sc_final, sc_final_outfile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# Now, we estimate the hyperparameters for the single cell reads.
hyper_params.df <- data.frame(ID = snv_final$ID, alpha = rep(1, n_snvs), beta = rep(1, n_snvs), delta0 = rep(0.5, n_snvs))
sc_final$b <- sc_final$d - sc_final$a
temp <- subset(sc_final, d > 0)
for (i in 1:n_snvs) {
    id <- snv_final$ID[i]
    temp2 <- subset(temp, ID == id)
    if (dim(temp2)[1] > 0) {
        hyper_params.df[i,2:4] <- EstimateHyperparameters(temp2$b, temp2$d)
    }
}
write.table(hyper_params.df, sc_hp_final_outfile, row.names = F, col.names = T, quote = F, sep="\t")
