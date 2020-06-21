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

config_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"

config_file <- args[1]
config <- read.table(config_file, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
SNV_PATH <- values[keys == "SNV_PATH"]
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
OUTPUT_PREFIX <- values[keys == "OUTPUT_PREFIX"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])

ssm_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "ssm.txt", sep="")
sc_raw_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "sc_raw.txt", sep="")
sc_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "sc.txt", sep="")
sc_hp_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "sc_hp.txt", sep="")

# Get all mutations
snvs <- read.table(SNV_PATH, header = T)
# Extract read counts from single cells.
if (!file.exists(sc_raw_outfile)) {
    sc <- CombineSingleCellReads(SC_READS_PATH)
    write.table(sc, file = sc_raw_outfile, row.names = F, col.names = T, quote = F, sep= "\t" )
} else {
    sc <- read.table(sc_raw_outfile, sep="\t", header=T)
}
dim(sc)
# Save the read counts.
n_cells <- length(unique(sc$Cell))
length(unique(sc$ID))
length(snvs$ID)

# Ensure at least MIN_CELLS are found with variant.
temp <- sc %>% group_by(ID) %>% dplyr::summarise(n_d = sum(d > 0), n_b = sum(d-a > 0))
temp2 <- temp[temp$n_b >= MIN_CELLS,]
dim(temp2)

snv_final <- subset(snvs, ID %in% temp2$ID)
dim(snv_final)
write.table(snv_final, ssm_final_outfile, row.names = F, col.names = T, quote = F, sep= "\t")

# Cell coverage by site?
summary(temp2$n_b/n_cells)
temp2[which.max(temp2$n_b/n_cells),]
idx <- which(temp2$n_b/temp2$n_d >= 0.9)
mut_ids <- temp2[idx,]$ID
head(snvs[snvs$ID %in% mut_ids,], 10)
subset(sc, ID == id)
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
