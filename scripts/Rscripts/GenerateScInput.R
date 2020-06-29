args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

print(args)

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

ssm_outfile <- paste(OUTPUT_PATH, "ssm.txt", sep="/")
sc_outfile <- paste(OUTPUT_PATH, "sc.txt", sep="/")
sc_hp_outfile <- paste(OUTPUT_PATH, "sc_hp.txt", sep="/")

ssm_trimmed_outfile <- paste(OUTPUT_PATH, "trimmed_ssm.txt", sep="/")
sc_trimmed_outfile <- paste(OUTPUT_PATH, "trimmed_sc.txt", sep="/")
sc_hp_trimmed_outfile <- paste(OUTPUT_PATH, "trimmed_sc_hp.txt", sep="/")

# Get all mutations
snv <- read.table(EXON_SNV_PATH, header = T)
# Extract read counts from single cells. This function will write the data to file.
sc <- CombineSingleCellReads(SC_READS_PATH)
n_cells <- length(unique(sc$Cell))
n_cells
length(unique(sc$ID))
length(snv$ID)

sc$ID <- factor(sc$ID, levels = snv$ID)

# Now, we estimate the hyperparameters for the single cell reads.
n_snvs <- dim(snv)[1]
hyper_params.df <- data.frame(ID = snv$ID, alpha = rep(1, n_snvs), beta = rep(1, n_snvs), delta0 = rep(0.5, n_snvs))
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
for (i in 1:n_snvs) {
    id <- snv$ID[i]
    temp2 <- subset(temp, ID == id)
    if (dim(temp2)[1] > 0) {
        hyper_params.df[i,2:4] <- EstimateHyperparameters(temp2$b, temp2$d)
    }
}

sc <- sc[,-which(names(sc) == "b")]
if (sum(names(snv) == "loc") > 0) {
    snv <- snv[,-which(names(snv) == "loc")]
}
write.table(snv, ssm_outfile, row.names = F, col.names = T, quote = F, sep="\t")
write.table(sc, sc_outfile, row.names = F, col.names = T, quote = F, sep="\t")
write.table(hyper_params.df, sc_hp_outfile, row.names = F, col.names = T, quote = F, sep="\t")

# Trim the data based on minimum number of cells.
cell_coverage_by_site <- sc %>% group_by(ID) %>% summarise(n_b = sum(d - a > 0))
sites_with_min_cells <- cell_coverage_by_site[cell_coverage_by_site$n_b >= MIN_CELLS,]
trimmed_snv <- subset(snv, ID %in% sites_with_min_cells$ID)
trimmed_sc <- subset(sc, ID %in% sites_with_min_cells$ID)
trimmed_sc_hp <- subset(hyper_params.df, ID %in% sites_with_min_cells$ID)
write.table(trimmed_snv, ssm_trimmed_outfile, row.names = F, col.names = T, quote = F, sep="\t")
write.table(trimmed_sc, sc_trimmed_outfile, row.names = F, col.names = T, quote = F, sep="\t")
write.table(trimmed_sc_hp, sc_hp_trimmed_outfile, row.names = F, col.names = T, quote = F, sep="\t")
