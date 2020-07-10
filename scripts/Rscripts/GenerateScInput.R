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
#CONFIG_FILE <- "/home/x_seoju/ScRNACloneEvaluation/Laks/Config.txt"
#CONFIG_FILE <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"

CONFIG_FILE <- args[1]
config <- read.table(CONFIG_FILE, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
SNV_PATH <- values[keys == "SNV_PATH"]
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])
MIN_READS <- as.numeric(values[keys == "MIN_READS"])

ssm_outfile <- paste(OUTPUT_PATH, "ssm.txt", sep="/")
sc_outfile <- paste(OUTPUT_PATH, "sc.txt", sep="/")
sc_hp_outfile <- paste(OUTPUT_PATH, "sc_hp.txt", sep="/")

# Get all mutations
snv <- read.table(SNV_PATH, header = T)
# Extract read counts from single cells. This function will write the data to file.
sc <- CombineSingleCellReads(SC_READS_PATH)
# Output the raw file first.
write.table(sc, paste(OUTPUT_PATH, "sc_raw.txt", sep="/"), quote=F, row.names = F, col.names = T)
n_cells <- length(unique(sc$Cell))
n_cells
length(unique(sc$ID))
length(snv$ID)

sc$ID <- factor(sc$ID, levels = snv$ID)

sc_min_reads <- subset(sc, d >= MIN_READS)
head(sc_min_reads)

# Trim the data based on minimum number of cell coverage.
cell_coverage_by_site <- sc_min_reads %>% group_by(ID) %>% summarise(n_b = sum(d - a > 0), n_d = sum(d > 0))
sites_with_min_cells <- cell_coverage_by_site[cell_coverage_by_site$n_b >= MIN_CELLS,]
sc <- subset(sc_min_reads, ID %in% sites_with_min_cells$ID)
head(sc)
length(unique(sc$Cell))

# Remove cells that do not contribute any meaningful signal.
temp <- sc %>% group_by(Cell) %>% summarise(n_b = sum(d - a > 0))
summary(temp$n_b)
cells_to_keep <- temp[temp$n_b >= 2,]$Cell
sc <- subset(sc, Cell %in% cells_to_keep)
dim(sc)
length(unique(sc$Cell))

snv <- subset(snv, ID %in% unique(sc$ID))

# Now, we estimate the hyperparameters for the single cell reads.
n_snvs <- dim(snv)[1]
hyper_params.df <- data.frame(ID = snv$ID, alpha = rep(1, n_snvs), beta = rep(1, n_snvs), delta0 = rep(0.5, n_snvs))
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
for (i in 1:n_snvs) {
    id <- snv$ID[i]
    temp2 <- subset(temp, as.character(ID) == as.character(id))
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

