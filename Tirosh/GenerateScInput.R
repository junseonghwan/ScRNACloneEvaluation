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

config_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/Tirosh/ConfigLocal.txt"
config_file <- args[1]

config <- read.table(config_file, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
SNV_PATH <- values[keys == "SNV_PATH"]
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])

ssm_final_outfile <- paste(OUTPUT_PATH, "/final_ssm.txt", sep="")
sc_raw_outfile <- paste(OUTPUT_PATH, "/final_sc_raw.txt", sep="")
sc_final_outfile <- paste(OUTPUT_PATH, "/final_sc.txt", sep="")
sc_hp_final_outfile <- paste(OUTPUT_PATH, "/final_sc_hp.txt", sep="")

# Get all mutations
snvs <- read.table(SNV_PATH, header = T)
#snv_raw <- read.table(paste(OUTPUT_PATH, "exon_raw_ssm.txt", sep="/"), header = T)
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
sc_min_cells <- temp[temp$n_b >= MIN_CELLS,]

sc_final <- subset(sc, ID %in% sc_min_cells$ID)
snv_final <- subset(snvs, ID %in% sc_min_cells$ID)

write.table(snv_final, ssm_final_outfile, row.names = F, col.names = T, quote = F, sep= "\t")
write.table(sc_final, sc_final_outfile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# How many SNVs are covered on average by each cell?
temp <- sc_final %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
n_snvs <- dim(snv_final)[1]
summary(temp$n_d/n_snvs)
summary(temp$n_b/n_snvs)

temp <- sc_final %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
summary(temp$n_b)

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

# Retrieve region specific SNVs.
sample_names <- as.character(sc_final$SampleName)
regions <- sapply(strsplit(sample_names, "-"), "[[", 2)
snv_by_regions <- list()
region_ids <- c(1, 3, 4, 5)

for (region in region_ids) {
    region_name <- paste("p", region, sep="")
    idx <- (regions == region_name)
    sc_region <- sc_final[idx,]
    reads_by_region <- sc_region %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
    snv_by_regions[[region]] <- reads_by_region[which(reads_by_region$n_b > 0),]
}

region_specific_snvs <- list()
for (region in region_ids) {
    region_snvs <- as.character(snv_by_regions[[region]]$ID)
    other_regions <- region_ids[region_ids != region]
    other_snv_ids <- c()
    for (other_region in other_regions) {
        other_snv_ids <- c(other_snv_ids, as.character(snv_by_regions[[other_region]]$ID))
    }
    other_snv_ids <- unique(other_snv_ids)
    region_specific_snvs[[region]] <- region_snvs[!(region_snvs %in% other_snv_ids)]
}

hist(sc_final[which(sc_final$b > 0),"b"], breaks = 50)

sc_min_cells[order(sc_min_cells$n_b, decreasing = T),]
subset(sc_final, ID == "s2879")

temp <- sc_final %>% group_by(ID) %>% dplyr::summarise(n_b = sum(d-a >= 5))
temp2 <- temp[temp$n_b >= 1,]

temp2[order(nchar(as.character(temp2$ID)), as.character(temp2$ID)),]$ID

