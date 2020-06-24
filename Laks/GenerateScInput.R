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
library(plyr)
library(ScRNAClone)
library(TailRank)

#config_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"

config_file <- args[1]
config <- read.table(config_file, header=F, as.is = TRUE)
print(config)
keys <- as.character(config$V1)
values <- as.character(config$V2)
ssm_file <- values[keys == "SNV_PATH"]
OUTPUT_PATH <- values[keys == "OUTPUT_PATH"]
OUTPUT_PREFIX <- values[keys == "OUTPUT_PREFIX"]
SC_READS_PATH <- values[keys == "SC_READS_PATH"]
MIN_CELLS <- as.numeric(values[keys == "MIN_CELLS"])

ssm_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "ssm.txt", sep="")
sc_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "sc.txt", sep="")
sc_hp_final_outfile <- paste(OUTPUT_PATH, "/", OUTPUT_PREFIX, "sc_hp.txt", sep="")

# Get all mutations
snvs <- read.table(ssm_file, header = T)
# Extract read counts from single cells. This function will write the data to file.
sc <- CombineSingleCellReads(SC_READS_PATH)
n_cells <- length(unique(sc$Cell))
length(unique(sc$ID))
length(snvs$ID)

# Cell coverage by site?
#sc_final <- subset(sc, ID %in% snv_final$ID)
temp <- sc %>% group_by(ID) %>% dplyr::summarise(n_d = sum(d > 0), n_b = sum(d-a > 0))
temp2 <- temp[temp$n_b >= MIN_CELLS,]
dim(temp2)
summary(temp2$n_b/n_cells)
temp2[which.max(temp2$n_b/n_cells),]
temp2[(temp2$n_b > 0) & (temp2$n_d - temp2$n_b > 0),]
temp2$ID

# Get the final SNV set.
snv_final <- subset(snvs, ID %in% temp2$ID)
head(snv_final)

# Because our data is written in a multi-region format, we will have to unwind and merge it into one.
ret <- laply(strsplit(as.character(snv_final$b), ","), as.numeric)
snv_final$b <- rowSums(ret)
ret <- laply(strsplit(as.character(snv_final$d), ","), as.numeric)
snv_final$d <- rowSums(ret)
ret <- laply(strsplit(as.character(snv_final$MajorCN), ","), as.numeric)
snv_final$MajorCN <- round(rowMeans(ret))
ret <- laply(strsplit(as.character(snv_final$MinorCN), ","), as.numeric)
snv_final$MinorCN <- round(rowMeans(ret))
head(snv_final)

loc_col_idx <- which(names(snv_final) == "loc")
write.table(snv_final[,-loc_col_idx], ssm_final_outfile, row.names = F, col.names = T, quote = F, sep= "\t")

sum(rowSums(ret > 0) == 1) # 25 region specific mutations.
sum(rowSums(ret > 1) == 1) # Make that 30 region specific mutations. Variant read of 1 is likely due to sequencing error.
snv_final[which(rowSums(ret > 1) == 1),] # These are the sites.
snv_final[which(rowSums(ret > 1) == 2),] # These are the mutations on two samples.
snv_final[which(rowSums(ret > 1) == 3),] # These are the mutations on all three samples.

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

#sc_final %>% group_by(ID) %>% summarise(n_b = sum(b >= 1))
temp <- subset(sc_final, ID == "s11")
temp[(temp$b > 0),]

cell <- subset(sc_final, Cell == "c254")
cell[which(cell$b > 0),]

cell_clone_membership <- read.table("/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_clusters.csv", sep=",", header=T)

clone_ids <- sort(unique(cell_clone_membership$clone_id))
cell_prev <- rep(0, length(clone_ids))
for (i in 1:length(clone_ids)) {
    clone_id <- clone_ids[i]
    cell_prev[i] <- mean(cell_clone_membership$clone_id == clone_id)
}

library(sabre)
predicted <- read.table("/Users/seonghwanjun/data/cell-line/phylo/joint/tree0/cluster_labels.txt", header=F, sep=",")
truth <- read.table("/Users/seonghwanjun/data/cell-line/phylo/ov2295_exon_ground_truth_curated.txt", header=T)
vmeasure(truth$clone_name, predicted$V2, B = 1)
cbind(truth$clone_name, predicted$V2)
