library(dplyr)
library(rjson)
library(Rcpp)
library(sabre)
library(ggplot2)
library(pdfCluster)
library(VAFclusterEM)
source("R/FindClones.R")
sourceCpp("src/rcpp_hello_world.cpp")

#gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_gt_curated.txt"
gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_gt.txt"
gt <- read.table(gt_file, header=T)

# Take subset of gt to be valid set.
valid_clusters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                    "A_B", "C_D", "E_F", "G_H",
                    "G_H_I", "A_B_C_D", "E_F_G_H_I",
                    "A_B_C_D_E_F_G_H_I")
gt_valid <- gt[(gt$clone_name %in% valid_clusters),]
gt[!(gt$clone_name %in% valid_clusters),]

# ddClone.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/ddClone/results.txt", header=T)
pred <- pred[pred$mutID %in% gt_valid$ID,]
mean(as.character(gt_valid$ID) == as.character(pred$mutID))
vmeasure(x = as.numeric(gt_valid$clone_name), y = as.numeric(pred$clusterID))
pdfCluster::adj.rand.index(as.numeric(gt_valid$clone_name), pred$clusterID)
ret <- pred %>% group_by(clusterID) %>% summarise(n_ = n())
subset(pred, clusterID %in% ret[ret$n_ > 1,]$clusterID)

# B-SCITE.
valid_cluster_idx <- which(gt$clone_name %in% valid_clusters)
vafs <- gt$b/gt$d
bscite_output <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/B-SCITE/bscite.matrices"
bscite_cluster_labels <- GetClones(vafs, bscite_output)
bscite_cluster_labels_valid <- bscite_cluster_labels[valid_cluster_idx]
length(bscite_cluster_labels_valid) == length(gt_valid$ID)
sabre::vmeasure(as.numeric(gt_valid$ID), bscite_cluster_labels_valid)
pdfCluster::adj.rand.index(as.numeric(gt_valid$clone_name), cluster_labels_valid)
cbind(as.character(gt_valid$clone_name), cluster_labels_valid)

# With single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/w_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/w_cells_w_total_cn/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt_valid$ID,]
mean(as.character(gt_valid$ID) == as.character(pred$V1))
vmeasure(x = as.numeric(gt_valid$clone_name), y = pred$V2)
pdfCluster::adj.rand.index(as.numeric(gt_valid$clone_name), pred$V2)

datum2node<-read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/w_cells/joint/tree0/datum2node.txt", sep=",")
datum2node$gt <- gt$clone_name
datum2node_valid <- subset(datum2node, V1 %in% gt_valid$ID)
datum2node_valid
length(unique(pred$V2))
length(unique(gt_valid$clone_name))

sample_path <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/w_cells/states/"
burn_in <- 10
n_states <- 110 - burn_in
vmeasures <- rep(0, n_states)
for (n in 1:n_states) {
    tree_path <- paste(sample_path, "tree", (n + burn_in), "/cluster_labels.txt", sep="")
    pred <- read.table(tree_path, header=F, sep=",")
    pred <- pred[pred$V1 %in% gt_valid$ID,]
    mean(as.character(gt_valid$ID) == as.character(pred$V1))
    vmeas <- vmeasure(x = as.numeric(gt_valid$clone_name), y = pred$V2)
    vmeasures[n] <- vmeas$v_measure
}
boxplot(vmeasures)

# Without single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/wo_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt_valid$ID,]
mean(as.character(gt_valid$ID) == as.character(pred$V1))
vmeasure(x = as.numeric(gt_valid$clone_name), y = pred$V2)
pdfCluster::adj.rand.index(as.numeric(gt_valid$clone_name), pred$V2)

length(unique(pred$V2))
length(unique(gt_valid$clone_name))

sample_path <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/wo_cells/states/"
burn_in <- 10
n_states <- 110 - burn_in
vmeasures_wo_cells <- rep(0, n_states)
for (n in 1:n_states) {
    tree_path <- paste(sample_path, "tree", (n + burn_in), "/cluster_labels.txt", sep="")
    pred <- read.table(tree_path, header=F, sep=",")
    pred <- pred[pred$V1 %in% gt_valid$ID,]
    mean(as.character(gt_valid$ID) == as.character(pred$V1))
    vmeas <- vmeasure(x = as.numeric(gt_valid$clone_name), y = pred$V2)
    vmeasures_wo_cells[n] <- vmeas$v_measure
}
boxplot(vmeasures_wo_cells)

# PhyloWGS
pgs_mcmc_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/witness/data/cell-line-results/mcmc_samples.txt"
pgs_tree_path <- "/Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/witness/data/cell-line-results/cell-line.mutass/"
mcmc_samples <- read.table(pgs_mcmc_file, header=T, sep="\t")
mcmc_samples <- mcmc_samples[mcmc_samples$Iteration >= 0,]
best_idx <- which.max(mcmc_samples$LLH)

json<-fromJSON(file=paste(pgs_tree_path, best_idx, ".json", sep=""))
n_clusters<-length(json$mut_assignments)
phylowgs_clusters<-data.frame()
for (i in 1:n_clusters) {
    temp.df <- data.frame(V1=json$mut_assignments[[i]]$ssms, stringsAsFactors = FALSE)
    temp.df$V2 <- i
    phylowgs_clusters <- rbind(phylowgs_clusters, temp.df)
}

id_map <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/PWGS/id_map.txt", header=T)
phylowgs_clusters <- phylowgs_clusters[order(nchar(phylowgs_clusters$V1), phylowgs_clusters$V1),]
names(phylowgs_clusters) <- c("pwgs_id", "cluster")
predicted <- left_join(phylowgs_clusters, id_map)
predicted <- predicted[predicted$original_id %in% as.character(gt_valid$ID),]
mean(predicted$original_id == gt_valid$ID)
vmeasure(x=as.numeric(gt_valid$clone_name), y=predicted$cluster, B = 1)
adj.rand.index(as.numeric(gt_valid$clone_name), predicted$cluster)

# Canopy
predicted <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/Canopy/predicted.csv", sep=",", header=T)
predicted <- predicted[predicted$ID %in% as.character(gt_valid$ID),]
mean(predicted$ID == gt_valid$ID)
vmeasure(x=as.numeric(gt_valid$clone_name), y=as.numeric(predicted$CloneID), B = 1)
pdfCluster::adj.rand.index(as.numeric(gt_valid$clone_name), as.numeric(predicted$CloneID))
