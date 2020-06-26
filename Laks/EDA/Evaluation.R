library(rjson)
library(sabre)

gt <- read.table("/Users/seonghwanjun/data/cell-line/bulk/A90554A/ground_truth.txt", header=T)
truth <- as.numeric(gt$clone_name)

# Without single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/A90554A/wo_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt$ID,]
mean(as.character(gt$ID) == as.character(pred$V1))
unique(truth)
unique(pred$V2)
vmeasure(x = truth, y = pred$V2) # Terrible performance.

# PhyloWGS.
pgs_mcmc_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/mcmc_samples.txt"
pgs_tree_path <- "/Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/witness/data/cell-line-results/"
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

phylowgs_clusters <- phylowgs_clusters[order(nchar(phylowgs_clusters$V1), phylowgs_clusters$V1),]
phylowgs_clusters <- phylowgs_clusters[phylowgs_clusters$V1 %in% gt$ID,]
mean(phylowgs_clusters$V1 == gt$ID)
vmeasure(x=truth, y=phylowgs_clusters$V2, B = 1)

# With single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/A90554A/w_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt$ID,]
mean(as.character(gt$ID) == as.character(pred$V1))
unique(truth)
unique(pred$V2)
vmeasure(x = truth, y = pred$V2) # Terrible performance.
cbind(truth, pred$V2)
