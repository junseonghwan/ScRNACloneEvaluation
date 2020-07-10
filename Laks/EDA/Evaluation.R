library(rjson)
library(sabre)

trimmed_gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm_gt.txt"
trimmed_pwgs_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/PWGS/trimmed_gt.txt"
trimmed_gt <- read.table(trimmed_gt_file, header=T)
trimmed_pwgs_gt <- read.table(trimmed_pwgs_file, header=T)

# PhyloWGS.
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

phylowgs_clusters <- phylowgs_clusters[order(nchar(phylowgs_clusters$V1), phylowgs_clusters$V1),]
phylowgs_clusters <- phylowgs_clusters[phylowgs_clusters$V1 %in% as.character(trimmed_pwgs_gt$id),]
mean(phylowgs_clusters$V1 == trimmed_pwgs_gt$id)
cbind(phylowgs_clusters$V1, as.character(trimmed_pwgs_gt$id))
vmeasure(x=as.numeric(trimmed_pwgs_gt$clone_name), y=phylowgs_clusters$V2, B = 1) # 0.04??? What happened?

# Without single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/Combined/wo_cells_trimmed/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% trimmed_gt$ID,]
mean(as.character(trimmed_gt$ID) == as.character(pred$V1))
vmeasure(x = trimmed_gt$clone_name, y = pred$V2) # ~0.36

# With single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/Combined//w_cells_trimmed/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% trimmed_gt$ID,]
mean(as.character(trimmed_gt$ID) == as.character(pred$V1))
vmeasure(x = trimmed_gt$clone_name, y = pred$V2) # ~0.64, nearly 2-fold increase.

# ddClone.
ddclone_path <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/ddClone/results.txt"
if (file.exists(ddclone_path)) {
    ddclone <- read.table(ddclone_path, header=T, as.is = TRUE)
    ddclone <- ddclone[ddclone$mutID %in% trimmed_gt$ID,]
    mean(as.character(ddclone$mutID) == as.character(trimmed_gt$ID))
    length(unique(ddclone$clusterID))
    vmeasure(ddclone$clusterID, trimmed_gt$clone_name, B=1) # 0.49 -- 76 clones detected.
}
