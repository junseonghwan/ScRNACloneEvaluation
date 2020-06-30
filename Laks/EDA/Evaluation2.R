library(rjson)
library(sabre)

gt_file <- "/Users/seonghwanjun/data/cell-line/bulk/Combined/ssm_gt.txt"
gt <- read.table(gt_file, header=T)

# PhyloWGS.
pgs_mcmc_file <- "/Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/witness/data/cell-line-results/cell-line.mutass/mcmc_samples.txt"
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
phylowgs_clusters <- phylowgs_clusters[phylowgs_clusters$V1 %in% as.character(pwgs_gt$id),]
mean(phylowgs_clusters$V1 == pwgs_gt$id)
vmeasure(x=as.numeric(pwgs_gt$clone_name), y=phylowgs_clusters$V2, B = 1)

# Without single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/Combined/wo_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt$ID,]
mean(as.character(gt$ID) == as.character(pred$V1))
vmeasure(x = gt$clone_name, y = pred$V2) # 0.26

# With single cells.
pred <- read.table("/Users/seonghwanjun/data/cell-line/bulk/Combined//w_cells/joint/tree0/cluster_labels.txt", header=F, sep=",")
pred <- pred[pred$V1 %in% gt$ID,]
mean(as.character(gt$ID) == as.character(pred$V1))
vmeasure(x = gt$clone_name, y = pred$V2)
