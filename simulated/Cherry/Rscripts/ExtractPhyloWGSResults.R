args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

PWGS_PATH <- as.character(args[1])
#PWGS_PATH <- "/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/phylowgs/"
print(PWGS_PATH)

library(rjson)

# PhyloWGS.
pwgs_mcmc_file <- paste(PWGS_PATH, "mcmc_samples.txt", sep="/")
pwgs_results_path <- paste(PWGS_PATH, "results", sep="/")
pwgs_tree_path <- paste(PWGS_PATH, "results/trees", sep="/")

mcmc_samples <- read.table(pwgs_mcmc_file, header=T, sep="\t")
mcmc_samples <- mcmc_samples[mcmc_samples$Iteration >= 0,]
best_idx <- which.max(mcmc_samples$LLH)

json<-fromJSON(file=paste(pwgs_tree_path, "/", best_idx, ".json", sep=""))
n_clusters<-length(json$mut_assignments)
phylowgs_clusters<-data.frame()
for (i in 1:n_clusters) {
    temp.df <- data.frame(V1=json$mut_assignments[[i]]$ssms, stringsAsFactors = FALSE)
    temp.df$V2 <- i
    phylowgs_clusters <- rbind(phylowgs_clusters, temp.df)
}

phylowgs_clusters <- phylowgs_clusters[order(nchar(phylowgs_clusters$V1), phylowgs_clusters$V1),]
pgws_output_file <- "/Users/seonghwanjun/data/cell-line/bulk/multi-region/reps/rep1/pwgs/predicted.txt"
names(phylowgs_clusters) <- c("ID", "Cluster")

summ.json<-fromJSON(file=paste(pwgs_results_path, "results.summ.json", sep="/"))
tree <- summ.json$trees[[best_idx]]$structure

node_count <- length(tree) + 1
parental_matrix <- matrix(0, node_count, node_count)
for (i in 1:length(tree)) {
    parental_matrix[i,tree[[i]]+1] <- 1
}

phylowgs_clusters$Cluster <- phylowgs_clusters$Cluster + 1
clusters <- sort(unique(phylowgs_clusters$Cluster))

# Compute ancestral matrix.
snv_count <- dim(phylowgs_clusters)[1]
A <- matrix(0, snv_count, snv_count)
for (i in 1:snv_count) {
    clone_i <- phylowgs_clusters$Cluster[i]
    for (j in 1:snv_count) {
        clone_j <- phylowgs_clusters$Cluster[j]
        A[i,j] <- parental_matrix[clone_i,clone_j]
    }
}

# We expect all antries of A to be 0.
error <- mean(A)
print(paste("Error:", error))

# Write the ancestral matrix.
pwgs_ancestral_matrix_path <- paste(PWGS_PATH, "ancestral_matrix.txt", sep="/")
write.table(A, pwgs_ancestral_matrix_path, quote=F, row.names = F, col.names = F)
