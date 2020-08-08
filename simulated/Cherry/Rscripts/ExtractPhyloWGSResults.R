args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

PWGS_PATH <- as.character(args[1])
#PWGS_PATH <- "/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/phylowgs/"
print(PWGS_PATH)

CHAIN_COUNT <- 4

library(rjson)

# PhyloWGS.
pwgs_results_path <- paste(PWGS_PATH, "results", sep="/")
pwgs_tree_path <- paste(pwgs_results_path, "trees", sep="/")
A_truth <- as.matrix(read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/ancestral_matrix.csv", sep=","))

# Compute ancestral matrix.
ComputeAncestralMatrix <- function(phylowgs_clusters, parental_matrix) {
    snv_count <- dim(phylowgs_clusters)[1]
    A <- matrix(0, snv_count, snv_count)
    for (i in 1:snv_count) {
        clone_i <- phylowgs_clusters$Cluster[i]
        for (j in 1:snv_count) {
            clone_j <- phylowgs_clusters$Cluster[j]
            A[i,j] <- parental_matrix[clone_i,clone_j]
        }
    }
    return(A)
}

ProcessSample <- function(pwgs_tree_path, tree) {
    json<-fromJSON(file=pwgs_tree_path)
    n_clusters<-length(json$mut_assignments)
    phylowgs_clusters<-data.frame()
    for (i in 1:n_clusters) {
        temp.df <- data.frame(V1=json$mut_assignments[[i]]$ssms, stringsAsFactors = FALSE)
        temp.df$V2 <- i
        phylowgs_clusters <- rbind(phylowgs_clusters, temp.df)
    }

    phylowgs_clusters <- phylowgs_clusters[order(nchar(phylowgs_clusters$V1), phylowgs_clusters$V1),]
    names(phylowgs_clusters) <- c("ID", "Cluster")

    node_count <- length(unique(phylowgs_clusters$Cluster)) + 1
    parental_matrix <- matrix(0, node_count, node_count)
    for (i in 1:length(tree)) {
        parental_matrix[i,tree[[i]]+1] <- 1
    }

    #clusters <- sort(unique(phylowgs_clusters$Cluster))
    return (ComputeAncestralMatrix(phylowgs_clusters, parental_matrix))
}

summ.json <- fromJSON(file=paste(pwgs_results_path, "results.summ.json", sep="/"))

# Select 100 samples.
pwgs_trees <- list.files(pwgs_tree_path, pattern = "*.json")
pwgs_trees <- pwgs_trees[order(nchar(pwgs_trees), pwgs_trees)]
tree_count <- length(pwgs_trees)

sample_idxs <- seq(1, tree_count, 10)
sample_count <- length(sample_idxs)
# Compute the error in the ancestral matrix.
errors <- rep(0, sample_count)
for (i in 1:length(sample_idxs)) {
    idx <- sample_idxs[i]
    tree <- summ.json$trees[[idx]]$structure
    A_idx <- ProcessSample(paste(pwgs_tree_path, pwgs_trees[idx], sep="/"), tree)
    errors[i] <- mean(abs(A_idx - A_truth))
}

write.table(x = errors, file = paste(PWGS_PATH, "ancestral_error.txt", sep=""), quote=F, row.names = F, col.names = F)

llh <- rep(0, tree_count)
for (i in 1:tree_count) {
    llh[i] <- summ.json$trees[[i]]$llh
}
best_idx <- which.max(llh)
summ.json$trees[[best_idx+1]]$populations
best_tree <- summ.json$trees[[best_idx+1]]$structure
best_A <- ProcessSample(paste(pwgs_tree_path, pwgs_trees[best_idx], sep="/"), tree)
mean(abs(best_A - A_truth))
write.table(best_A, "/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/phylowgs/ancestral_matrix.txt", row.names = F, quote = F, col.names = F)
