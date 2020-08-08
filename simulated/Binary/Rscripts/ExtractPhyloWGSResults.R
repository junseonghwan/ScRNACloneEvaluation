args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

PWGS_PATH <- as.character(args[1])
#PWGS_PATH <- "/Users/seonghwanjun/data/simulation/binary/case0/sim0/rep1/phylowgs/"
print(PWGS_PATH)

library(rjson)

# PhyloWGS.
pwgs_mcmc_file <- paste(PWGS_PATH, "mcmc_samples.txt", sep="/")
pwgs_results_path <- paste(PWGS_PATH, "results", sep="/")
pwgs_tree_path <- paste(PWGS_PATH, "results/trees", sep="/")

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

ProcessSample <- function(pwgs_tree_path, idx, tree) {
    json<-fromJSON(file=paste(pwgs_tree_path, "/", idx, ".json", sep=""))
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
    A <- ComputeAncestralMatrix(phylowgs_clusters, parental_matrix)
    return(list("ancestral_matrix"=A, "clustering"=phylowgs_clusters))
}

mcmc_samples <- read.table(pwgs_mcmc_file, header=T, sep="\t")
mcmc_samples <- mcmc_samples[mcmc_samples$Iteration >= 0,]
summ.json <- fromJSON(file=paste(pwgs_results_path, "results.summ.json", sep="/"))
best_idx <- which.max(mcmc_samples$LLH)

# Write the ancestral matrix.
tree <- summ.json$trees[[best_idx+1]]$structure
ret <- ProcessSample(pwgs_tree_path, best_idx, tree)
A_best <- ret$ancestral_matrix
pwgs_clustering <-ret$clustering
pwgs_ancestral_matrix_path <- paste(PWGS_PATH, "ancestral_matrix.txt", sep="/")
write.table(A_best, pwgs_ancestral_matrix_path, quote=F, row.names = F, col.names = F)
pwgs_clustering_path <- paste(PWGS_PATH, "clustering.txt", sep="/")
write.table(pwgs_clustering, pwgs_clustering_path, quote=F, row.names = F, col.names = F)

# Select 100 samples.
sample_count <- 100
step_size <- length(mcmc_samples$Iteration) / sample_count
iterations <- mcmc_samples$Iteration + 1
sample_idxs <- iterations[seq(0, length(iterations), step_size)] - 1
# Compute the error in the ancestral matrix.
errors <- rep(0, sample_count)
for (i in 1:length(sample_idxs)) {
    idx <- sample_idxs[i]
    tree <- summ.json$trees[[idx+1]]$structure
    A_idx <- ProcessSample(pwgs_tree_path, idx, tree)
    errors[i] <- mean(A_idx)
}

write.table(x = errors, file = paste(PWGS_PATH, "ancestral_error.txt", sep=""), quote=F, row.names = F, col.names = F)


