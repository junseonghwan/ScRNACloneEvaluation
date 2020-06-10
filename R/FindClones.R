ExtractChains <- function(file_path, mutation_count) {
    if (!file.exists(file_path)) {
        stop(paste(file_path, "does not exist."))
    }
    ret_str <- ReadParentVector(file_path, mutation_count)
    parent_vec_str <- strsplit(strsplit(ret_str, "\t")[[1]][2], " ")[[1]]
    parent_vec <- as.numeric(parent_vec_str)
    chains <- GetChains(parent_vec)
}

GetClones <- function(vafs,
                      file_path,
                      seed = 1,
                      coverage = 1000,
                      minK = 1,
                      iteration_count = 100) {
    vafs <- vafs + 1e-6 # In case vaf = 0.
    mutation_count <- length(vafs)
    chains <- ExtractChains(file_path, mutation_count)
    df <- data.frame(vaf=vafs, chain=chains)
    cluster_labels <- rep(0, mutation_count)
    cluster_id <- 0
    for (chain_ in unique(chains)) {
        idxs <- which(df$chain == chain_)
        datas <- subset(df, chain == chain_)$vaf
        if (length(datas) > 1) {
            maxK <- min(10,length(datas))
            AICsearch <- bestAICsearch(dataVec = datas, minK = minK, maxK = maxK, coverage = coverage, startseed = seed, nIterations = iteration_count, breakOnIncrease=TRUE, verbose=FALSE)
            AICs <- rep(0, length(AICsearch))
            for (i in 1:length(AICsearch)) {
                AICs[i] <- AICsearch[[i]]$AIC
            }
            cluster_labels[idxs] <- cluster_id + AICsearch[[which.min(AICs)]]$newclustermembership
            cluster_id <- cluster_id + max(AICsearch[[which.min(AICs)]]$newclustermembership)
        } else {
            cluster_labels[idxs] <- cluster_id + 1
            cluster_id <- cluster_id + 1
        }
        #print("Cluster ids assigned: ")
        #print(unique(sort(cluster_labels[idxs],decreasing = F)))
    }
    cluster_labels
}

