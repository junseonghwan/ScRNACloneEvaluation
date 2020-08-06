args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

REP_PATH <- as.character(args[1])
prediction_path <- paste(REP_PATH, "genotype", "joint", "tree0", sep="/")

datum2node <- read.table(paste(prediction_path, "datum2node.tsv", sep="/"), sep="\t", header=F)
names(datum2node) <- c("ID", "Cluster")
n_snvs <- dim(datum2node)[1]

A <- matrix(0, n_snvs, n_snvs)
for (i in 1:n_snvs) {
    clone_i <- datum2node$Cluster[i]
    for (j in 1:n_snvs) {
        clone_j <- datum2node$Cluster[j]
        if (clone_i != clone_j & grepl(clone_i, clone_j, fixed = TRUE)) {
            A[i,j] <- 1
        }
    }
}

output_path <- paste(REP_PATH, "predicteed_ancestral_matrix.txt", sep="/")
write.table(A, output_path, quote=F, col.names = F, row.names = F)
