library(matrixStats)
library(Rcpp)
library(TailRank)
sourceCpp("src/rcpp_hello_world.cpp")

datum2node <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/w_cells/joint/tree0/datum2node.txt", sep=",", header=F, as.is=T)
names(datum2node) <- c("ID", "Node")
nodes <- unique(datum2node$Node)
nodes <- nodes[order(nodes, nodes)]

Config <- GetConfigMatrix(datum2node, nodes)
n_clones <- dim(Config)[2]
rownames(Config) <- datum2node$ID
colnames(Config) <- paste("Clone", 1:n_clones, sep="")

sc_reads <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/sc.txt", header=T)
n_cells <- length(unique(sc_reads$Cell))
n_snvs <- length(unique(sc_reads$ID))
# Prepare a matrix of variant reads and total reads for each cell.
# Matrix dimension is N x C, where N is the number of mutations and C is the number of cells.
sc_reads$b <- sc_reads$d - sc_reads$a
sc_reads$ID <- factor(sc_reads$ID, levels = rownames(output.tree$Z))
B <- reshape2::dcast(sc_reads, formula = ID ~ Cell, value.var = "b")
D <- reshape2::dcast(sc_reads, formula = ID ~ Cell, value.var = "d")
rownames(B) <- as.character(B[,1])
rownames(D) <- as.character(D[,1])
B.mat <- as.matrix(B[,-1])
D.mat <- as.matrix(D[,-1])

colnames(B.mat)
colnames(D.mat)
assignments <- clone_id(B.mat, D.mat, Config = Config)
prob_heatmap(assignments$prob)
df <- assign_cells_to_clones(assignments$prob, threshold = 0.5)
table(df$clone)
nodes

# Assign cells to tree by computing likelihood of cell assignment to each clone.
sc_hp <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/sc_hp.txt", header=T)
alpha_hps <- sc_hp$alpha
beta_hps <- sc_hp$beta
delta0_hps <- sc_hp$delta0
delta0_hps[delta0_hps == 0] <- 1e-6
delta0_hps[delta0_hps == 1] <- 1 - 1e-6

scRNA_likelihood <- function(b, d, mut_profile,
                             alpha_hps, beta_hps, delta0_hps,
                             seq_err = 0.001, alpha0 = 0.01, beta0 = 0.01) {
    n_snvs <- length(mut_profile)

    idx_presence <- (mut_profile == 1 & !is.na(b))
    log_lik1 <- dbb(b[idx_presence], d[idx_presence], alpha_hps[idx_presence], beta_hps[idx_presence], log=TRUE) + log(1 - delta0_hps[idx_presence])
    log_lik2 <- dbb(b[idx_presence], d[idx_presence], alpha0, beta0, log=TRUE) + log(delta0_hps[idx_presence])
    log_val1 <- apply(cbind(log_lik1, log_lik2), 1, matrixStats::logSumExp)

    idx_absence <- (mut_profile == 0 & !is.na(b))
    log_val2 <- dbb(b[idx_absence], d[idx_absence], seq_err, 1 - seq_err, log=TRUE)
    log_vals <- c(log_val1, log_val2)
    if (length(log_vals) != sum(!is.na(b))) {
        stop("Error: bug in single cell likelihood calculation.")
    }
    return(sum(log_vals))
}

log_liks <- matrix(-Inf, nrow = n_cells, ncol = n_clones)
for (i in 1:n_clones) {
    for (j in 1:n_cells) {
        log_liks[j,i] <- scRNA_likelihood(B.mat[,j], D.mat[,j], Config[,i], alpha_hps, beta_hps, delta0_hps)
    }
}

normalize <- function(log_vals) {
    log_norm <- logSumExp(log_vals)
    log_vals <- log_vals - log_norm
    return(exp(log_vals))
}
dim(log_liks)
assign_probs <- t(apply(log_liks, 1, normalize))
rowSums(assign_probs)
cell_assign <- apply(assign_probs, 1, which.max)
table(cell_assign)
colnames(B.mat)[(cell_assign == 6)]
