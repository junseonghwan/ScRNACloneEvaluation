library(ggplot2)
library(dplyr)
library(rjson)
source("R/CellAssign.R")

dat <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/genotype_ssm.txt", header=T, sep="\t")
sc <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/simul_sc.txt", header=T, sep="\t")
sc_hp <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/simul_sc_hp.txt", header=T, sep="\t")

snv_count <- dim(dat)[1]

# Plot cell mutation profile.
sc$b <- sc$d - sc$a

cells <- unique(sc$Cell)
ids <- dat$ID
sc$Cell <- factor(sc$Cell, levels = cells)
sc$ID <- factor(sc$ID, levels = ids)

base_size <- 11
p <- ggplot(sc, aes(Cell, ID, fill = b)) + geom_tile(colour = "white")
p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
p <- p + ylab("Loci") + xlab("Cells")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(legend.position = "none", axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/cherry_sc.pdf", p, width = 6, height=6, units="in")

cluster_labels <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/genotype/joint/tree0/cluster_labels.tsv", header=F, sep="\t")
names(cluster_labels) <- c("ID", "Cluster")
cluster_labels$Cluster <- as.character(cluster_labels$Cluster)

datum2node <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/genotype/joint/tree0/datum2node.tsv", header=F, sep="\t")

# Assign cell to nodes.
dropout_hp <- list(alpha=0.01, beta=1)
bursty_hp <- list(alpha=1, beta=0.01)
cell_assignment_nodes <- AssignCells(sc, datum2node, dropout_hp, bursty_hp, sc_hp[,2:3])
cell_assignment_truth <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/cell2node.txt", header=F, sep="\t")
names(cell_assignment_truth) <- c("Cell", "Node")

cell_assignment <- data.frame(Cell=cells, Node=cell_assignment_nodes)
cell_assignment$Node <- as.character(cell_assignment$Node)

ids <- cluster_labels[order(cluster_labels$Cluster),"ID"]
cells <- cell_assignment[order(nchar(cell_assignment$Node), cell_assignment$Node),"Cell"]
sc$Cell <- factor(sc$Cell, levels = cells)
sc$ID <- factor(sc$ID, levels = ids)

p <- ggplot(sc, aes(Cell, ID, fill = b)) + geom_tile(colour = "white")
p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
p <- p + ylab("Loci") + xlab("Cells")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(legend.position = "none", axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/cherry_inferred.pdf", p, width = 6, height=6, units="in")

# Generate posterior credible interval plot for the cellular prevalence.
# Mean squared error plot for the cellular prevalences + ancestral distance metric.
cell_prevs_truth <- read.csv("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/cellular_prev.csv", header=F, sep="\t")
A_truth <- as.matrix(read.csv("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/ancestral_matrix.csv", header=F))

path_to_states <- "/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/genotype/states/"
mcmc_iter <- 800
burn_in <- 300
thinning <- 10
begin <- burn_in / thinning + 1
end <- (mcmc_iter + burn_in) / thinning
tree_nos <- begin:end
cell_prevs.df <- data.frame()
A_pred <- matrix(0, snv_count, snv_count)
errs <- rep(0, length(tree_nos))
for (j in 1:length(tree_nos)) {
    tree_path <- paste(path_to_states, "tree", tree_nos[j], sep="")
    cell_prevs <- read.csv(paste(tree_path, "cellular_prev.tsv", sep="/"), sep="\t", header=F)
    cell_prevs$sample <- j
    cell_prevs.df <- rbind(cell_prevs.df, cell_prevs)

    A_j <- as.matrix(read.csv(paste(tree_path, "ancestral_matrix.csv", sep="/"), sep=",", header=F))
    A_pred <- A_pred + A_j
    errs[j] <- mean(abs(A_j - A_truth))
}
names(cell_prevs.df) <- c("ID", "CellFraction", "Sample")
dat <- cell_prevs.df %>% group_by(ID) %>% summarise(mean = mean(CellFraction), std = sd(CellFraction))
names(cell_prevs_truth) <- c("ID", "CellFraction")
dat <- dplyr::full_join(dat, cell_prevs_truth, "ID")
dat$ID <- factor(dat$ID, levels = ids)
p <- ggplot(dat, aes(ID, CellFraction)) + geom_point(col='red') + geom_errorbar(aes(ymin = mean - 2*std, ymax = mean + 2*std), col='black')
p <- p + xlab("Loci") + ylab("Cell Fraction Posterior Credible Interval") + theme_bw()
p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/cherry_CI.pdf", p, height=7, units = "in")

# Ancestral matrix plot: not very informative
A_pred <- A_pred/length(tree_nos)
rownames(A_pred) <- cluster_labels$ID
colnames(A_pred) <- cluster_labels$ID
A_pred.df <- reshape2::melt(A_pred)
names(A_pred.df)
A_pred.df$Var1 <- factor(A_pred.df$Var1, levels = ids)
A_pred.df$Var2 <- factor(A_pred.df$Var2, levels = ids)
p <- ggplot(A_pred.df, aes(Var1, Var2, fill=value)) + geom_tile(color="white")
p <- p + theme_bw()
p <- p + scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5)

# Retrieve the ancestral errors computed from PhyloWGS samples and compare using boxplot.
pwgs_err <- read.table("/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/phylowgs/ancestral_error.txt", header=F)
pwgs_err.df <- data.frame(Method="PhyloWGS", Error=pwgs_err$V1)
err.df <- data.frame(Method="Our Method", Error=errs)
anc_err.df <- rbind(err.df, pwgs_err.df)
p <- ggplot(anc_err.df) + geom_boxplot(aes(x = Method, y=Error, fill = Method)) + theme_bw()
p <- p + theme(legend.position = "none")
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/cherry_structure_inference.pdf", p)

# Get the best tree from PhyloWGS.
PWGS_PATH <- "/Users/seonghwanjun/data/simulation/cherry/sim0/rep0/phylowgs/"
pwgs_results_path <- paste(PWGS_PATH, "results", sep="/")
summ.json <- fromJSON(file=paste(pwgs_results_path, "results.summ.json", sep="/"))
tree_count <- length(summ.json$trees)
llh <- rep(0, tree_count)
for (i in 1:tree_count) {
    llh[i] <- summ.json$trees[[i]]$llh
}
best_idx <- which.max(llh)
# Best tree structure and cell prevalences.
summ.json$trees[[best_idx + 1]]$structure
summ.json$trees[[best_idx + 1]]$populations
