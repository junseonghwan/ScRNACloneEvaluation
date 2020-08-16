library(ggplot2)
library(dplyr)
library(rjson)
library(Rcpp)
sourceCpp("src/rcpp_hello_world.cpp")
source("R/CellAssign.R")
source("R/EvaluationFunctions.R")
source("R/FindClones.R")

# Figure 1: clustering and ancestral metric against the baseline.

# With single cells:
df0 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/case0.csv", header=T)
df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/case1.csv", header=T)
df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/case2.csv", header=T)
df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/case3.csv", header=T)
df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/case4.csv", header=T)

# Without single cells: PhyloWGS and Canopy. Canopy performed worse, so we will just use PhyloWGS.
pwgs <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/phylowgs.csv", header=T)
#canopy <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/canopy.csv", header=T)

pwgs$CellCount <- 0
df1$CellCount <- 100
df2$CellCount <- 200
df3$CellCount <- 400
df4$CellCount <- 800

pwgs$Method <- "PhyloWGS"
df1$Method <- "OurMethod"
df2$Method <- "OurMethod"
df3$Method <- "OurMethod"
df4$Method <- "OurMethod"

df <- rbind(pwgs, df1, df2, df3, df4)

base_size <- 11

p <- ggplot(df, aes(as.factor(CellCount), VMeasure, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("V-Measure") + theme(legend.position = "none")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_vmeasure.pdf", p, height=8, units = "in")

p <- ggplot(df, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Ancestral Metric Mean Absolute Error") + theme(legend.position = "none")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_ancestral.pdf", p, height = 8, units = "in")

ProcessCells <- function(sc, output_file) {
    sc$b <- sc$d - sc$a
    p <- ggplot(sc, aes(Cell, ID, fill = b)) + geom_tile(colour = "white")
    p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
    p <- p + xlab("Cell") + ylab("Loci")
    p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
    p <- p + theme(legend.position = "none", axis.ticks = element_blank())
    p <- p + theme(axis.title.x =element_text(size = base_size * 2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2))
    p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(output_file, p, width = 12, height=4, units="in")
}


# Plot single cell mutation profile.
dropout_hp <- list(alpha=0.5, beta=1)
bursty_hp <- list(alpha=1, beta=0.01)
for (case in 1:4) {
    case_path <- paste("/Users/seonghwanjun/data/simulation/binary4/case", case, sep="")
    for (rep in 0:9) {
        rep_path <- paste(case_path, "/sim0/rep", rep, sep="")
        dat <- read.table(paste(rep_path, "genotype_ssm.txt", sep="/"), header=T, sep="\t")
        sc <- read.table(paste(rep_path, "simul_sc.txt", sep="/"), header=T, sep="\t")
        sc_hp <- read.table(paste(rep_path, "simul_sc_hp.txt", sep="/"), header=T, sep="\t")
        datum2node <- read.table(paste(rep_path, "genotype/joint/tree0/datum2node.tsv", sep="/"), header=F, sep="\t")
        names(datum2node) <- c("ID", "Node")
        ProcessCells(sc, paste(rep_path, "genotype/joint/tree0", "sc_before.pdf", sep="/"))

        # Assign cell to nodes.
        cell_assignment <- AssignCells(sc, datum2node, dropout_hp, bursty_hp, sc_hp[,2:3])
        cell.df <- data.frame(Cell = unique(sc$Cell), Node=cell_assignment)
        cell.df$Node <- as.character(cell.df$Node)
        cell_order <- cell.df[order(nchar(cell.df$Node), cell.df$Node), "Cell"]

        # Plot single cell data after clustering by mutations and cells.
        datum2node$Node <- as.character(datum2node$Node)
        datum2node.sorted <- datum2node[order(nchar(datum2node$Node), datum2node$Node),]

        sc_copy <- sc
        sc_copy$Cell <- factor(sc_copy$Cell, levels = cell_order)
        sc_copy$ID <- factor(sc_copy$ID, levels = datum2node.sorted$ID)

        sc_copy.df <- dplyr::full_join(sc_copy, cell.df)
        sc_copy.df$b <- sc_copy.df$d - sc_copy.df$a
        sc_copy.df$Cell <- as.numeric(sc_copy.df$Cell)
        sc_copy.df$ID <- as.numeric(sc_copy.df$ID)
        ProcessCells(sc_copy, paste(rep_path, "genotype/joint/tree0", "sc_after.pdf", sep="/"))

        write.table(cell.df, file = paste(rep_path, "genotype/joint/tree0/cell_assignment.txt", sep="/"), quote=F, row.names = F, col.names = F)
    }
}

# Figure 2: Comparison to existing joint analysis approaches.

# Generate the boxplot.
ddclone_df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/ddClone_case1.csv", header=T)
ddclone_df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/ddClone_case2.csv", header=T)
ddclone_df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/ddClone_case3.csv", header=T)
ddclone_df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/ddClone_case4.csv", header=T)

ddclone_df1$CellCount <- 100
ddclone_df2$CellCount <- 200
ddclone_df3$CellCount <- 400
ddclone_df4$CellCount <- 800

ddclone_df1$Method <- "ddClone"
ddclone_df2$Method <- "ddClone"
ddclone_df3$Method <- "ddClone"
ddclone_df4$Method <- "ddClone"

bscite_df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/bscite_case1.csv", header=T)
bscite_df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/bscite_case2.csv", header=T)
bscite_df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/bscite_case3.csv", header=T)
bscite_df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary4/bscite_case4.csv", header=T)

bscite_df1$CellCount <- 100
bscite_df2$CellCount <- 200
bscite_df3$CellCount <- 400
bscite_df4$CellCount <- 800

bscite_df1$Method <- "B-SCITE"
bscite_df2$Method <- "B-SCITE"
bscite_df3$Method <- "B-SCITE"
bscite_df4$Method <- "B-SCITE"

ddclone_df <- rbind(ddclone_df1, ddclone_df2, ddclone_df3, ddclone_df4)
ddclone_df$AncestralMetric <- NA
bscite_df <- rbind(bscite_df1, bscite_df2, bscite_df3, bscite_df4)
our_df <- rbind(df1, df2, df3, df4)
df <- rbind(our_df, bscite_df, ddclone_df)
names(df)

p <- ggplot(df, aes(as.factor(CellCount), VMeasure, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("V-Measure")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(legend.text = element_text(size = base_size*1.2), legend.title = element_blank())
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_vmeasure.pdf", p, height = 8, units = "in")

p <- ggplot(df, aes(as.factor(CellCount), AdjMutualInformation, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Adjusted Mutual Information")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(legend.text = element_text(size = base_size*1.2), legend.title = element_blank())
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_mutual_information.pdf", p, height = 8, units = "in")

df_a <- rbind(our_df, bscite_df)
p <- ggplot(df_a, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Ancestral Metric Mean Absolute Error")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(legend.text = element_text(size = base_size*1.2), legend.title = element_blank())
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_ancestral.pdf", p, height = 8, units = "in")


# TODO: Plot cellular prevalence C.I. and error.
