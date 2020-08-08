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
df0 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/case0.csv", header=T)
df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/case1.csv", header=T)
df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/case2.csv", header=T)
df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/case3.csv", header=T)
df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/case4.csv", header=T)

# Without single cells: PhyloWGS and Canopy. Canopy performed worse, so we will just use PhyloWGS.
pwgs <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/phylowgs.csv", header=T)
canopy <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/canopy.csv", header=T)

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
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_vmeasure.pdf", p)

p <- ggplot(df, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Mean Absolute Error")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_ancestral.pdf", p)

ProcessCells <- function(sc, output_path) {
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
    p
    ggsave(paste(output_path, "binary_sc.pdf", sep="/"), p, width = 12, height=4, units="in")
}

# Plot single cell mutation profile.
dropout_hp <- list(alpha=0.01, beta=1)
bursty_hp <- list(alpha=1, beta=0.01)
for (case in 1:4) {
    case_path <- paste("/Users/seonghwanjun/data/simulation/binary/case", case, sep="")
    for (rep in 0:9) {
        rep_path <- paste(case_path, "/sim0/rep", rep, sep="")
        dat <- read.table(paste(rep_path, "genotype_ssm.txt", sep="/"), header=T, sep="\t")
        sc <- read.table(paste(rep_path, "simul_sc.txt", sep="/"), header=T, sep="\t")
        datum2node <- read.table(paste(rep_path, "genotype/joint/tree0/datum2node.tsv", sep="/"), header=F, sep="\t")

        # Assign cell to nodes.
        cell_assignment <- AssignCells(sc, datum2node, dropout_hp, bursty_hp, sc_hp[,2:3])
        write.table(cell_assignment, file = paste(rep_path, "genotype/joint/tree0/cell_assignment.txt", sep="/"), quote=F, row.names = F, col.names = F)
    }
}

# Figure 2: Comparison to existing joint analysis approaches.

# First, generate process B-SCITE results.
for (case in 1:4) {
    case_path <- paste("/Users/seonghwanjun/data/simulation/binary/case", case, sep="")
    sim_path <- paste(case_path, "sim0", sep="/")
    GetBSciteResults(sim_path, rep_end = 10)
}

# Then, go run ComparisonToOtherMethods from Jupyter notebook and run evaluation, output as dataframe.
# This is annoying but somehow R does not have a reliable library for computing clustering metrics.

# Generate the boxplot.
ddclone_df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/ddClone_case1.csv", header=T)
ddclone_df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/ddClone_case2.csv", header=T)
ddclone_df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/ddClone_case3.csv", header=T)
ddclone_df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/ddClone_case4.csv", header=T)

ddclone_df1$CellCount <- 100
ddclone_df2$CellCount <- 200
ddclone_df3$CellCount <- 400
ddclone_df4$CellCount <- 800

ddclone_df1$Method <- "ddClone"
ddclone_df2$Method <- "ddClone"
ddclone_df3$Method <- "ddClone"
ddclone_df4$Method <- "ddClone"

bscite_df1 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/bscite_case1.csv", header=T)
bscite_df2 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/bscite_case2.csv", header=T)
bscite_df3 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/bscite_case3.csv", header=T)
bscite_df4 <- read.csv("/Users/seonghwanjun/ScRNACloneEvaluation/data/simulated/binary/bscite_case4.csv", header=T)

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
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_vmeasure.pdf", p)

p <- ggplot(df, aes(as.factor(CellCount), AdjMutualInformation, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Adjusted Mutual Information")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_mutual_information.pdf", p)

df_a <- rbind(our_df, bscite_df)
p <- ggplot(df_a, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
p <- p + theme_bw() + xlab("Cell Count") + ylab("Mean Absolute Error")
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/binary_comparison_ancestral.pdf", p)


# TODO: Plot cellular prevalence C.I. and error.
