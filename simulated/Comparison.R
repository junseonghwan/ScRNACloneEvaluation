rm(list=ls())

library(devtools)
library(dplyr)
library(ggplot2)
library(Rcpp)
library(sabre)
sourceCpp("src/rcpp_hello_world.cpp")
source("R/FindClones.R")
source("R/EvaluationFunctions.R")

n_cases <- 4 # not inclusive
sim_path <- "/Users/seonghwanjun/data/simulation/binary"
expression_levels <- c(0, 10, 20, 40)
df <- data.frame()
for (case in 1:n_cases) {
    sim_case_path <- paste(sim_path, "/case", (case-1), "/sim0", sep="")
    # Check that our simulation runs have completed.
    CheckSimulationCompleted(sim_case_path)
    vmeas <- ComputeVMeasure(sim_case_path)
    cp_err <- GetCellPrevError(sim_case_path)
    df <- rbind(df, data.frame(Method="Ours", ExprLevel=expression_levels[case], VMeasure = vmeas, CellPrevErr = cp_err))

    # ddClone
    ddClone_results <- GetDDCloneResults(sim_case_path)
    df <- rbind(df, data.frame(Method="ddClone", ExprLevel=expression_levels[case], VMeasure = ddClone_results$vmeasure, CellPrevErr = ddClone_results$mean_abs_err))

    #bscite_results <- GetBSciteResults(sim_case_path)
    #df <- rbind(df, data.frame(Method="B-SCITE", ExprLevel=expression_levels[case], VMeasure = bscite_results$vmeasure, CellPrevErr = bscite_results$mean_abs_err))
}

names(df)
p <- ggplot(na.omit(df), aes(x = as.factor(ExprLevel), y = VMeasure)) + geom_boxplot() + theme_bw()
p <- p + facet_wrap(~ Method) + xlab("Expression Levels")
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/V-measure-comparison.pdf", p)

p <- ggplot(na.omit(df), aes(x = as.factor(ExprLevel), y = CellPrevErr)) + geom_boxplot() + theme_bw()
p <- p + facet_wrap(~ Method) + xlab("Expression Levels") + ylab("Mean Abs Error")
ggsave("/Users/seonghwanjun/Dropbox/Research/papers/sc-bulk-phylo/figures/mean-abs-error.pdf", p)
