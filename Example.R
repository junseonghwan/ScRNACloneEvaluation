rm(list=ls())
library(devtools)
library(dplyr)
library(sabre)
library(Rcpp)
library(VAFclusterEM)
#devtools::install_github("junseonghwan/VAFclusterEM", force=TRUE)
#install.packages("/Users/seonghwanjun/VAFclusterEM/", repos=NULL, type="source")

# Cluster the data points by chain.
ssm <- read.table("/Users/seonghwanjun/data/simulation/binary/case4/sim0/rep0/genotype_ssm.txt", header=T)
vafs <- ssm$b/ssm$d
BScite_output <- "/Users/seonghwanjun/data/simulation/binary/case4/sim0/rep0/bscite/bscite.matrices"
cluster_labels <- VAFclusterEM::GetClones(vafs, BScite_output)
truth <- read.table("/Users/seonghwanjun/data/simulation/binary/case4/sim0/rep0/cluster_labels.txt", sep=",", header=F)
vmeasure(x = truth$V2, y = cluster_labels, B = 1)
