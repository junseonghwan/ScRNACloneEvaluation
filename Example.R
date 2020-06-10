rm(list=ls())
library(devtools)
library(dplyr)
library(sabre)
library(Rcpp)
library(VAFclusterEM)
#devtools::install_github("junseonghwan/VAFclusterEM", force=TRUE)
#remove.packages("VAFclusterEM2", lib = "/Library/Frameworks/R.framework/Resources/library")
#load_all()

# Cluster the data points by chain.
ssm <- read.table("/Users/seonghwanjun/data/binary_cn/case5/sim0/rep19/genotype_ssm.txt", header=T)
vafs <- ssm$b/ssm$d
BScite_output <- "/Users/seonghwanjun/data/binary_cn/case5/sim0/rep19/bscite/bscite.matrices"
cluster_labels <- VAFclusterEM::GetClones(vafs, BScite_output)
VAFclusterEM::ExtractChains
truth <- read.table("/Users/seonghwanjun/data/binary_cn/case5/sim0/rep19/cluster_labels.txt", sep=",", header=F)
vmeasure(x = truth$V2, y = cluster_labels, B = 1)
