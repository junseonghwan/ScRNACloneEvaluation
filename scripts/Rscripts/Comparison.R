rm(list=ls())

library(devtools)
library(dplyr)
library(Rcpp)
library(sabre)
sourceCpp("src/rcpp_hello_world.cpp")
source("R/FindClones.R")

data_path <- "/Users/seonghwanjun/data/simulation/binary/case4/sim0/rep0"
true_clusters <- read.table(paste(data_path, "cluster_labels.txt", sep="/"), header=F, sep=",")
cell_prevs <- read.table(paste(data_path, "cellular_prev.csv", sep="/"), header=F, sep=",")

length(unique(true_clusters$V2))

# ddClone
ddclone <- read.table(paste(data_path, "ddClone", "results.txt", sep="/"), header=T)
vmeasure(true_clusters$V2, ddclone$clusterID)
length(unique(ddclone$clusterID))
mean(abs(ddclone$phi - cell_prevs$V2))

# B-SCITE
ssms <- read.table(paste(data_path, "genotype_ssm.txt", sep="/"), header=T)
vafs <- ssms$b/ssms$d
mutation_count <- length(vafs)
bscite_output <- paste(data_path, "bscite", "bscite.matrices", sep="/")
bscite_clones <- GetClones(vafs, bscite_output)
length(unique(bscite_clones))
vmeasure(true_clusters$V2, bscite_clones)

# Our method
ours <- read.table(paste(data_path, "genotype", "joint", "tree0", "cluster_labels.txt", sep="/"), header=F, sep=",")
ours_cp <- read.table(paste(data_path, "genotype", "joint", "tree0", "cellular_prev.csv", sep="/"), header=F, sep=",")
vmeasure(true_clusters$V2, ours$V2)
length(unique(ours$V2))
mean(abs(ours_cp$V2 - cell_prevs$V2))
