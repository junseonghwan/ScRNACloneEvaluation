args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
#SEED <- 157
REP_PATH <- args[2]
#REP_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/"
MCMC_ITER <- as.numeric(args[3])
#MCMC_ITER <- 2000
THRESHOLD <- as.numeric(args[4])
#THRESHOLD <- 1

OUTPUT_PATH <- paste(REP_PATH, "ddClone", sep="/")
#OUTPUT_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/ddClone"
setwd(OUTPUT_PATH)

library(Canopy)

canopy_path <- "/Users/seonghwanjun/data/cell-line/bulk/multi-region/reps/rep1/Canopy"
if (dir.exists(canopy_path)) {
    dir.create(canopy_path)
}

R <- read.table("/Users/seonghwanjun/data/cell-line/bulk/multi-region/reps/rep1/R.txt", header = T, as.is = T)
X <- read.table("/Users/seonghwanjun/data/cell-line/bulk/multi-region/reps/rep1/X.txt", header = T, as.is = T)
R <- as.matrix(R)
X <- as.matrix(X)
rownames(R)
colnames(R)

K <- 3:10
numchain <- 3
min_iter <- 10000
max_iter <- 100000
projname <- "cell-line-rep1"

sampchain <- canopy.sample.nocna(R, X, K, numchain,
                                 max.simrun = max_iter, min.simrun = min_iter,
                                 writeskip = 200, projectname = projname, cell.line=TRUE, plot.likelihood=F)

burnin <- 10
thin <- 5
bic <- canopy.BIC(sampchain = sampchain, projectname = projname, K = K,
                  numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK <- K[which.max(bic)]
post <- canopy.post(sampchain = sampchain, projectname = projname, K = K,
                    numchain = numchain, burnin = burnin, thin = thin,
                    optK = optK, post.config.cutoff = 0.01)

samptreethin <- post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik <- post[[2]]   # likelihoods of trees in samptree
config <- post[[3]]
config.summary <- post[[4]]
print(config.summary)

config.i <- config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree <- canopy.output(post, config.i, C=NULL)
pdf.name <- paste(canopy_path, 'highest_likelihood.pdf', sep='/')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), Cluster=clone_name)
write.table(predicted, "/Users/seonghwanjun/data/cell-line/bulk/multi-region/reps/rep1/Canopy/predicted.csv", sep=",", col.names=T, row.names=F, quote=F)
