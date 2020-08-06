args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
#SEED <- 157
REP_PATH <- args[2]
#REP_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/"
K_BEGIN <- as.numeric(args[3])
#K_BEGIN <- 4
K_END <- as.numeric(args[4])
#K_END <- 6
SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
OUTPUT_PATH <- paste(REP_PATH, "canopy", sep="/")
#OUTPUT_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/canopy"
setwd(OUTPUT_PATH)

library(Canopy)

set.seed(SEED)

if (dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH)
}

bulk <- read.table(SNV_PATH, header=T, sep="\t")
X <- as.matrix(bulk$d - bulk$b)
R <- as.matrix(bulk$b)
rownames(X) <- bulk$ID
rownames(R) <- bulk$ID
colnames(X) <- "Sample1"
colnames(R) <- "Sample1"

K <- K_BEGIN:K_END
numchains <- 3
min_iter <- 10000
max_iter <- 100000
#projname <- "canopy"

sampchain <- canopy.sample.nocna(R, X, K, numchains,
                                 max.simrun = max_iter, min.simrun = min_iter,
                                 writeskip = 200, cell.line=TRUE, plot.likelihood=F)

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
write.table(predicted, paste(OUTPUT_PATH, "predicted.csv", sep="/"), sep=",", col.names=T, row.names=F, quote=F)

