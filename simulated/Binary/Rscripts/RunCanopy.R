args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
REP_PATH <- args[2]
K_BEGIN <- as.numeric(args[3])
K_END <- as.numeric(args[4])
SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
OUTPUT_PATH <- paste(REP_PATH, "canopy", sep="/")

#SEED <- 157
#REP_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/"
#K_BEGIN <- 4
#K_END <- 12
#OUTPUT_PATH <- "/Users/seonghwanjun/data/simulation/binary/case1/sim0/rep0/canopy"

library(Canopy)

set.seed(SEED)

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH)
}
setwd(OUTPUT_PATH)

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
projname <- "canopy"

sampchain <- canopy.sample.nocna(R, X, K, numchains, projectname = projname,
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
pdf.name <- paste(OUTPUT_PATH, 'highest_likelihood.pdf', sep='/')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), Cluster=clone_name)
write.table(predicted, paste(OUTPUT_PATH, "predicted.csv", sep="/"), sep=",", col.names=T, row.names=F, quote=F)

snv_count <- dim(output.tree$Z)[1]
A <- matrix(0, snv_count, snv_count)
for (i in 1:snv_count) {
    clone_i <- sort(which(output.tree$Z[i,] == 1))
    for (j in 1:snv_count) {
        clone_j <- sort(which(output.tree$Z[j,] == 1))
        if (!setequal(clone_i, clone_j)) {
            A[i,j] <- as.numeric(all(clone_j %in% clone_i))
        }
    }
}

write.table(A, paste(OUTPUT_PATH, "ancestral_matrix.csv", sep="/"), row.names = F, quote=F, col.names = F)
