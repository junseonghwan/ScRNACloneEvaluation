library(cardelino)
library(Canopy)

# Run Canopy to get configuration matrix.
# We will do a pilot run without any copy number calls.
ssm <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_wo_cnv.txt", header=T)
n_snvs <- dim(ssm)[1]

R <- matrix(ssm$b, ncol = 1)
X <- matrix(ssm$d, ncol = 1)
#rownames(R) <- paste("s", 1:n_snvs, sep="")
#rownames(X) <- paste("s", 1:n_snvs, sep="")
rownames(R) <- ssm$ID
rownames(X) <- ssm$ID
colnames(R) <- c("sample1")
colnames(X) <- c("sample1")
K <- 2:10
numchain <- 5
min_iter <- 10000
max_iter <- 100000
projname <- "cell-line"
# Slow, should be parallelized and run on server.
sampchain <- canopy.sample.nocna(R, X, K, numchain,
                                max.simrun = max_iter, min.simrun = min_iter,
                                writeskip = 200, projectname = projname, cell.line=TRUE, plot.likelihood=T)

burnin <- 10
thin <- 5
bic <- canopy.BIC(sampchain = sampchain, projectname = projname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK <- K[which.max(bic)]
post <- canopy.post(sampchain = sampchain, projectname = projname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin,
                   optK = optK, post.config.cutoff = 0.05)

samptreethin <- post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik <- post[[2]]   # likelihoods of trees in samptree
config <- post[[3]]
config.summary <- post[[4]]
print(config.summary)

config.i <- config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree <- canopy.output(post, config.i, C=NULL)
pdf.name <- paste(projname, '_config_highest_likelihood.pdf', sep='')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), CloneID=clone_name)
write.table(predicted, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/Canopy/predicted.csv", sep=",", col.names=T, row.names=F, quote=F)

# Now, we are ready to run Cardelino.
# First, get the reads from single cell data.
sc_reads <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/sc.txt", header=T)
head(sc_reads)
n_cells <- length(unique(sc_reads$Cell))
n_snvs <- length(unique(sc_reads$ID))
# Prepare a matrix of variant reads and total reads for each cell.
# Matrix dimension is N x C, where N is the number of mutations and C is the number of cells.
sc_reads$b <- sc_reads$d - sc_reads$a
sc_reads$ID <- factor(sc_reads$ID, levels = rownames(output.tree$Z))
B <- reshape2::dcast(sc_reads, formula = ID ~ Cell, value.var = "b")
D <- reshape2::dcast(sc_reads, formula = ID ~ Cell, value.var = "d")
rownames(B) <- as.character(B[,1])
rownames(D) <- as.character(D[,1])
B.mat <- as.matrix(B[,-1])
D.mat <- as.matrix(D[,-1])

# Run:
assignments <- clone_id(B.mat, D.mat, Config = output.tree$Z)
names(assignments)
prob_heatmap(assignments$prob)

df <- assign_cells_to_clones(assignments$prob)
table(df$clone)

# Use Cardelino's functionality for inferring the tree?
lsf.str("package:cardelino")
assignments_relax <- clone_id(B.mat, D.mat, Config = output.tree$Z, relax_Config = TRUE)
# Not sure if Cardelino modifies the tree in any way other than to add/eliminate an edge?

