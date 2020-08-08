library(dplyr)
library(ggplot2)

CheckSimulationCompleted <- function(sim_case_path, rep_end = 20,
                                     mcmc_iter = 1000, burn_in = 100, thinning = 10) {
    done <- TRUE
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        ours <- paste(rep_path, "genotype", "states", sep="/")
        n_trees <- length(list.files(path = paste(rep_path, "genotype", "states", sep="/"), pattern="tree", include.dirs = TRUE, full.names = TRUE))
        expected_n_trees <- (mcmc_iter + burn_in)/thinning + 1 # +1 is for the initial tree.
        #print(n_trees)
        #print(expected_n_trees)
        if (expected_n_trees != n_trees) {
            print(paste(rep_path, "did not complete."))
            done <- FALSE
        }
    }
    return(done)
}

ComputeVMeasure <- function(sim_case_path, rep_begin = 0, rep_end = 20)
{
    vmeasures <- rep(0, rep_end)
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        true_clusters <- read.table(paste(rep_path, "cluster_labels.txt", sep="/"), header=F, sep=",")$V2
        best_tree_path <- paste(rep_path, "genotype", "joint", "tree0", "cluster_labels.txt", sep="/")
        predicted_clustering <- read.csv(best_tree_path, header=F)$V2
        ret <- vmeasure(x = true_clusters, y = predicted_clustering, B = 1)
        vmeasures[rep_no + 1] <- ret$v_measure
    }
    return(vmeasures)
}

GetCellPrevError <- function(sim_case_path, rep_begin = 0, rep_end = 20) {
    mean_abs_err <- rep(0, rep_end)
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        true_cell_prevs <- read.table(paste(rep_path, "cellular_prev.csv", sep="/"), header=F, sep=",")
        predicted_cell_prev_path <- paste(rep_path, "genotype", "joint", "tree0", "cellular_prev.csv", sep="/")
        predicted_cell_prevs <- read.csv(predicted_cell_prev_path, header=F)
        mean_abs_err[rep_no + 1] <- mean(abs(true_cell_prevs$V2 - predicted_cell_prevs$V2))
    }
    return(mean_abs_err)
}

CellPrevCoverage <- function(sim_case_path,
                             rep_end = 20,
                             mcmc_iter = 1000,
                             burn_in = 100,
                             thinning = 10)
{
    coverage <- rep(0, rep_end)
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        cell_prevs <- read.table(paste(rep_path, "cellular_prev.csv", sep="/"), header=F, sep=",")
        names(cell_prevs) <- c("ID", "TrueCellPrev")
        n_snvs <- dim(cell_prevs)[1]
        ours <- paste(rep_path, "genotype", "joint", sep="/")
        trees <- list.files(path = paste(rep_path, "genotype", "states", sep="/"), pattern="tree", include.dirs = TRUE, full.names = TRUE)
        trees <- trees[order(nchar(trees), trees)]
        trees <- trees[-1] # remove tree0.
        n_trees <- mcmc_iter/thinning
        first_tree_idx <- burn_in/thinning + 1
        last_tree_idx <- burn_in/thinning + n_trees
        mcmc_sample_ids <- first_tree_idx:last_tree_idx
        cell_prev_posterior_samples <- matrix(0, nrow = n_snvs, ncol = n_trees)
        for (i in 1:length(mcmc_sample_ids)) {
            cell_prev_path <- paste(trees[mcmc_sample_ids[i]], "cellular_prev.csv", sep="/")
            predicted_cell_prevs <- read.csv(cell_prev_path, header=F, as.is = TRUE)
            names(predicted_cell_prevs) <- c("ID", "CellPrev")
            cell_prev_posterior_samples[,i] <- predicted_cell_prevs$CellPrev
        }
        cp_mean <- rowMeans(cell_prev_posterior_samples)
        cp_sd <- apply(cell_prev_posterior_samples, 1, sd)
        cp.df <- data.frame(ID = as.character(cell_prevs$ID), cellular_prev = cp_mean, sd = cp_sd)

        ll <- cp.df$cellular_prev - 2*cp.df$sd
        uu <- cp.df$cellular_prev + 2*cp.df$sd
        coverage[rep_no+1] <- mean(cell_prevs$TrueCellPrev >= ll & cell_prevs$TrueCellPrev <= uu)

        # Sort the SNVs based on their cluster grouping.
        cp.df$ID <- factor(cp.df$ID, levels = cell_prevs[order(cell_prevs$TrueCellPrev),"ID"])
        cp.df <- left_join(cp.df, cell_prevs, "ID")
        p <- ggplot(cp.df, aes(ID, cellular_prev))
        p <- p + geom_errorbar(aes(ymin = cellular_prev - 2*sd, ymax = cellular_prev + 2*sd), width=0.05)
        p <- p + geom_point(aes(ID, TrueCellPrev), col='red', size=0.5)
        p <- p + theme_bw() + xlab("") + ylab("")
        p <- p + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = 9, angle = 90, hjust = 0, colour = "grey50"))
        ggsave(filename = paste(rep_path, "genotype", "cell_prev.pdf", sep="/"), p)
    }
}

GetDDCloneResults<-function(sim_case_path,
                            rep_end = 20)
{
    mean_abs_err <- rep(NA, rep_end)
    vmeasures <- rep(NA, rep_end)
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        cell_prevs <- read.table(paste(rep_path, "cellular_prev.csv", sep="/"), header=F, sep=",")
        cluster_labels <- read.table(paste(rep_path, "cluster_labels.txt", sep="/"), header=F, sep=",")

        ddclone_path <- paste(rep_path, "ddClone", "results.txt", sep="/")
        if (file.exists(ddclone_path)) {
            ddclone <- read.table(ddclone_path, header=T)
            mean_abs_err[rep_no+1] <- mean(abs(ddclone$phi - cell_prevs$V2))
            ret <- vmeasure(ddclone$clusterID, cluster_labels$V2, B=1)
            vmeasures[rep_no+1] <- ret$v_measure
            mean_abs_err[rep_no+1] <- mean(abs(ddclone$phi - cell_prevs$V2))
        }
    }
    return(data.frame(vmeasure=vmeasures, mean_abs_err=mean_abs_err))
}

GetBSciteResults <- function(sim_case_path, rep_end = 20)
{
    mean_abs_err <- rep(NA, rep_end)
    ancestral_metric <- rep(NA, rep_end)
    for (rep_no in 0:(rep_end-1)) {
        rep_path <- paste(sim_case_path, "/rep", rep_no, sep="")
        ssms <- read.table(paste(rep_path, "genotype_ssm.txt", sep="/"), header=T, as.is=TRUE)
        vafs <- ssms$b/ssms$d
        mutation_count <- length(vafs)

        cell_prevs <- read.table(paste(rep_path, "cellular_prev.csv", sep="/"), header=F, sep="\t")
        cluster_labels <- read.table(paste(rep_path, "cluster_labels.txt", sep="/"), header=F, sep=",")

        bscite_output <- paste(rep_path, "bscite", "bscite.matrices", sep="/")
        if (file.exists(bscite_output)) {
            bscite_clones <- GetClones(vafs, bscite_output)
            df <- data.frame(ID=ssms$ID, Cluster=bscite_clones, VAF=vafs)
            mean_abs_err[rep_no+1] <- mean(abs(vafs - cell_prevs$V2))
            #ret <- vmeasure(cluster_labels$V2, bscite_clones)
            #vmeasures[rep_no+1] <- ret$v_measure
            write.table(data.frame(ID=ssms$ID, Cluster=bscite_clones), paste(rep_path, "/bscite/results.txt", sep=""), quote=F, row.names=F)
        }
    }
    #return(data.frame(vmeasure=vmeasures, mean_abs_err=mean_abs_err))
}
