# Assign cells to nodes.
# 1. Read in cell reads and mutations harboured by each node.
# 2. For each cell, compute the likelihood that it should be assigned to a node.
# 3. Compute the marginalization to obtain the cell assignment probabilities.
# 4. Choose the one with maximal probability.

library(dplyr)
library(devtools)
library(matrixStats)
library(Rcpp)
library(RcppGSL)
library(reshape2)
sourceCpp("src/rcpp_hello_world.cpp")

# dropout_hp: named list containing 'alpha' and 'beta'.
# bursty_hp: named list containing 'alpha' and 'beta'.
# biallelic_hp: numeric matrix with the first column the alpha and the second column the beta
AssignCells <- function(cell_data,
                        datum2node,
                        dropout_hp,
                        bursty_hp,
                        biallelic_hp) {
    names(datum2node) <- c("ID", "Node")
    snv_count <- length(datum2node$ID)

    nodes <- as.character(unique(datum2node$Node))
    nodes_breadth_first <- nodes[order(nchar(nodes), nodes)]

    mut_ids <- unique(as.character(cell_data$ID))
    mut_ids <- mut_ids[order(nchar(mut_ids), mut_ids)]
    cells <- unique(as.character(cell_data$Cell))
    cells <- cells[order(nchar(cells), cells)]

    cell_data$ID <- factor(cell_data$ID, levels = mut_ids)
    cell_data$Cell <- factor(cell_data$Cell, levels = cells)
    if (!("b" %in% names(cell_data))) {
        cell_data$b <- cell_data$d - cell_data$a
    }
    var_reads <- dcast(cell_data, Cell ~ ID, value.var = "b")
    var_reads[is.na(var_reads)] <- 0
    total_reads <- dcast(cell_data, Cell ~ ID, value.var = "d")
    total_reads[is.na(total_reads)] <- 0

    dropout_hp_mat <- matrix(c(dropout_hp$alpha, dropout_hp$beta), nrow = snv_count, ncol=2, byrow=T)
    bursty_hp_mat <- matrix(c(bursty_hp$alpha, bursty_hp$beta), nrow = snv_count, ncol=2, byrow=T)

    log_unnorm_liks <- IdentifyCellMutationStatus(datum2node,
                                                  nodes_breadth_first,
                                                  mut_ids,
                                                  as.matrix(var_reads[,-1]),
                                                  as.matrix(total_reads[,-1]),
                                                  as.matrix(dropout_hp_mat),
                                                  as.matrix(bursty_hp_mat),
                                                  as.matrix(biallelic_hp))
    cells_with_no_reads <- (rowSums(log_unnorm_liks) == 0)
    log_unnorm_liks2 <- log_unnorm_liks[!cells_with_no_reads,]
    #dim(log_unnorm_liks2)
    log_norms <- apply(log_unnorm_liks2, 1, logSumExp)
    cell_assign_probs <- exp(log_unnorm_liks2 - log_norms)
    err <- (rowSums(cell_assign_probs) - 1) > 1e-6
    if (sum(err) > 0) {
        print(paste("Row does not sum to 1 for ", which(err)))
    }
    cell_assignment <- nodes_breadth_first[apply(cell_assign_probs, 1, which.max)]
    return(cell_assignment)
}
