
get_overlap <- function(ssm_file, strelka_pass_file) {
    df1 <- read.table(ssm_file, header=T)
    df2 <- read.table(strelka_pass_file, header=T)
    df1.gr <- ConstructGranges(df1$CHR, df1$POS, 0)
    df2.gr <- ConstructGranges(df2$CHR, df2$POS, 0)
    ret <- findOverlaps(df1.gr, df2.gr)
    return(df1[ret@from,])
}

find_overlaps <- function(strelka_pass_files) {
    n_samples <- length(strelka_pass_files)
    pass_snvs <- read.table(strelka_pass_files[1], header=T)
    pass.gr <- ConstructGranges(pass_snvs$CHROM, start=pass_snvs$POS, width=0)
    if (n_samples >= 2) {
        for (n in 2:n_samples) {
            pass_snvs_n <- read.table(strelka_pass_files[n], header=T)
            pass.gr_n <- ConstructGranges(pass_snvs_n$CHROM, start=pass_snvs_n$POS, width=0)
            pass.gr <- ConstructGranges(pass_snvs$CHROM, start=pass_snvs$POS, width=0)
            ret <- findOverlaps(pass.gr, pass.gr_n)
            pass_snvs <- pass_snvs[ret@from,]
        }
    }
    return(pass_snvs)
}

strelka_pass_files <- c("/Users/seonghwanjun/data/cell-line/phylo/strelka_pass1.txt",
                        "/Users/seonghwanjun/data/cell-line/phylo/strelka_pass2.txt",
                        "/Users/seonghwanjun/data/cell-line/phylo/strelka_pass3.txt")
ssm_files <- c("/Users/seonghwanjun/data/cell-line/phylo/final_1_ssm.txt",
               "/Users/seonghwanjun/data/cell-line/phylo/final_2_ssm.txt",
               "/Users/seonghwanjun/data/cell-line/phylo/final_3_ssm.txt")
ssm_all_files <- c("/Users/seonghwanjun/data/cell-line/phylo/ssm1.txt",
                   "/Users/seonghwanjun/data/cell-line/phylo/ssm2.txt",
                   "/Users/seonghwanjun/data/cell-line/phylo/ssm3.txt")
ret1 <- get_overlap(ssm_file[1], strelka_pass_file[1])
ret2 <- get_overlap(ssm_file[2], strelka_pass_file[2])
ret3 <- get_overlap(ssm_file[3], strelka_pass_file[3])
ret1$LOC <- paste(ret1$CHR, ret1$POS, sep=":")
ret2$LOC <- paste(ret2$CHR, ret2$POS, sep=":")
ret3$LOC <- paste(ret3$CHR, ret3$POS, sep=":")
sum(ret1$LOC %in% ret2$LOC)
sum(ret1$LOC %in% ret3$LOC)
sum(ret2$LOC %in% ret3$LOC)

ret1[which(ret1$LOC %in% ret3$LOC),1:5]
ret3[which(ret3$LOC %in% ret1$LOC),1:5]
ret2[which(ret2$LOC %in% ret3$LOC),1:5]

# These are mutations found across all samples.
ancestral_snvs <- find_overlaps(strelka_pass_files)
ancestral_snvs.gr <- ConstructGranges(ancestral_snvs$CHROM, start = ancestral_snvs$POS, width = 0)

sc_raw_files <- c("/Users/seonghwanjun/data/cell-line/phylo/final_1_sc_raw.txt",
                  "/Users/seonghwanjun/data/cell-line/phylo/final_2_sc_raw.txt",
                  "/Users/seonghwanjun/data/cell-line/phylo/final_3_sc_raw.txt")
sc1_raw <- read.table(sc_raw_files[1], header=T)
sc2_raw <- read.table(sc_raw_files[2], header=T)
sc3_raw <- read.table(sc_raw_files[3], header=T)

# Find sample specific mutations.
ssm1 <- read.table(strelka_pass_files[1], header=T)
ssm1.gr <- ConstructGranges(ssm1$CHR, start = ssm1$POS, width = 0)
ret1 <- findOverlaps(ssm1.gr, ancestral_snvs.gr)
dim(ssm1)[1] - dim(ancestral_snvs)[1]

ssm2 <- read.table(strelka_pass_files[2], header=T)
ssm2.gr <- ConstructGranges(ssm2$CHR, start = ssm2$POS, width = 0)
ret2 <- findOverlaps(ssm2.gr, ancestral_snvs.gr)
dim(ssm2)[1] - dim(ancestral_snvs)[1]

ssm3 <- read.table(strelka_pass_files[3], header=T)
ssm3.gr <- ConstructGranges(ssm3$CHR, start = ssm3$POS, width = 0)
ret3 <- findOverlaps(ssm3.gr, ancestral_snvs.gr)
dim(ssm3)[1] - dim(ancestral_snvs)[1]

# For each cell, check for expression of ancestral mutations.
# Look up the chromosome and position for the single cells
ssm_1_all <- read.table(ssm_all_files[1], header=TRUE)
temp1 <- dplyr::left_join(sc1_raw, ssm_1_all, "ID")
length(unique(temp1$ID))
temp1.gr <- ConstructGranges(temp1$CHR, temp1$POS, width=0)
sc_1_ancestral_ret <- findOverlaps(temp1.gr, ancestral_snvs.gr)
sc_1_ancestral <- temp1[sc_1_ancestral_ret@from,]
head(sc_1_ancestral)
head(ancestral_snvs)
coverage_1_ancestral <- sc_1_ancestral %>% group_by(ID) %>% summarise(n_d = sum(d.x > 0), n_b = sum(d.x - a > 0))
# How many of the ancestral sites are covered by single cells?
sum(coverage_1_ancestral$n_b >= 2)
sum(coverage_1_ancestral$n_d >= 2)

# For each cell, check for expression of sample specific mutations.
# Start with sc1.
ssm1_specific_snvs <- ssm1[-ret1@from,]
ssm1_specific_snvs.gr <- ConstructGranges(ssm1_specific_snvs$CHROM, ssm1_specific_snvs$POS, width = 0)
sc_1_specific_ret <- findOverlaps(temp1.gr, ssm1_specific_snvs.gr)
sc_1_specific <- temp1[sc_1_specific_ret@from,]
coverage_1_specific <- sc_1_specific %>% group_by(ID) %>% summarise(n_d = sum(d.x > 0), n_b = sum(d.x - a > 0))
# How many of the ancestral sites are covered by single cells?
sum(coverage_1_specific$n_b >= 2)
sum(coverage_1_specific$n_d >= 2)

# sc2
ssm_2_all <- read.table(ssm_all_files[2], header=TRUE)
temp2 <- dplyr::left_join(sc2_raw, ssm_2_all, "ID")
length(unique(temp2$ID))
temp2.gr <- ConstructGranges(temp2$CHR, temp2$POS, width=0)
sc_2_ancestral_ret <- findOverlaps(temp2.gr, ancestral_snvs.gr)
sc_2_ancestral <- temp2[sc_2_ancestral_ret@from,]

ssm2_specific_snvs <- ssm2[-ret2@from,]
ssm2_specific_snvs.gr <- ConstructGranges(ssm2_specific_snvs$CHROM, ssm2_specific_snvs$POS, width = 0)
sc_2_specific_ret <- findOverlaps(temp2.gr, ssm2_specific_snvs.gr)
sc_2_specific <- temp2[sc_2_specific_ret@from,]
coverage_2_specific <- sc_2_specific %>% group_by(ID) %>% summarise(n_d = sum(d.x > 0), n_b = sum(d.x - a > 0))
# How many of the ancestral sites are covered by single cells?
sum(coverage_2_specific$n_b >= 2)
sum(coverage_2_specific$n_d >= 2)

# sc3
ssm_3_all <- read.table(ssm_all_files[3], header=TRUE)
temp3 <- dplyr::left_join(sc3_raw, ssm_3_all, "ID")
length(unique(temp3$ID))
length(ssm_3_all$ID)
temp3.gr <- ConstructGranges(temp3$CHR, temp3$POS, width=0)
sc_3_ancestral_ret <- findOverlaps(temp3.gr, ancestral_snvs.gr)
sc_3_ancestral <- temp3[sc_3_ancestral_ret@from,]

ssm3_specific_snvs <- ssm3[-ret3@from,]
ssm3_specific_snvs.gr <- ConstructGranges(ssm3_specific_snvs$CHROM, ssm3_specific_snvs$POS, width = 0)
sc_3_specific_ret <- findOverlaps(temp3.gr, ssm3_specific_snvs.gr)
sc_3_specific <- temp3[sc_3_specific_ret@from,]
coverage_3_specific <- sc_3_specific %>% group_by(ID) %>% summarise(n_d = sum(d.x > 0), n_b = sum(d.x - a > 0))
# How many of the ancestral sites are covered by single cells?
sum(is.na(coverage_3_specific))
sum(coverage_3_specific$n_b >= 2)
sum(coverage_3_specific$n_d >= 2)


ssm_1_all.gr <- ConstructGranges(ssm_1_all$CHR, ssm_1_all$POS, width = 0)
ssm_2_all.gr <- ConstructGranges(ssm_2_all$CHR, ssm_2_all$POS, width = 0)
ssm_3_all.gr <- ConstructGranges(ssm_3_all$CHR, ssm_3_all$POS, width = 0)
# Get overlap across all 3 samples.
ret12 <- findOverlaps(ssm_1_all.gr, ssm_2_all.gr)
ssm_12 <- ssm_1_all[ret12@from,]
ssm_12.gr <- ConstructGranges(ssm_12$CHR, ssm_12$POS, width = 0)
ret123 <- findOverlaps(ssm_12.gr, ssm_3_all.gr)
ssm_123 <- ssm_12[ret123@from,]
ssm_123.gr <- ConstructGranges(ssm_123$CHR, ssm_123$POS, width = 0)

# Sample specific mutations:
ret <- findOverlaps(ssm_1_all.gr, ssm_2_all.gr)
ssm_1_specific <- ssm_1_all[-ret@from,]
ssm_1_specific.gr <- ConstructGranges(ssm_1_specific$CHR, ssm_1_specific$POS, width = 0)
ret <- findOverlaps(ssm_1_specific.gr, ssm_3_all.gr)
ssm_1_specific <- ssm_1_specific[-ret@from,]
ssm_1_specific.gr <- ConstructGranges(ssm_1_specific$CHR, ssm_1_specific$POS, width = 0)
# Get single cell reads at ssm_1_specific
sc1_specific <- subset(sc1_raw, ID %in% ssm_1_specific$ID)
length(unique(sc1_specific$ID))
stats1_specific <- sc1_specific %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
sum(stats1_specific$n_b >= 2) # 209 sites have variant expression on sample 1 specific sites (>= 2 cells).
# For each cell, how many of the ssm_1_specific sites are expressed.
cell_1_specific <- sc1_specific %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
summary(cell_1_specific$n_b)

# Sample specific mutations:
ret <- findOverlaps(ssm_2_all.gr, ssm_1_all.gr)
ssm_2_specific <- ssm_2_all[-ret@from,]
ssm_2_specific.gr <- ConstructGranges(ssm_2_specific$CHR, ssm_2_specific$POS, width = 0)
ret <- findOverlaps(ssm_2_specific.gr, ssm_3_all.gr)
ssm_2_specific <- ssm_2_specific[-ret@from,]
ssm_2_specific.gr <- ConstructGranges(ssm_2_specific$CHR, ssm_2_specific$POS, width = 0)
length(ssm_2_specific.gr)
# Get single cell reads at ssm_1_specific
sc2_specific <- subset(sc2_raw, ID %in% ssm_2_specific$ID)
length(unique(sc2_specific$ID))
stats2_specific <- sc2_specific %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
sum(stats2_specific$n_b >= 2) # 359 sites have variant expression on sample 2 specific sites (>= 2 cells).
# For each cell, how many of the ssm_1_specific sites are expressed.
cell_2_specific <- sc2_specific %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
summary(cell_2_specific$n_b)

# Sample specific mutations:
ret <- findOverlaps(ssm_3_all.gr, ssm_1_all.gr)
ssm_3_specific <- ssm_3_all[-ret@from,]
ssm_3_specific.gr <- ConstructGranges(ssm_3_specific$CHR, ssm_3_specific$POS, width = 0)
ret <- findOverlaps(ssm_3_specific.gr, ssm_2_all.gr)
ssm_3_specific <- ssm_3_specific[-ret@from,]
ssm_3_specific.gr <- ConstructGranges(ssm_3_specific$CHR, ssm_3_specific$POS, width = 0)
length(ssm_3_specific.gr)
# Get single cell reads at ssm_1_specific
sc3_specific <- subset(sc3_raw, ID %in% ssm_3_specific$ID)
length(unique(sc3_specific$ID))
stats3_specific <- sc3_specific %>% group_by(ID) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
sum(stats3_specific$n_b >= 2) # 447 sites have variant expression on sample 3 specific sites (>= 2 cells).
# For each cell, how many of the ssm_1_specific sites are expressed.
cell_3_specific <- sc3_specific %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(d - a > 0))
summary(cell_3_specific$n_b)
