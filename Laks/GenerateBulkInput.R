library(dplyr)
library(GenomicRanges)
library(reshape2)
library(Rsamtools)
library(ScRNAClone)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)
print(args)
CONFIG_FILE <- as.character(args[1])
CONFIG_FILE <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"


# Scripts to prepare bulk input for scRNA + bulk inference.
###################
# 1. Identify SNVs on exons.
# 2. Extract genotype information at the identified SNVs from TitanCNA.
###################
# 3. Extract reads from single cell BAM files.
###################
# 4. Combine single cell reads extracted.
# 5. Generate hyperparameters.
###################

exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")
chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")

config <- read.table(CONFIG_FILE, header=F, as.is = T)
print(config)
names(config) <- c("Key", "Value")
strelka_vcf_files <- as.character(config$Value[config$Key == "STRELKA_VCF"])
strelka_vcf_files <- strsplit(strelka_vcf_files, ",")[[1]]
titan_cna_files <- as.character(config$Value[config$Key == "TITAN_CNA"])
titan_cna_files <- strsplit(titan_cna_files, ",")[[1]]
germline_vcf <- as.character(config$Value[config$Key == "GERMLINE_VCF"])
sc_bam_path <- as.character(config$Value[config$Key == "SC_BAM_PATH"])
MIN_DEPTH <- as.numeric(config$Value[config$Key == "MIN_DEPTH"])
MIN_CELLS <- as.numeric(config$Value[config$Key == "MIN_CELLS"])
output_path <- as.character(config$Value[config$Key == "OUTPUT_PATH"])

n_samples <- length(strelka_vcf_files)
if (n_samples != length(titan_cna_files)) {
    stop("Number of SNV files does not match the number of CNV files provided.")
}

ssm_outfile <- paste(output_path, "combined_ssm.txt", sep="/")

bulk_samples <- list()
bulk.grs <- list()
for (n in 1:n_samples) {
    # First, we will retrieve all SNVs that fall on exons (no filter is used).
    strelka_exon <- FilterSNVByExon(strelka_vcf_files[n], exon_file, chrs)
    n_snvs <- dim(strelka_exon)[1]
    print(paste("Total number of SNVS: ", n_snvs))
    # Generate SSM data for input.
    bulk <- GenerateSSMInput(strelka_exon, min_depth = MIN_DEPTH, phylo_wgs = FALSE)
    # Generate CNV data.
    cnv <- GenerateCNVInput(titan_cna_files[n], bulk)
    # Combine cnv data into bulk.
    bulk$MajorCN <- cnv$MajorCN
    bulk$MinorCN <- cnv$MinorCN
    # Write to file.
    out_file <- paste(output_path, "/ssm", n, ".txt", sep="")
    write.table(bulk, out_file, row.names = F, col.names = T, quote = F, sep= "\t")
    # Write the raw VCF as a text file.
    out_file <- paste(output_path, "/strelka", n, ".txt", sep="")
    write.table(strelka_exon, out_file, row.names = F, col.names = T, quote = F, sep= "\t")
    out_file <- paste(output_path, "/strelka_pass", n, ".txt", sep="")
    write.table(subset(strelka_exon, FILTER == "PASS"), out_file, row.names = F, col.names = T, quote = F, sep= "\t")

    bulk_samples[[n]] <- bulk
    bulk.grs[[n]] <- ConstructGranges(bulk_samples[[n]]$CHR, start = bulk_samples[[n]]$POS, width = 0)
}

# Find intersection across the samples -- do a table join.
bulk.gr <- bulk.grs[[1]]
bulk <- bulk_samples[[1]]
for (n in 2:n_samples) {
    ret <- findOverlaps(bulk.gr, bulk.grs[[n]])
    temp1 <- bulk[ret@from,]
    temp2 <- bulk_samples[[n]][ret@to,]
    temp1$b <- paste(temp1$b, temp2$b, sep=",")
    temp1$d <- paste(temp1$d, temp2$d, sep=",")
    temp1$MajorCN <- paste(temp1$MajorCN, temp2$MajorCN, sep=",")
    temp1$MinorCN <- paste(temp1$MinorCN, temp2$MinorCN, sep=",")
    bulk <- temp1
    bulk.gr <- ConstructGranges(bulk$CHR, start = bulk$POS, width = 0)
}
dim(bulk)

# Load VCF from HaplotypeCaller.
# We will eliminate SNVs that appear as germline.
# ret <- read.vcfR(germline_vcf)
# fix <- getFIX(ret)
# fix.df <- data.frame(fix)
# fix.df$POS <- as.numeric(as.character(fix.df$POS))
# fix.df$CHROM <- as.character(fix.df$CHROM)
# germline.gr <- ConstructGranges(fix.df$CHROM, start = fix.df$POS, width = 0)
# somatic.gr <- ConstructGranges(bulk$CHR, start = bulk$POS, width = 0)
# ret2 <- findOverlaps(somatic.gr, germline.gr)
# bulk2 <- bulk[-ret2@from,]

#n_snvs <- dim(bulk2)[1]
#bulk2$ID <- paste("s", 1:n_snvs, sep="")
#write.table(bulk2, ssm_outfile, row.names = F, col.names = T, quote = F, sep= "\t")

n_snvs <- dim(bulk)[1]
bulk$ID <- paste("s", 1:n_snvs, sep="")
write.table(bulk, ssm_outfile, row.names = F, col.names = T, quote = F, sep= "\t")

# At this point, we need to process single cells at the loci identified in bulk.
# The script to use is in scripts/BatchExtractReadCounts.py.
# It is a python script that will launch SLURM jobs, one for each cell.
# Each job will execute ProcessSingleCellBAMScript.R.
# Once completed, use GenerateScInput.R to combine single cell reads and estimates hyperparams.
