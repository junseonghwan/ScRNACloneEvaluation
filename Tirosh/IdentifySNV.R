library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)
library(Rsamtools)
library(ScRNAClone)
library(vcfR)

CONFIG_FILE <- as.character(args[1])
CONFIG_FILE <- "/Users/seonghwanjun/ScRNACloneEvaluation/Tirosh/ConfigLocal.txt"
config <- read.table(CONFIG_FILE, header=F, as.is = T)
names(config) <- c("Key", "Value")

strelka_vcf_files <- as.character(config$Value[config$Key == "STRELKA_VCF"])
strelka_vcf_files <- strsplit(strelka_vcf_files, ",")[[1]]
bulk_bam_files <- as.character(config$Value[config$Key == "BULK"])
bulk_bam_files <- strsplit(bulk_bam_files, ",")[[1]]
titan_cna_files <- as.character(config$Value[config$Key == "TITAN_CNA"])
titan_cna_files <- strsplit(titan_cna_files, ",")[[1]]
germline_vcf <- as.character(config$Value[config$Key == "GERMLINE_VCF"])
sc_bam_path <- as.character(config$Value[config$Key == "SC_BAM_PATH"])
MIN_DEPTH <- as.numeric(config$Value[config$Key == "MIN_DEPTH"])
MIN_CELLS <- as.numeric(config$Value[config$Key == "MIN_CELLS"])
output_path <- as.character(config$Value[config$Key == "OUTPUT_PATH"])

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")
exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")
exons <- read.table(exon_file, header=F)
names(exons) <- c("CHR", "START", "END", "GENE")
exons$CHR <- as.character(exons$CHR)

n_samples <- length(bulk_bam_files)
if (n_samples != length(titan_cna_files)) {
    stop("Number of SNV files does not match the number of CNV files provided.")
}

ssm_outfile <- paste(output_path, "exon_ssm.txt", sep="/")

# First, we will retrieve all SNVs that fall on exons (no filter is used).
strelka_exon <- FilterSNVByExon(strelka_vcf_files[1], exon_file, chrs)
sum(strelka_exon$FILTER == "PASS") # SNVs on exon that pass the filter.
n_snvs <- dim(strelka_exon)[1]
print(paste("Total number of SNVS: ", n_snvs))
# Generate SSM data for input.
bulk <- GenerateSSMInput(strelka_exon, min_depth = MIN_DEPTH, phylo_wgs = FALSE)
# Generate CNV data.
cnv <- GenerateCNVInput(titan_cna_files[1], bulk)
# Combine cnv data into bulk
bulk$MajorCN <- cnv$MajorCN
bulk$MinorCN <- cnv$MinorCN
# Write to file.
bulk$ID <- paste("s", 1:n_snvs, sep="")
write.table(bulk, ssm_outfile, row.names = F, col.names = T, quote = F, sep= "\t")
