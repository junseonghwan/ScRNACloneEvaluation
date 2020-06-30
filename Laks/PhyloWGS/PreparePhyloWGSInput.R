args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

print(args)

library(GenomicRanges)
library(ScRNAClone)
library(vcfR)

library(dplyr)
library(reshape2)
library(Rsamtools)

CONFIG_FILE <- as.character(args[1])
#CONFIG_FILE <- "/home/x_seoju/ScRNACloneEvaluation/Laks/PhyloWGS/ConfigPWGS.txt"
#CONFIG_FILE <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/PhyloWGS/ConfigPWGSLocal.txt"
config <- read.table(CONFIG_FILE, header=F, as.is = T)
names(config) <- c("Key", "Value")
VCF_PATH <- as.character(config$Value[config$Key == "VCF_PATH"])
LAKS_SNV_FILE <- as.character(config$Value[config$Key == "LAKS_SNV_PATH"])
TITAN_CNA_FILE <- as.character(config$Value[config$Key == "TITAN_CNA"])
OUTPUT_PATH <- as.character(config$Value[config$Key == "OUTPUT_PATH"])
vcf_exon_outfile <- paste(OUTPUT_PATH, "PWGS", "somatic.snvs.exon.vcf", sep="/")
trimmed_vcf_exon_outfile <- paste(OUTPUT_PATH, "PWGS", "somatic.snvs.exon.trimmed.vcf", sep="/")
titan_outfile <- paste(OUTPUT_PATH, "PWGS", "titan.modified.segs.txt", sep="/")

if (!dir.exists(paste(OUTPUT_PATH, "PWGS", sep="/"))) {
    dir.create(paste(OUTPUT_PATH, "PWGS", sep="/"))
}

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")
exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")

# Read exon bed file.
# Read Strelka vcf file and filter it by exon.s
# Read SNV list file from Laks et al.
# Take intersection to get SNVs on exons from Laks et al.
# Subset VCF file and output it.
# Outside of this script, run PhyloWGS's parser to generate input files for inference.

# Load the exons.
exons <- read.table(exon_file, header=F)
names(exons) <- c("CHR", "START", "END", "GENE")
exons$CHR <- as.character(exons$CHR)

# Load the VCF file.
# Reading will take few minutes.
vcf.df <- read.table(VCF_PATH, header=F)
names(vcf.df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR")

# Load the SNVs from Laks.
Laks <- read.table(LAKS_SNV_FILE, header=T, sep=",")
Laks$loc <- paste(Laks$chrom, Laks$coord, sep=":")
Laks <- Laks[!duplicated(Laks$loc),]
sum(duplicated(Laks$loc))
dim(Laks)

# We want to include the SNVs that fall on an exon.
# We want exons on Laks and use that to trim the vcf file and then output the resulting data.

# Get exon from Laks.
Laks.gr <- ConstructGranges(Laks$chrom, Laks$coord, width = 0)
exons.gr <- ConstructGranges(exons$CHR, exons$START, width = exons$END - exons$START)
ret <- findOverlaps(Laks.gr, exons.gr)
Laks.exon.snv <- Laks[ret@from,]
# Note: duplicates can exist since an exon may fall on two genes with overlapping boundaries.
Laks.exon.snv <- Laks.exon.snv[!duplicated(Laks.exon.snv$loc),]
dim(Laks.exon.snv)

# Now, let's trim the vcf file.
vcf.df.gr <- ConstructGranges(vcf.df$CHROM, vcf.df$POS, width = 0)
Laks.exon.snv.gr <- ConstructGranges(Laks.exon.snv$chrom, Laks.exon.snv$coord, width = 0)
ret <- findOverlaps(Laks.exon.snv.gr, vcf.df.gr)
vcf.df.exon <- vcf.df[ret@to,]
dim(vcf.df.exon)

# Write to file -- this file can be process by PhyloWGS's parser.
write.table(vcf.df.exon, vcf_exon_outfile, quote=F, row.names=F)

# We will clean TitanCNA output. Sometimes TitanCNA output contains rows where start pos = end pos.
# tumor_sample_1	20	902060	902060	...
# PhyloWGS's parser doesn't like it when the start and end position are the same -- probably because they use half intervals.
# We will change the end position to be +1.
titan <- read.table(TITAN_CNA_FILE, header=T)
titan$End_Position.bp. <- titan$End_Position.bp. + 1
write.table(titan, titan_outfile, quote=F, row.names = F, sep="\t")

# Now, we should run the parser included in PhyloWGS to generate SSM and CNV data.

# We will also generate a trimmed VCF file based on SNVs that have single cell coverage.
trimmed_snvs <- read.table("/Users/seonghwanjun/data/cell-line/bulk/Combined/trimmed_ssm.txt", header=T)
trimmed_snvs.gr <- ConstructGranges(trimmed_snvs$CHR, trimmed_snvs$POS, width = 0)
overlaps <- findOverlaps(vcf.df.gr, trimmed_snvs.gr)
vcf.df.exon.trimmed <- vcf.df[overlaps@from,]
write.table(vcf.df.exon.trimmed, trimmed_vcf_exon_outfile, quote=F, row.names=F)
