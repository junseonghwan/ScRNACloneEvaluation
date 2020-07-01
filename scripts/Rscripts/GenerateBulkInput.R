args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

print(args)

library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)
library(Rsamtools)
library(ScRNAClone)
library(vcfR)

CONFIG_FILE <- as.character(args[1])
#CONFIG_FILE <- "/home/x_seoju/ScRNACloneEvaluation/Laks/Config.txt"
#CONFIG_FILE <- "/Users/seonghwanjun/ScRNACloneEvaluation/Laks/ConfigLocal.txt"
config <- read.table(CONFIG_FILE, header=F, as.is = T)
names(config) <- c("Key", "Value")

#SNV_FILE <- as.character(config$Value[config$Key == "SNV_PATH"])
EXON_SNV_FILE <- as.character(config$Value[config$Key == "EXON_SNV_PATH"])
BULK_BAM_FILE <- as.character(config$Value[config$Key == "BULK"])
TITAN_CNA_FILE <- as.character(config$Value[config$Key == "TITAN_CNA"])
#MIN_DEPTH <- as.numeric(config$Value[config$Key == "MIN_DEPTH"])
BULK_OUTPUT_FILE <- as.character(config$Value[config$Key == "SNV_PATH"])

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")

snvs.exon <- read.table(EXON_SNV_FILE, header=T, sep=",")
snvs.exon$loc <- paste(snvs.exon$chrom, snvs.exon$coord, sep=":")
snvs.exon$ID <- paste("s", 1:dim(snvs.exon)[1], sep="")
sum(duplicated(snvs.exon$loc)) # Check that there aren't any duplicates.

# exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")
# exons <- read.table(exon_file, header=F)
# names(exons) <- c("CHR", "START", "END", "GENE")
# exons$CHR <- as.character(exons$CHR)

# TODO: Provide feature for specifying min_base_quality, min_mapq, and min_depth in Config.txt.
# Get read counts from the BAM files.
somatic.gr <- ConstructGranges(snvs.exon$chrom, snvs.exon$coord, width = 0)
sbp <- ScanBamParam(which = somatic.gr)
p_param <- PileupParam(distinguish_nucleotides = TRUE,
                       distinguish_strands = FALSE,
                       include_insertions = FALSE,
                       min_base_quality = 20, min_mapq = 20,
                       max_depth = 100000)

snvs.exon.df <- data.frame(ID=snvs.exon$ID, CHR=snvs.exon$chrom, POS=snvs.exon$coord,
                           REF=snvs.exon$ref, ALT=snvs.exon$alt,
                           a = 0, b = 0, d = 0, MajorCN = 0, MinorCN = 0)
res_somatic <- pileup(file=BULK_BAM_FILE,
                      scanBamParam = sbp,
                      pileupParam = p_param)
res_somatic$loc <- paste(res_somatic$seqnames, res_somatic$pos, sep=":")

# TODO: pileup returns the results ordered by chr then position but in case it doesn't, the below code may not work.
res_somatic$loc <- factor(res_somatic$loc, levels = unique(res_somatic$loc))
temp <- dcast(res_somatic, loc ~ nucleotide, value.var = "count")
temp2 <- temp[,1:5]
temp2[is.na(temp2)] <- 0
snvs.exon.df$loc <- paste(snvs.exon.df$CHR, snvs.exon.df$POS, sep=":")
# Check that the locations match -- if this value is not 1, there is a problem.
mean(snvs.exon.df$loc == temp2$loc)
for (base in nucleotides) {
    idx <- which(snvs.exon.df$REF == base)
    snvs.exon.df[idx,"a"] <- temp2[idx,base]
    idx <- which(snvs.exon.df$ALT == base)
    snvs.exon.df[idx,"b"] <- temp2[idx,base]
}
snvs.exon.df$d <- snvs.exon.df$a + snvs.exon.df$b

# Get copy number profile.
# Generate CNV data.
cnv <- GenerateCNVInput(TITAN_CNA_FILE, snvs.exon.df)
snvs.exon.df$MajorCN <- cnv$MajorCN
snvs.exon.df$MinorCN <- cnv$MinorCN

# Remove column "a".
snvs.exon.df <- snvs.exon.df[,-which(names(snvs.exon.df) == "a")]
# This will output bulk file for all exonic snvs.
# We will retrieve reads from single cells and remove any SNVs without any cell coverage.
write.table(snvs.exon.df, BULK_OUTPUT_FILE, row.names = F, col.names = T, quote = F, sep= "\t")
