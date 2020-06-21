args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args) != 3) {
  stop("Path to somatic vcf, sorted single cell bam file, and output path are necessary to run the program.")
}

print("Begin extracting reads in R.")

library(dplyr)
library(GenomicAlignments)
library(GenomicRanges)
library(reshape2)
library(Rsamtools)
library(vcfR)
library(ScRNAClone)

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")

#loci_file <- "/Users/seonghwanjun/data/cell-line/bulk/phylo/ov2295_clone_snvs_exon.csv"
#bam_file <- "/Users/seonghwanjun/data/cell-line/smartseq2/P15451_1094_S95_L001_R.Aligned.sortedByCoord.out.bam"
#sc_output_file <- "/Users/seonghwanjun/data/cell-line/smartseq2/Reads/P15451_1094_S95_L001_R.Aligned.sortedByCoord.out.reads.txt"

loci_file <- args[1]
bam_file <- args[2]
sc_output_file <- args[3]
print(bam_file)

snv_loci <- read.table(loci_file, header=T, sep="\t")
snv_loci$REF_COUNT <- 0
snv_loci$ALT_COUNT <- 0
print(dim(snv_loci))
print(names(snv_loci))
snv_loci$CHR <- as.character(snv_loci$CHR)

# if (startsWith(snv_loci$CHR[1], "chr")) {
#   snv_loci$CHR <- unlist(lapply(as.character(snv_loci$CHR), function(str) {
#     substr(str, 4, nchar(str))
#   }))
# }

somatic.gr <- ConstructGranges(snv_loci$CHR, snv_loci$POS, width = 0)
print(length(somatic.gr))

# get counts at somatic SNV loci
sbp <- ScanBamParam(which = somatic.gr)
p_param <- PileupParam(distinguish_nucleotides = TRUE,
                       distinguish_strands = FALSE,
                       include_insertions = TRUE)
res_somatic <- pileup(file=bam_file,
                      scanBamParam = sbp,
                      pileupParam = p_param)
if (dim(res_somatic)[1] > 0) {
  df <- data.frame(seqnames = as.character(res_somatic[,1]),
                  pos = as.numeric(as.character(res_somatic[,2])),
                  nucleotide = as.character(res_somatic[,3]),
                  count = as.numeric(as.character(res_somatic[,4])),
                  stringsAsFactors = F)
  df$nucleotide <- factor(df$nucleotide, nucleotides)
  df_somatic <- dcast(df, seqnames + pos ~ nucleotide, value.var = "count", fun.aggregate = sum, drop = FALSE)
  df_somatic$loc <- paste(df_somatic$seqnames, df_somatic$pos, sep=":")
  df$loc <- paste(df$seqnames, df$pos, sep=":")
  df_somatic <- semi_join(df_somatic, df, by = "loc")
  df_somatic <- df_somatic[,-which(names(df_somatic) == "loc")]

  if (sum(colnames(df_somatic) == "-") > 0) {
    col_to_remove<-which(colnames(df_somatic) == "-")
    df_somatic<-df_somatic[,-col_to_remove]
  }
  n <- dim(df_somatic)[2]
  df_somatic$Total <- rowSums(df_somatic[,3:n])

  nucleotide_counts <- rowSums(df_somatic[,3:n] > 0)
  print("How many biallelic sites?")
  print(paste(sum(nucleotide_counts > 1), sum(nucleotide_counts > 0), sep="/"))

  # get REF and VAR counts
  nn <- dim(df_somatic)[1]
  idx <- paste(snv_loci$CHR, snv_loci$POS, sep=":") %in% paste(df_somatic$seqnames, df_somatic$pos, sep=":")
  for (i in 1:nn) {
    idx <- which(snv_loci$CHR == df_somatic[i,"seqnames"] & snv_loci$POS == df_somatic[i,"pos"])
    if (length(idx) > 0) {
      alt <- as.character(snv_loci[idx,"ALT"])
      ref <- as.character(snv_loci[idx,"REF"])
      if ((ref %in% nucleotides) & (alt %in% nucleotides)) {
        alt_idx <- which(alt == nucleotides)
        ref_idx <- which(ref == nucleotides)
        snv_loci[idx,"REF_COUNT"] <- df_somatic[i,2+ref_idx]
        snv_loci[idx,"ALT_COUNT"] <- df_somatic[i,2+alt_idx]
      }
    }
  }
}
write.table(snv_loci[,c("ID", "REF_COUNT", "ALT_COUNT")], file=sc_output_file, sep=",", row.names=F, quote=F)
