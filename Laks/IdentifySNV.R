library(ScRNAClone)
library(GenomicRanges)

CONFIG_FILE <- as.character(args[1])
CONFIG_FILE <- "/home/x_seoju/ScRNACloneEvaluation/"
config <- read.table(CONFIG_FILE, header=F, as.is = T)
names(config) <- c("Key", "Value")

snv_file <- as.character(config$Value[config$Key == "SNV"])
bulk_bam_files <- as.character(config$Value[config$Key == "BULK"])
bulk_bam_files <- strsplit(bulk_bam_files, ",")[[1]]
titan_cna_files <- as.character(config$Value[config$Key == "TITAN_CNA"])
titan_cna_files <- strsplit(titan_cna_files, ",")[[1]]
germline_vcf <- as.character(config$Value[config$Key == "GERMLINE_VCF"])
sc_bam_path <- as.character(config$Value[config$Key == "SC_BAM_PATH"])
MIN_DEPTH <- as.numeric(config$Value[config$Key == "MIN_DEPTH"])
MIN_CELLS <- as.numeric(config$Value[config$Key == "MIN_CELLS"])
output_path <- as.character(config$Value[config$Key == "OUTPUT_PATH"])

snvs <- read.table(snv_file, header=T, sep=",")
snvs$loc <- paste(snvs$chrom, snvs$coord, sep=":")
snvs <- snvs[!duplicated(snvs$loc),]
sum(duplicated(snvs$loc))

exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")
exons <- read.table(exon_file, header=F)
names(exons) <- c("CHR", "START", "END", "GENE")
exons$CHR <- as.character(exons$CHR)
chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")

snvs.gr <- ConstructGranges(snvs$chrom, snvs$coord, width = 0)
exons.gr <- ConstructGranges(exons$CHR, exons$START, width = exons$END - exons$START)
ret <- findOverlaps(snvs.gr, exons.gr)
snvs.exon <- snvs[ret@from,]
snvs.exon$gene <- exons[ret@to,"GENE"]
# TODO: What about the SNVs that fall on two genes?
# In theory, this doesn't seem possible but depending on the annotation, it could happen.
# For now, we will just ignore this.
snvs.exon <- snvs.exon[!duplicated(snvs.exon$loc),]
dim(snvs.exon)
head(snvs.exon)
snvs.exon$ID <- paste("s", 1:dim(snvs.exon)[1], sep="")
# This is not quite complete because we actually need to grab read counts from the BAMs.
snvs.exon.df <- data.frame(ID=snvs.exon$ID, CHR=snvs.exon$chrom, POS=snvs.exon$coord,
                           GENE=snvs.exon$gene, REF=snvs.exon$ref, ALT=snvs.exon$alt,
                           b = snvs.exon$alt_counts, d = snvs.exon$total_counts)

outfile <- "/Users/seonghwanjun/data/cell-line/bulk/phylo/ov2295_clone_snvs_exon.csv"
write.table(snvs.exon.df, outfile, row.names = F, col.names = T, quote = F, sep= "\t")
