library(GenomicRanges)
library(ScRNAClone)
library(vcfR)

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")
exon_file <- system.file("extdata", "exons.bed", package = "ScRNAClone")
vcf_file <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/somatic.snvs.vcf.gz"
vcf_exon_out <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/somatic.snvs.exon.vcf"
laks_file <- "/Users/seonghwanjun/data/cell-line/phylo/ov2295_clone_snvs.csv"
titan_file <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/hmm/optimalClusterSolution/tumor_sample_1_cluster1.segs.txt"
titan_outfile <- "/Users/seonghwanjun/data/cell-line/bulk/A90554A/hmm/optimalClusterSolution/tumor_sample_1_cluster1.segs.edited.txt"

# Read Strelka vcf file for sample 1, filter by exon.
# Read exon bed file.
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
vcf.df <- read.table(vcf_file, header=F)
names(vcf.df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR")

# Load the SNVs from Laks.
Laks <- read.table(laks_file, header=T, sep=",")
Laks$loc <- paste(Laks$chrom, Laks$coord, sep=":")
Laks <- Laks[!duplicated(Laks$loc),]
sum(duplicated(Laks$loc))

# We first want to include the SNVs that fall on an exon.
# We want exons on Laks and use that to trim the vcf file and output it.

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
write.table(vcf.df.exon, vcf_exon_out, quote=F, row.names=F)

# We will clean TitanCNA output. Sometimes TitanCNA output contains rows like this:
# tumor_sample_1	20	902060	902060	...
# PhyloWGS's parser doesn't like it when the start and end position are the same -- probably because they use half intervals.
# We can effectively change the end position to be +1.
titan <- read.table(titan_file, header=T)
titan$End_Position.bp. <- titan$End_Position.bp. + 1
write.table(titan, titan_outfile, quote=F, row.names = F, sep="\t")

# Now, we will run the parser included in PhyloWGS to generate SSM and CNV data.

