library(GenomicRanges)
library(ScRNAClone)

# TODO: Location to GTF file to be provided as an argument.
gtf <- read.table("/Users/seonghwanjun/data/references/Homo_sapiens.GRCh37.75.gtf", header=F, sep="\t")
names(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
exon.gtf <- subset(gtf, feature == "exon")
dim(exon.gtf)

laks_snv <- read.table("data/Laks/ov2295_clone_snvs.csv", header=T, sep=",")
laks_snv$loc <- paste(laks_snv$chrom, laks_snv$coord, sep=":")
laks_snv <- laks_snv[!duplicated(laks_snv$loc),]
dim(laks_snv)

laks_snv.gr <- ConstructGranges(laks_snv$chrom, laks_snv$coord, width = 0)
exon.gtf.gr <- ConstructGranges(exon.gtf$seqname, exon.gtf$start, width = exon.gtf$end - exon.gtf$start)
overlaps <- findOverlaps(laks_snv.gr, exon.gtf.gr)
laks_snv_exon <- laks_snv[overlaps@from,]
laks_snv_exon <- laks_snv_exon[!duplicated(laks_snv_exon$loc),]
length(unique(laks_snv_exon$loc))

head(laks_snv_exon)

# Output laks_snv_exon.
# We will get read counts from the bulk and scRNA-seq at these locations.
exonic_snv_file <- "data/Laks/ov2295_clone_exonic_snvs.csv"
write.table(laks_snv_exon, exonic_snv_file, quote = F, row.names = F, col.names = T, sep=",")

# Generate BED file for Strelka SNV.
head(laks_snv_exon)
bed <- laks_snv_exon[,c("chrom", "coord")]
colnames(bed) <- c("chrom", "chromStart")
bed$chromEnd <- bed$chromStart
write.table(bed, "data/Laks/ov2295_clone_exonic_snvs.bed", quote=F, row.names = F)
