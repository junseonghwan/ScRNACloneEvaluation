library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)
library(Rsamtools)
library(ScRNAClone)
library(vcfR)

CONFIG_FILE <- as.character(args[1])
CONFIG_FILE <- "/home/x_seoju/ScRNACloneEvaluation/Laks/Config.txt"
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
OUTPUT_PATH <- as.character(config$Value[config$Key == "OUTPUT_PATH"])

n_samples <- length(bulk_bam_files)

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

# Get read counts from the BAM files.
somatic.gr <- ConstructGranges(snvs.exon$chrom, snvs.exon$coord, width = 0)
sbp <- ScanBamParam(which = somatic.gr)
p_param <- PileupParam(distinguish_nucleotides = TRUE,
                       distinguish_strands = FALSE,
                       include_insertions = FALSE,
                       max_depth = 100000)

snvs.exon.dfs <- list()
for(n in 1:n_samples) {
    snvs.exon.df <- data.frame(ID=snvs.exon$ID, CHR=snvs.exon$chrom, POS=snvs.exon$coord,
                               GENE=snvs.exon$gene, REF=snvs.exon$ref, ALT=snvs.exon$alt,
                               a = 0, b = 0, d = 0, MajorCN = 0, MinorCN = 0)
    res_somatic <- pileup(file=bulk_bam_files[n],
                      scanBamParam = sbp,
                      pileupParam = p_param)
    res_somatic$loc <- paste(res_somatic$seqnames, res_somatic$pos, sep=":")
    # TODO: pileup returns the results in ordered fashion but in case it doesn't the below code will
    # be errerneous.
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
    cnv <- GenerateCNVInput(titan_cna_files[n], snvs.exon.df)
    snvs.exon.df$MajorCN <- cnv$MajorCN
    snvs.exon.df$MinorCN <- cnv$MinorCN
    snvs.exon.dfs[[n]] <- snvs.exon.df
}

# Now, merge snvs.exons.dfs into one data frame.
snvs.exon.df <- snvs.exon.dfs[[1]]
for (n in 2:n_samples) {
    snvs.exon.df$a <- paste(snvs.exon.df$a, snvs.exon.dfs[[n]]$a, sep=",")
    snvs.exon.df$b <- paste(snvs.exon.df$b, snvs.exon.dfs[[n]]$b, sep=",")
    snvs.exon.df$d <- paste(snvs.exon.df$d, snvs.exon.dfs[[n]]$d, sep=",")
    snvs.exon.df$MajorCN <- paste(snvs.exon.df$MajorCN, snvs.exon.dfs[[n]]$MajorCN, sep=",")
    snvs.exon.df$MinorCN <- paste(snvs.exon.df$MinorCN, snvs.exon.dfs[[n]]$MinorCN, sep=",")
}

outfile <- paste(OUTPUT_PATH, "ov2295_clone_snvs_exon.csv", sep="/")
# Remove column "a".
snvs.exon.df <- snvs.exon.df[,-which(names(snvs.exon.df) == "a")]
write.table(snvs.exon.df, outfile, row.names = F, col.names = T, quote = F, sep= "\t")
