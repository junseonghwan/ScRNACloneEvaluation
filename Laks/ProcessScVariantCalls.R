
library(GenomicRanges)
library(ScRNAClone)
library(vcfR)

MUT_MATRIX_OUTPUT_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/mut_matrix.txt"

# Get exonic SNVs.
EXON_SNV_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_wo_cn.txt"
snv.exon <- read.table(EXON_SNV_PATH, header=T)
snv.exon.gr <- ConstructGranges(snv.exon$CHR, snv.exon$POS, width = 0)
n_snvs <- dim(snv.exon)[1]

VCF_PATH <- "/Users/seonghwanjun/data/cell-line/smart-seq3/VariantCalls/"
VCF_FILES <- list.files(VCF_PATH, pattern=".vcf", full.names = FALSE)
n_cells <- length(VCF_FILES)
cell_mut_matrix <- matrix(0, n_cells, n_snvs)
cell_names <- rep("", n_cells)
for (i in 1:n_cells) {
    vcf_full_path <- paste(VCF_PATH, VCF_FILES[i], sep="/")
    vcf <- read.vcfR(vcf_full_path)
    fix.df <- as.data.frame(vcf@fix)
    if (dim(vcf@fix)[1] == 0) {
        next
    }
    fix.df$CHROM <- as.character(fix.df$CHROM)
    fix.df$POS <- as.numeric(as.character(fix.df$POS))
    vcf.gr <- ConstructGranges(fix.df$CHROM, fix.df$POS, width = 0)
    overlaps <- suppressWarnings(findOverlaps(snv.exon.gr, vcf.gr))
    if (length(overlaps) > 0) {
        cell_mut_matrix[i,overlaps@from] <- 1
    }
    cell_names[i] <- strsplit(VCF_FILES[i], split = ".vcf")[[1]]
}

rownames(cell_mut_matrix) <- cell_names
cell_mut_matrix <- as.data.frame(cell_mut_matrix)
colnames(cell_mut_matrix) <- snv.exon$ID
head(cell_mut_matrix)

remove_idx <- rowSums(cell_mut_matrix) > 0
cell_mut_matrix <- cell_mut_matrix[remove_idx,]
dim(cell_mut_matrix)
head(cell_mut_matrix)
write.table(cell_mut_matrix, MUT_MATRIX_OUTPUT_PATH, row.names = T, col.names = T, quote=F)

# We can prepare input for ddClone and B-SCITE.
# See RunDDClone.R for generating input and running ddClone.
# See GenerateInputForBScite.R for generating input for B-SCITE.

