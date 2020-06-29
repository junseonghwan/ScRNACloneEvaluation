#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J AlleleCounts
#SBATCH --mem 8G

echo "$(date) Running on: $(hostname)"

BAM_FILE=$1
OUTPUT_PATH=$2

/home/x_seoju/bcftools-1.10.2/bcftools mpileup \
	-Ou \
	-f /proj/sc_ml/Tirosh/reference/GRCh37/Ensembl/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	${BAM_FILE_NAME} | bcftools call -P 0.1 -mv -Oz -o ${OUTPUT_PATH}.vcf.gz

echo "$(date) Finished variants."
