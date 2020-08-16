#!/bin/bash -l

#SBATCH -A snic2019-3-292
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J Strelka
#SBATCH --mem 16G

echo "$(date) Running on: $(hostname)"

SAMPLE_DIR=$1
TUMOR_SAMPLE=$2
NORMAL_SAMPLE=$3
BED_FILE=$4

/home/x_seoju/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $SAMPLE_DIR/$NORMAL_SAMPLE \
    --tumorBam $SAMPLE_DIR/$TUMOR_SAMPLE \
    --referenceFasta /proj/sc_ml/regev/reference/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
    --runDir $SAMPLE_DIR \
    --calllRegions $BED_FILE \
	--exome
wait

$SAMPLE_DIR/runWorkflow.py -m local -j 8

echo "$(date) Finished running Strelka."
