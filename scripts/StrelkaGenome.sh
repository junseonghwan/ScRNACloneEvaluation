#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 8
#SBATCH -t 100:00:00
#SBATCH -J StrelkaGenome
#SBATCH --mem 32G

echo "$(date) Running on: $(hostname)"

SAMPLE_DIR=$1
TUMOR_BAM=$2
NORMAL_BAM=$3

/home/x_seoju/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $NORMAL_BAM \
    --tumorBam $TUMOR_BAM \
    --referenceFasta /proj/sc_ml/Tirosh/reference/GRCh37/Ensembl/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --runDir $SAMPLE_DIR
wait

$SAMPLE_DIR/runWorkflow.py -m local -j 8

echo "$(date) Finished running Strelka."
