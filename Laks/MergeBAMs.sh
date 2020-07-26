#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 16
#SBATCH -t 100:00:00
#SBATCH -J MergeBAM
#SBATCH --mem 32G

echo "$(date) Running on: $(hostname)."

BAM_LIST=$1
OUTPUT_BAM=$2

samtools merge -c -f -b ${BAM_LIST} ${OUTPUT_BAM} -@ 15

echo "$(date) Finished merging BAMs."
