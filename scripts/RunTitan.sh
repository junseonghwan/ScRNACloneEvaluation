#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 16
#SBATCH -t 168:00:00
#SBATCH -J Titan
#SBATCH --mem 32G

TitanCNASnakefile=$1
ClusterSlurm=$2

echo "$(date) Running on: $(hostname)"

module load R/3.6.3-nsc1-gcc-7.3.0
wait

conda activate snakemake
wait

snakemake -s ${TitanCNASnakefile} --cluster-config ${ClusterSlurm} --cores 16 -j 16

echo "$(date)."
