#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 16
#SBATCH -t 100:00:00
#SBATCH -J Ichor
#SBATCH --mem 32G

IchorCNASnakefile=$1
ClusterSlurm=$2

echo "$(date) Running on: $(hostname)"

module load R/3.6.3-nsc1-gcc-7.3.0
wait

conda activate snakemake
wait

snakemake -s ${IchorCNASnakefile} --cluster-config ${ClusterSlurm} --cores 16 -j 16

echo "$(date)."
