#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J PrepareInputData
#SBATCH --mem 4G

echo "$(date) Running on: $(hostname)"

module load R/3.6.3-nsc1-gcc-7.3.0
wait

R_SCRIPT=$1
CONFIG_FILE=$2

Rscript --vanilla \
        ${R_SCRIPT} \
        ${CONFIG_FILE} 
wait

echo "$(date) Finished running script ${R_SCRIPT}."
