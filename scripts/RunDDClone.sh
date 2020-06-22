#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J ddClone

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SIMUL_PATH=$1
MCMC_ITER=$2
SEED=$3

Rscript --vanilla Rscripts/RunDDClone.R \
	${SIMUL_PATH} ${MCMC_ITER} ${SEED}

echo "$(date) Run finished."
