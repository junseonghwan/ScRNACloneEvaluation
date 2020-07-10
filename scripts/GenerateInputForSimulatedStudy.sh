#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem 8G
#SBATCH -J GenerateInputForSimulation

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SIMUL_PATH=$1

Rscript --vanilla Rscripts/GenerateInputForSimulatedStudy.R \
	${SIMUL_PATH}

echo "$(date) Run finished."
