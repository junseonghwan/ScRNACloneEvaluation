#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 3
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J scRNA

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

CONFIG_FILE=$1

./../BulkScRNAClone/run -c ${CONFIG_FILE}

echo "$(date) Run finished."
