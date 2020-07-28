#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem 4G
#SBATCH -J ProcessPhyloWGS

module load Python/2.7.14-nsc1-gcc-2018a-eb

REP_PATH=$1

mkdir ${REP_PATH}/phylowgs/results/

python2 /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/phylowgs/write_results.py results ${REP_PATH}/phylowgs/trees.zip ${REP_PATH}/phylowgs/results/results.summ.json.gz ${REP_PATH}/phylowgs/results/results.muts.json.gz ${REP_PATH}/phylowgs/results/results.mutass.zip
wait

mkdir ${REP_PATH}/phylowgs/results/trees
unzip -d ${REP_PATH}/phylowgs/results/trees ${REP_PATH}/phylowgs/results/results.mutass.zip
wait

gunzip ${REP_PATH}/phylowgs/results/results.summ.json.gz
wait

Rscript --vanilla Rscripts/ExtractPhyloWGSResults.R ${REP_PATH}/phylowgs

echo "End."
