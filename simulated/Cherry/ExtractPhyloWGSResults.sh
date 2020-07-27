#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem 4G
#SBATCH -J ProcessPhyloWGS

module load Python/2.7.14-nsc1-gcc-2018a-eb

REP_PATH=$1

cd ${REP_PATH}/phylowgs
mkdir results/
cd results
python2 /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/phylowgs/write_results.py results ../trees.zip results.summ.json.gz results.muts.json.gz results.mutass.zip
mkdir trees
unzip -d trees/ results.mutass.zip