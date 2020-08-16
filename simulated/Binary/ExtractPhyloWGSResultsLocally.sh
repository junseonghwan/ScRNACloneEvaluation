#!/bin/bash

REP_PATH=$1

# Run below line from the terminal:
#conda activate py27


mkdir ${REP_PATH}/phylowgs/results/

python2 /Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/write_results.py results ${REP_PATH}/phylowgs/trees.zip ${REP_PATH}/phylowgs/results/results.summ.json.gz ${REP_PATH}/phylowgs/results/results.muts.json.gz ${REP_PATH}/phylowgs/results/results.mutass.zip
wait

mkdir ${REP_PATH}/phylowgs/results/trees
unzip -d ${REP_PATH}/phylowgs/results/trees ${REP_PATH}/phylowgs/results/results.mutass.zip
wait

gunzip ${REP_PATH}/phylowgs/results/results.summ.json.gz
wait

echo "End."
