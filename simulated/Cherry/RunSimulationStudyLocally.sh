#!/bin/bash

echo "$(date) Run begins."

SIMUL_PATH=$1
CONFIG_FILE=$2

Rscript --vanilla Rscripts/GenerateInputForSimulatedStudy.R \
	${SIMUL_PATH} 
wait

/Users/seonghwanjun/ScRNACloneEvaluation/BulkScRNAClone/run -c ${CONFIG_FILE}

echo "$(date) Run finished."
