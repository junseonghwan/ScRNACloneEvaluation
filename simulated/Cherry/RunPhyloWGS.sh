#!/bin/bash

#module load Python/3.7.0-anaconda-5.3.0-extras-nsc1
module load Python/2.7.14-nsc1-gcc-2018a-eb

conda activate py27
wait

python /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/phylowgs/evolve.py $1 $2 -O $3 -B $4 -s $5 -i $6 -r $7
