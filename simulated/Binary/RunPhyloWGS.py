# script to launch PhyloWGS experiments
import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

parser = argparse.ArgumentParser(description='Script to batch submit simulation studies.')
parser.add_argument('seed', help="Seed for random.", type=int)
parser.add_argument('case_path', help="Path to simulation path.", type=str)
parser.add_argument('--sbegin', help="Sim begin.", type=int, default = 0)
parser.add_argument('--send', help="Sim end (not inclusive).", type=int, default = 10)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 10)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 2500)
parser.add_argument('-p', '--mh', help="Number of MH iterations.", type=int, default = 2000)
parser.add_argument('-u', '--burn_in', help="Burn_in.", type=int, default = 1000)
parser.add_argument('--local', help="Run it locally.", action="store_true")
args = parser.parse_args()

if args.local:
    PATH_TO_EXECUTABLE  =  "RunPhyloWGSLocally.sh" # path to executable
else:
    PATH_TO_EXECUTABLE  =  "RunPhyloWGS.sh" # path to executable
SEED = args.seed
CASE_PATH = args.case_path
SIM_BEGIN = args.sbegin
SIM_END = args.send
REP_BEGIN = args.rbegin
REP_END = args.rend
N_MCMC_ITER = args.mcmc
N_MH_ITER = args.mh
BURN_IN = args.burn_in

np.random.seed(SEED)

for sim_no in range(SIM_BEGIN, SIM_END):
    SIMUL_PATH = CASE_PATH + "/sim" + str(sim_no)
    for rep_no in range(REP_BEGIN, REP_END):
        REP_PATH = SIMUL_PATH + "/rep" + str(rep_no) + "/"
        rep_seed = np.random.randint(100000000)
        rep_pwgs_bulk_path = REP_PATH + "pwgs_snv.txt"
        rep_pwgs_cnv_path = REP_PATH + "pwgs_cnv.txt"
        rep_output_path = REP_PATH + "phylowgs/"

        if not os.path.exists(rep_output_path):
            os.makedirs(rep_output_path)

        if args.local:
            run_command = "./" + PATH_TO_EXECUTABLE + " "
        else:
            run_command = "sbatch "
            run_command += PATH_TO_EXECUTABLE + " "

        run_command += REP_PATH + " "
        run_command += rep_pwgs_bulk_path + " "
        run_command += rep_pwgs_cnv_path + " "
        run_command += rep_output_path + " "
        run_command += str(BURN_IN) + " "
        run_command += str(N_MCMC_ITER) + " "
        run_command += str(N_MH_ITER) + " "
        run_command += str(rep_seed) + " "

        #print(run_command)
        os.system(run_command)
