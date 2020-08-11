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
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 1000)
parser.add_argument('-p', '--mh', help="Number of MH iterations.", type=int, default = 2000)
parser.add_argument('-t', '--thinning', help="Thinning interval.", type=int, default = 10)
parser.add_argument('-u', '--burn_in', help="Burn_in.", type=int, default = 100)
parser.add_argument('-o', '--overwrite', help="Overwrite main.config if it exists.", action="store_true")
parser.add_argument('--geometric', help="Use geometric mean for single cell likelihood.", action="store_true")
parser.add_argument('--genotype', help="Use genotype information.", action="store_true")
parser.add_argument('--local', help="Run it locally.", action="store_true")
args = parser.parse_args()

if args.local:
    PATH_TO_EXECUTABLE = "RunSimulationStudyLocally.sh"
else:
    PATH_TO_EXECUTABLE  =  "RunSimulationStudy.sh" # path to executable
SEED = args.seed
CASE_PATH = args.case_path
SIM_BEGIN = args.sbegin
SIM_END = args.send
REP_BEGIN = args.rbegin
REP_END = args.rend
N_MCMC_ITER = args.mcmc
N_MH_ITER = args.mh
THINNING = args.thinning
BURN_IN = args.burn_in

# Default values.
ALPHA0_MAX = 10 
LAMBDA_MAX = 0.8
GAMMA_MAX = 0.3
ALPHA0_MIN = 0.1
LAMBDA_MIN = 0.1
GAMMA_MIN = 0.1
SEQ_ERR = 0.001
VAR_CP_PROB = 0.5
SC_DROPOUT_ALPHA0 = 0.5
SC_DROPOUT_BETA0 = 1
SC_BURSTY_ALPHA0 = 1
SC_BURSTY_BETA0 = 0.01
GEOMETRIC_MEAN = (1 if args.geometric else 0)

np.random.seed(SEED)

for sim_no in range(SIM_BEGIN, SIM_END):
    SIMUL_PATH = CASE_PATH + "/sim" + str(sim_no)
    for rep_no in range(REP_BEGIN, REP_END):
        REP_PATH = SIMUL_PATH + "/rep" + str(rep_no) + "/"
        # write main.config file to REP_PATH
        config_file_path = REP_PATH + "/main.config"
        if not os.path.exists(config_file_path) or args.overwrite:
            print("Writing a new main.config file for " + REP_PATH)
            rep_seed = np.random.randint(100000000)
            if args.genotype:
                rep_bulk_path = REP_PATH + "genotype_ssm.txt"
                output_path = REP_PATH + "/genotype/"
            else:
                rep_bulk_path = REP_PATH + "simul_ssm.txt"
                output_path = REP_PATH + "/total_cn/"
            rep_scRNA_path = REP_PATH + "simul_sc.txt"
            rep_hyper_params_path = REP_PATH + "simul_sc_hp.txt"

            f = open(config_file_path, "w")
            f.write("seed: " + str(rep_seed) + "\n")
            f.write("bulk_data_path: " + rep_bulk_path + "\n")
            f.write("sc_rna_data_path: " + rep_scRNA_path + "\n")
            f.write("hyperparams_path: " + rep_hyper_params_path + "\n")
            f.write("output_path: " + output_path + "\n")
            f.write("n_mcmc_iter: " + str(N_MCMC_ITER) + "\n")
            f.write("n_mh_iter: " + str(N_MH_ITER) + "\n")
            f.write("thinning: " + str(THINNING) + "\n")
            f.write("burn_in: " + str(BURN_IN) + "\n")
            f.write("alpha0_max: " + str(ALPHA0_MAX) + "\n")
            f.write("lambda_max: " + str(LAMBDA_MAX) + "\n")
            f.write("gamma_max: " + str(GAMMA_MAX) + "\n")
            f.write("alpha0_min: " + str(ALPHA0_MIN) + "\n")
            f.write("lambda_min: " + str(LAMBDA_MIN) + "\n")
            f.write("gamma_min: " + str(GAMMA_MIN) + "\n")
            f.write("seq_err: " + str(SEQ_ERR) + "\n")
            f.write("var_cp_prob: " + str(VAR_CP_PROB) + "\n")
            f.write("sc_dropout_alpha0: " + str(SC_DROPOUT_ALPHA0) + "\n")
            f.write("sc_dropout_beta0: " + str(SC_DROPOUT_BETA0) + "\n")
            f.write("sc_bursty_alpha0: " + str(SC_BURSTY_ALPHA0) + "\n")
            f.write("sc_bursty_beta0: " + str(SC_BURSTY_BETA0) + "\n")
            f.write("geometric_mean: " + str(GEOMETRIC_MEAN) + "\n")
            f.close()

        if args.local:
            run_command = "./" + PATH_TO_EXECUTABLE + " "
        else:
            run_command = "sbatch "
            run_command += PATH_TO_EXECUTABLE + " "

        run_command += REP_PATH + " "
        run_command += config_file_path
        #print(run_command)
        os.system(run_command)
