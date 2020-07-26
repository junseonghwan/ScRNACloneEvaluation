import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "RunLaks.sh" # path to executable
parser = argparse.ArgumentParser(description='Script to batch submit replicates generated from Laks data.')
parser.add_argument('seed', help="Seed for random.", type=int)
parser.add_argument('data_path', help="Path to replicates.", type=str)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 1)
parser.add_argument('--rend', help="Rep end (inclusive).", type=int, default = 20)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 2000)
parser.add_argument('-p', '--mh', help="Number of MH iterations.", type=int, default = 2000)
parser.add_argument('-t', '--thinning', help="Thinning interval.", type=int, default = 20)
parser.add_argument('-u', '--burn_in', help="Burn_in.", type=int, default = 100)
parser.add_argument('-o', '--overwrite', help="Overwrite main.config if it exists.", action="store_true")
parser.add_argument('--w_cells', help="Use cells.", action="store_true")
args = parser.parse_args()

SEED = args.seed
DATA_PATH = args.data_path
REP_BEGIN = args.rbegin
REP_END = args.rend
N_MCMC_ITER = args.mcmc
N_MH_ITER = args.mh
THINNING = args.thinning
BURN_IN = args.burn_in

# Default values.
ALPHA0_MAX = 50 
LAMBDA_MAX = 1
GAMMA_MAX = 10
SEQ_ERR = 0.001
VAR_CP_PROB = 0.25
SC_DROPOUT_ALPHA0 = 0.01
SC_DROPOUT_BETA0 = 0.01

np.random.seed(SEED)

for rep_no in range(REP_BEGIN, REP_END+1):
    REP_PATH = DATA_PATH + "/rep" + str(rep_no) + "/"
    # write main.config file to REP_PATH
    if args.w_cells:
        config_file_path = REP_PATH + "/main.config"
    else:
        config_file_path = REP_PATH + "/main_no_sc.config"
    if not os.path.exists(config_file_path) or args.overwrite:
        print("Writing a new " + config_file_path + " file at " + REP_PATH)
        rep_seed = np.random.randint(100000000)
        rep_bulk_path = REP_PATH + "snv.txt"
        
        f = open(config_file_path, "w")
        f.write("seed: " + str(rep_seed) + "\n")
        f.write("bulk_data_path: " + rep_bulk_path + "\n")
        if args.w_cells:
            output_path = REP_PATH + "/w_cells/"
            rep_scRNA_path = REP_PATH + "sc.txt"
            rep_hyper_params_path = REP_PATH + "sc_hp.txt"
            f.write("sc_rna_data_path: " + rep_scRNA_path + "\n")
            f.write("hyperparams_path: " + rep_hyper_params_path + "\n")
        else:
            output_path = REP_PATH + "/wo_cells/"

        f.write("output_path: " + output_path + "\n")
        f.write("n_mcmc_iter: " + str(N_MCMC_ITER) + "\n")
        f.write("n_mh_iter: " + str(N_MH_ITER) + "\n")
        f.write("thinning: " + str(THINNING) + "\n")
        f.write("burn_in: " + str(BURN_IN) + "\n")
        f.write("alpha0_max: " + str(ALPHA0_MAX) + "\n")
        f.write("lambda_max: " + str(LAMBDA_MAX) + "\n")
        f.write("gamma_max: " + str(GAMMA_MAX) + "\n")
        f.write("seq_err: " + str(SEQ_ERR) + "\n")
        f.write("var_cp_prob: " + str(VAR_CP_PROB) + "\n")
        f.write("sc_dropout_alpha0: " + str(SC_DROPOUT_ALPHA0) + "\n")
        f.write("sc_dropout_beta0: " + str(SC_DROPOUT_BETA0) + "\n")
        f.close()

    run_command = "sbatch "
    run_command += PATH_TO_EXECUTABLE + " "
    run_command += config_file_path
    print(run_command)
    #os.system(run_command)
