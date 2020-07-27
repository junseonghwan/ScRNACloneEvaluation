# script to launch PhyloWGS experiments
import sys
import os
import numpy as np

PATH_TO_EXECUTABLE  =  "RunPhyloWGS.sh" # path to executable
SEED = int(sys.argv[1])
SIMUL_PATH = sys.argv[2]
REP_BEGIN = int(sys.argv[3])
REP_END = int(sys.argv[4])
N_MCMC_ITER = int(sys.argv[5])
N_MH_ITER = int(sys.argv[6])
BURN_IN = int(sys.argv[7])

np.random.seed(SEED)

for rep_no in range(REP_BEGIN, REP_END):
    REP_PATH = SIMUL_PATH + "/rep" + str(rep_no) + "/"
    rep_seed = np.random.randint(100000000)
    rep_pwgs_bulk_path = REP_PATH + "pwgs_snv.txt"
    rep_pwgs_cnv_path = REP_PATH + "pwgs_cnv.txt"
    rep_output_path = REP_PATH + "phylowgs/"

    if not os.path.exists(rep_output_path):
        os.makedirs(rep_output_path)

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
