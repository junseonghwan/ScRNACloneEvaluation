import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "RunCanopy.sh" 
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('case_path', help="Path to simulation case.", type=str)
parser.add_argument('--sbegin', help="Sim begin.", type=int, default = 0)
parser.add_argument('--send', help="Sim end (not inclusive).", type=int, default = 10)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 10)
parser.add_argument('-kbegin', , help="Minimum number of clones.", type=int, default = 3)
parser.add_argument('-kend', , help="Minimum number of clones.", type=int, default = 12)
args = parser.parse_args()

CASE_PATH = args.case_path
SIM_BEGIN = args.sbegin
SIM_END = args.send
REP_BEGIN = args.rbegin
REP_END = args.rend

for sim_no in range(SIM_BEGIN, SIM_END):
	SIMUL_PATH = CASE_PATH + "/sim" + str(sim_no)
	for rep_no in range(REP_BEGIN, REP_END):
		REP_PATH = SIMUL_PATH + "/rep" + str(rep_no) + "/"
		rep_seed = np.random.randint(100000000)

		run_command = "sbatch "
		run_command += PATH_TO_EXECUTABLE + " "
		run_command += str(rep_seed) + " "
		run_command += REP_PATH + " "
		run_command += str(args.kbegin) + " "
		run_command += str(args.kend)
		#os.system(run_command)
		print(run_command)
