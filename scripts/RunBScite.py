import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "../B-SCITE/src/bscite.exe" # path to B-SCITE executable
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('case_path', help="Path to simulation case.", type=str)
parser.add_argument('--sbegin', help="Sim begin.", type=int, default = 0)
parser.add_argument('--send', help="Sim end (not inclusive).", type=int, default = 1)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 20)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 100000)
parser.add_argument('-o', '--output_prefix', help="Output prefix.", type=str, default = "bscite")
args = parser.parse_args()

OUTPUT_FILES_PREFIX = args.output_prefix # this is the prefix used for the three output files (.matrices, .gv and .newick) given as the output of B-SCITE
CASE_PATH = args.case_path
SIM_BEGIN = args.sbegin
SIM_END = args.send
REP_BEGIN = args.rbegin
REP_END = args.rend

fp = 0.001 	# esimated false positive rate of SCS experiment
fn = 0.2     # estimated false negative rate of SCS experiment
n  = 100      # number of mutations
m  = 100      # number of cells
r  = 1       # number of repeats
l  = args.mcmc    # number of loops

for sim_no in range(SIM_BEGIN, SIM_END):
	SIMUL_PATH = CASE_PATH + "/sim" + str(sim_no)
	for rep_no in range(REP_BEGIN, REP_END):
		REP_PATH = SIMUL_PATH + "/rep" + str(rep_no) + "/"
		OUTPUT_PATH = REP_PATH + "/bscite"
		if not os.path.exists(OUTPUT_PATH):
			os.makedirs(OUTPUT_PATH)
		sc_file = REP_PATH + "/simul_bscite.SC"
		bulk_file = REP_PATH + "/simul_bscite.bulk"

		# Run Rscript to generate input for B-SCITE.
		run_command = "Rscript --vanilla Rscripts/GenerateInputForBScite.R "
		run_command += REP_PATH
		os.system(run_command)
		#print(run_command)

		# Run B-SCITE.
		run_command = PATH_TO_EXECUTABLE + " "
		run_command += "-i "   + sc_file   + " "
		run_command += "-bulkFileLocation " + bulk_file + " "
		run_command += "-n "  + str(n)  + " "
		run_command += "-m "  + str(m)  + " "
		run_command += "-fd " + str(fp) + " "
		run_command += "-ad " + str(fn) + " "
		run_command += "-r "  + str(r)  + " "
		run_command += "-l "  + str(l)  + " "
		run_command += "-o "  + OUTPUT_PATH + "/" + OUTPUT_FILES_PREFIX + " "
		os.system(run_command)
		#print(run_command)

