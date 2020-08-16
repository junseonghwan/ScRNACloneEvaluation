# script to process PhyloWGS outputs.
import argparse
import sys
import os
import numpy as np

PATH_TO_EXECUTABLE  =  "ExtractPhyloWGSResultsLocally.sh" # path to executable
RSCRIPT_CMD = "Rscript --vanilla Rscripts/ExtractPhyloWGSResults.R"

parser = argparse.ArgumentParser(description='Script to batch process PhyloWGS results.')
parser.add_argument('case_path', help="Path to simulation case.", type=str)
parser.add_argument('--sbegin', help="Sim begin.", type=int, default = 0)
parser.add_argument('--send', help="Sim end (not inclusive).", type=int, default = 10)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 10)
parser.add_argument('--extract', help="Extract PhyloWGS outputs.", type=bool, default = False)
args = parser.parse_args()

SIMUL_PATH = args.case_path
SIM_BEGIN = args.sbegin
SIM_END = args.send
REP_BEGIN = args.rbegin
REP_END = args.rend
for sim_no in range(SIM_BEGIN, SIM_END):
	sim_path = os.path.join(SIMUL_PATH, "sim" + str(sim_no))
	for rep_no in range(REP_BEGIN, REP_END):
		rep_path = os.path.join(sim_path, "rep" + str(rep_no))

		if args.extract:
			run_command = "./"
			run_command += PATH_TO_EXECUTABLE + " "
			run_command += rep_path
			os.system(run_command)
		run_command = RSCRIPT_CMD + " " + rep_path
		print(run_command)
		os.system(run_command)


