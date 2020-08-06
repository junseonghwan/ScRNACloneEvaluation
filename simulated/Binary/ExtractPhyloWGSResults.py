# script to process PhyloWGS outputs.
import sys
import os
import numpy as np

PATH_TO_EXECUTABLE  =  "ExtractPhyloWGSResultsLocally.sh" # path to executable
SIMUL_PATH = sys.argv[1]
SIM_BEGIN = int(sys.argv[2])
SIM_END = int(sys.argv[3])
REP_BEGIN = int(sys.argv[4])
REP_END = int(sys.argv[5])

for sim_no in range(SIM_BEGIN, SIM_END):
	sim_path = os.path.join(SIMUL_PATH, "sim" + str(sim_no))
	for rep_no in range(REP_BEGIN, REP_END):
	    rep_path = os.path.join(sim_path, "rep" + str(rep_no))

	    #run_command = "sbatch "
	    run_command = "./"
	    run_command += PATH_TO_EXECUTABLE + " "
	    run_command += rep_path

	    #print(run_command)
	    os.system(run_command)
