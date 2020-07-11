import sys
import os

PATH_TO_EXECUTABLE  =  "../B-SCITE/src/bscite.exe" # path to B-SCITE executable
OUTPUT_FILES_PREFIX = "/Users/seonghwanjun/data/cell-line/bulk/OV2295/B-SCITE/bscite" # this is the prefix used for the three output files (.matrices, .gv and .newick) given as the output of B-SCITE

SCFile   = "/Users/seonghwanjun/data/cell-line/bulk/OV2295/B-SCITE/bscite.SC " # path input SC File
bulkFile = "/Users/seonghwanjun/data/cell-line/bulk/OV2295/B-SCITE/bscite.bulk"     # path input bulk file
fp = 0.01 # esimated false positive rate of SCS experiment
fn = 0.2     # estimated false negative rate of SCS experiment
n  = 45      # number of mutations
m  = 276      # number of cells
r  = 1       # number of repeats
l  = 2000    # number of loops

run_command =""
run_command += PATH_TO_EXECUTABLE + " "
run_command += "-i "   + SCFile   + " "
run_command += "-bulkFileLocation " + bulkFile + " "
run_command += "-n "  + str(n)  + " "
run_command += "-m "  + str(m)  + " "
run_command += "-fd " + str(fp) + " "
run_command += "-ad " + str(fn) + " "
run_command += "-r "  + str(r)  + " "
run_command += "-l "  + str(l)  + " "
run_command += "-o "  + OUTPUT_FILES_PREFIX + " "
os.system(run_command)
