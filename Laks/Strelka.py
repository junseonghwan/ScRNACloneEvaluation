#!/usr/bin/env python

import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

np.random.seed(1)

# code to process command line arguments
parser = argparse.ArgumentParser(description='Script to call variants using Strelka on pseudo-BAMs.')
parser.add_argument('sample_path', help="Path to the pseudo-bulk.", type=str)
parser.add_argument('normal_bam', help="Path to the normal BAM file.", type=str)
args = parser.parse_args()

command = "sbatch /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/Laks/MergeBAMs.sh "
command += os.path.join(args.sample_path, "sample.sorted.bam") + " "
command += args.normal_bam + " "
command += args.sample_path
print(command)
#os.system(command)
