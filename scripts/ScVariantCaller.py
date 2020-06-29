#!/usr/bin/env python

import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

# code to process command line arguments
parser = argparse.ArgumentParser(description='Script to call variants from scRNA BAM file.')
parser.add_argument('bam_path', help="Path to the BAM files.", type=str)
parser.add_argument('output_name', help="Directory name for output. The read counts will go to bam_path/output_name", type=str)
parser.add_argument('-f', '--file_pattern_suffix', help="BAM file suffix.", type=str, default=".bam")
args = parser.parse_args()

bam_path = os.path.join(args.bam_path)
output_path = os.path.join(bam_path, args.output_name)
if not os.path.exists(output_path):
	os.makedirs(output_path)
samples = glob.glob(args.bam_path + "/*" + args.file_pattern_suffix)
for sample in samples:
	sample_name = os.path.basename(sample)
	cell_name = sample_name.split(args.file_pattern_suffix)[0]
	output_name = os.path.join(output_path, cell_name)
	command = "sbatch ScVariantCaller.sh " + sample + " " + output_name
	#os.system(command)
	print(command)
