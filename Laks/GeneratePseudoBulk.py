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
parser = argparse.ArgumentParser(description='Script to generate pseudo bulk file from single cell BAMs.')
parser.add_argument('cell_to_clone_file', help="Path to cell to clone membership file.", type=str)
parser.add_argument('bam_paths', help="Path to the BAM files to be merged separated by comma.", type=str)
parser.add_argument('output_path', help="Path to output the merged BAM file.", type=str)
parser.add_argument('bam_out_name', help="BAM output file name prefix.", type=str)
parser.add_argument('num_samples', help="Number of samples.", type=int)
parser.add_argument('num_cells', help="Number of cells per region.", type=int)
parser.add_argument('-d', '--dirichlet_param', help="Dirichlet parameter.", type=float, default=1.0)
parser.add_argument('-s', '--suffix', help="BAM file suffix.", type=str, default=".bam")
args = parser.parse_args()

cell2clone = pd.read_csv(args.cell_to_clone_file)

bam_paths = args.bam_paths.split(",") # this is a list.
bam_files = []
suffix = args.suffix
for i in range(len(bam_paths)):
	bam_path = bam_paths[i]
	bam_files += glob.glob(os.path.join(bam_path) + "/*" + suffix)
cell_ids = list(map(lambda bam_file : os.path.basename(bam_file), bam_files))
sc_bams = pd.DataFrame(data = {'bam_path': bam_files, 'cell_id': cell_ids})
joined = sc_bams.join(cell2clone, lsuffix='', rsuffix='_other')
joined = joined[~joined["clone_id"].isin(["G","H","I",np.nan])]

clones = list(joined["clone_id"].unique())
num_clones = len(clones)
num_samples = args.num_samples
print(clones)
print(num_clones)
print(num_samples)
# Sample proportions vector for each clone from Dirichlet.
proportions = np.random.dirichlet(np.repeat(args.dirichlet_param, num_clones), num_samples)
for sample in range(num_samples):
	sample_outpath = os.path.join(args.output_path, "Sample" + str(sample))
	if not os.path.exists(args.output_path):
		os.makedirs(args.output_path)

	bam_list_file = os.path.join(sample_outpath, args.bam_out_name + ".txt")
	bam_out = os.path.join(sample_outpath, args.bam_out_name + ".bam")

	sample_props = proportions[sample]
	sampled_idxs = []
	lines = ""
	print("Sample " + str(sample))
	print(sample_props)
	for i in range(num_clones):
		n_files = int(round(sample_props[i] * args.num_cells))
		cell_subset = joined[joined["clone_id"] == clones[i]]
		print(str(n_files) + "/" + str(len(cell_subset)))
		sampled_idxs = np.random.choice(len(cell_subset), n_files, replace=True)

		# Sample bam files to include in forming the pseudo-bulk.	
		lines += "\n".join(list(cell_subset.iloc[sampled_idxs]["bam_path"]))
		lines += "\n"
	with open(bam_list_file, "w+") as f:
		f.write(lines)
	f.close()

	command = "sbatch /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/Laks/MergeBAMs.sh " 
	command += bam_list_file + " "
	command += bam_out
	#print(command)
	os.system(command)

