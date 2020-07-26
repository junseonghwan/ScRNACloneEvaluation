# Scripts and R code for analyzing data from Laks

## Pre-processing
- Merge single cell DLP data to generate pseudo-bulk data.
    + `GeneratePseudoBulk.py`.
    + `MergeBAMs.sh`.
- Run copy number analysis. We use [TitanCNA](https://github.com/gavinha/TitanCNA).

## Extracting data
- Extract exonic SNVs from the SNV data [here](https://zenodo.org/record/3445364/files/ov2295_clone_snvs.csv.gz.
- Get read counts from the pseudo bulk data and generate bulk data.
    + Do this for each region separately for multi-region analysis.
- Extract read counts from single cell RNA-seq data.
- Generate single cell data and single cell hyper parameters.

## Data for methods for comparison
- PhyloWGS: 
    + Generate VCF file at exonic SNVs from the pseudo-bulk file.
    + Run PhyloWGS' parser that takes the VCF and TitanCNA output to generate the input files.
- Generate single cell data matrices for B-SCITE and ddClone.
    + Custom script.
- Generate data for Canopy:
    + Need to run a separate copy number analysis software?

## Scripts for inference
- Our method
- PhyloWGS
- B-SCITE
- ddClone
- Canopy

## Scripts for evaluation
- Code written in Python notebook.
