#!/bin/bash
#SBATCH -J sm_mk
#SBATCH -o out_mk
#SBATCH -e err_mk
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

source activate snakemake

snakemake --snakefile workflow/Snakefile_vcf2mk --profile ./profiles/slurm
