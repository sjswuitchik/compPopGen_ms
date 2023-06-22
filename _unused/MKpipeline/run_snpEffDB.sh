#!/bin/bash
#SBATCH -J sm_db
#SBATCH -o out_db
#SBATCH -e err_db
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

source activate snakemake

snakemake --snakefile workflow/Snakefile_snpEffDB --profile ./profiles/slurm
