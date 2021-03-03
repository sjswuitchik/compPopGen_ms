#!/bin/bash
#SBATCH -J vcfqc
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

source activate vcfqc
snakemake --snakefile Snakefile_vcfqc --profile ./profiles/slurm
