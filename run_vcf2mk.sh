#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module purge
module load Anaconda/5.0.1-fasrc02

source activate mk
snakemake --snakefile Snakefile_vcf2mk 


# currently hitting some issues with SLURM architecture but runs without this flag (2021/02/09)
# --profile ./profiles/slurm
