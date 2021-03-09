#!/bin/bash
#SBATCH -J concatVCF
#SBATCH -e logs/slurm-%j.err
#SBATCH -o logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir 
# sbatch run_concatVCFs.sh spp_name

module load bcftools/1.5-fasrc02

bcftools concat ../$1/vcf/*.gz -O z -o $1_combined.vcf.gz -a 
bcftools index $1_combined.vcf.gz -t 
