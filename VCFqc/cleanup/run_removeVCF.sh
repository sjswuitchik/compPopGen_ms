#!/bin/bash
#SBATCH -J vcfRemove
#SBATCH -o gatherVCFs_dir/logs/out_%j
#SBATCH -e gatherVCFs_dir/logs/err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=4000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS
# sbatch run_removeVCF.sh spp_name

rm -r $1/vcf/*
