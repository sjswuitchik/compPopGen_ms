#!/bin/bash
#SBATCH -J bgMove
#SBATCH -o out_%j
#SBATCH -e err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS
# sbatch moveBedgraphs.sh spp_name

cd $1/dedup
mv *.bg /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/$1
