#!/bin/bash
#SBATCH -J bgMove
#SBATCH -o gatherVCFs_dir/coverage/logs/out_%j
#SBATCH -e gatherVCFs_dir/coverage/logs/err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=4000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS
# sbatch moveBedgraphs.sh spp_name

cd $1dedup
mv *.bg /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/$1
