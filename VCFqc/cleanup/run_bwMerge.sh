#!/bin/bash
#SBATCH -J bwmerge
#SBATCH -e logs/slurm-%j.err
#SBATCH -o logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch run_bwMerge.sh spp_name

ls $1/*.bw > $1/list
./bigWigMerge -inList $1/list $1.merge.bg
