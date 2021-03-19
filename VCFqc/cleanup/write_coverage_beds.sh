#!/bin/bash
#SBATCH -J covBeds
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch write_coverage_beds.sh spp_name genomewideavgcov

## untested as written here 

cd $1/
gzip -dc $1.merge.bg.gz | awk -v "avg=$2" "spp=$1" -f sum_cov.awk
