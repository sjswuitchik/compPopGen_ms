#!/bin/bash
#SBATCH -J bwSums
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=12000

# nb: mem=8000 is fine for most, some spp needed 10000 - 12000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch run_bwSummaries.sh spp_name

set -o errexit

cd $1/

awk 'BEGIN{FS=OFS="\t"}{print $1, 0, $2, $1}' $1.chrom.sizes > $1.genome.bed

.././bedGraphToBigWig $1.merge.bg $1.chrom.sizes $1.merge.bw

.././bigWigAverageOverBed $1.merge.bw $1.genome.bed $1.summary.tab

gzip $1.merge.bg
