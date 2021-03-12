#!/bin/bash
#SBATCH -J bedg2bw
#SBATCH -e coverage/logs/slurm-%j.err
#SBATCH -o coverage/logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch run_bedg2bw.sh spp_name

for file in $1/*.bg;
do
  ./bedGraphToBigWig $file.bg $1.chrom.sizes $file.bw
  ./brename -p ".dedup.sorted.bam.bg.bw" -r ".bw" -R
done
