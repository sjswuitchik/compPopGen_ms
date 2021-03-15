#!/bin/bash
#SBATCH -J bedg2bw
#SBATCH -e logs/slurm-%j.err
#SBATCH -o logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 05-00:00:00
#SBATCH --mem=10000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch run_bedg2bw.sh spp_name

for file in $1/*.bg;
do
  sort -k1,1 -k2,2n $file > $file.sorted.bg
  ./brename -p ".dedup.sorted.bam.bg.sorted.bg" -r ".bg" -R
  ./bedGraphToBigWig $file.sorted.bg $1.chrom.sizes $file.bw
done
