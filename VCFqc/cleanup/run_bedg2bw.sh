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

set -o errexit

for file in $1/*.bg;
do
  sort -k1,1 -k2,2n $file > $file.sorted
  ./brename -p ".dedup.sorted.bam.bg.sorted" -r ".bg" -R
done

mkdir -p $1/unsortedBG
mv $1/*.dedup.sorted.bam.bg $1/unsortedBG

for file in $1/*.bg;
do
  ./bedGraphToBigWig $file $1/$1.chrom.sizes $file.bw
done

cd $1/ 
ls *.bw > list
.././bigWigMerge -inList list $1.merge.bg
