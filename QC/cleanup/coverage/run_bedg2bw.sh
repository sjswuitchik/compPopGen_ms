#!/bin/bash
#SBATCH -J bedg2bw
#SBATCH -e logs/slurm-%j.err
#SBATCH -o logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=10000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch run_bedg2bw.sh spp_name

set -o errexit

cd $1/

for file in *.gz;
do
  zcat $file | sort -k1,1 -k2,2n - > $file.sorted
  .././brename -p ".dedup.sorted.bam.bg.gz.sorted" -r ".bg" -R
done

mkdir -p unsortedBG
mv *.dedup.sorted.bam.bg.gz unsortedBG/

for file in *.bg;
do
  .././bedGraphToBigWig $file $1.chrom.sizes $file.bw
done
 
ls *.bw > list
.././bigWigMerge -inList list $1.merge.bg
