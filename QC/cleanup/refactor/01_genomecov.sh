#!/bin/bash
#SBATCH -J genomecov
#SBATCH -e logs/slurm-%j.err
#SBATCH -o logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 05-00:00:00
#SBATCH --mem=10000

# sbatch 01_genomecov.sh spp_name

set -o errexit

source activate qc

for file in $1/*.bam;
do
  bedtools genomecov -bga -ibam $file -g $1/$1.chrom.sizes > $file.bg
done

for file in $1/*.bg;
do
  sort -k1,1 -k2,2n $file > $file.sorted
done

rename 's/\.dedup\.sorted\.bam\.bg\.sorted/\.bg/' $1/*.sorted 
mkdir -p $1/unsortedBG
mv $1/*.dedup.sorted.bam.bg $1/unsortedBG/

for file in $1/*.bg;
do
  bedGraphToBigWig $file $1/$1.chrom.sizes $file.bw
done

gzip $1/*.bg
