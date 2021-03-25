#!/bin/bash
#SBATCH -J genomecov
#SBATCH -e logs/err_%j
#SBATCH -o logs/out_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# sbatch 01_genomecov.sh spp_name

set -o errexit

for file in $1/*.bam;
do
  bedtools genomecov -bga -ibam $file -g $1/$1.chrom.sizes > $file.bg
done

for file in $1/*.bg;
do
  sort -k1,1 -k2,2n - > $file.sorted
done

for file in $1/*.bw;
do
  bedGraphToBigWig $file $1/$1.chrom.sizes $file.bw
done

gzip $1/*.bg
