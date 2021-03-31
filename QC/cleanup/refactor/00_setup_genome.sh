#!/bin/bash
#SBATCH -J genome_set
#SBATCH -o logs/out_%j.out
#SBATCH -e logs/err_%j.err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=4000

# sbatch 00_setup_genome.sh spp_name

faToTwoBit -long $1/$1.fa $1/$1.2bit
twoBitInfo $1/$1.2bit stdout | sort -k2rn > $1/$1.chrom.sizes
awk 'BEGIN{FS=OFS="\t"}{print $1, 0, $2, $1}' $1/$1.chrom.sizes > $1/$1.genome.bed
