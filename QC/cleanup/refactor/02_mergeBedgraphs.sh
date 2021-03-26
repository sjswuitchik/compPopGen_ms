#!/bin/bash
#SBATCH -J mergeBedg
#SBATCH -e logs/err_%j
#SBATCH -o logs/out_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# sbatch 02_mergeBedgraphs.sh spp_name

set -o errexit

source activate qc

ls $1/*.bw > $1/list
bigWigMerge -inList $1/list $1/$1.merge.bg

bedGraphToBigWig $1/$1.merge.bg $1/$1.chrom.sizes $1/$1.merge.bw

bigWigAverageOverBed $1/$1.merge.bw $1/$1.genome.bed $1/$1.summary.tab

gzip $1/$1.merge.bg
