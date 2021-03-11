#!/bin/bash
#SBATCH -J bedgraphs
#SBATCH -e coverage/logs/slurm-%j.err
#SBATCH -o coverage/logs/slurm-%j.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir 
# sbatch run_genomecov.sh spp_name

module load bedtools2/2.26.0-fasrc01

for file in /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/$1/dedup/*.dedup.sorted.bam;
do
  coverage/./faToTwoBit /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/$1/genome/$1.fa coverage/$1.2bit
  coverage/./twoBitInfo coverage/$1.2bit stdout | sort -k2rn > coverage/$1.chrom.sizes
  bedtools genomecov -bga -ibam $file -g coverage/$1.chrom.sizes > $file.bg
  bedtools genomecov -ibam $file -g coverage/$1.chrom.sizes > $file.hist
done
