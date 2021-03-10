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
  bedtools genomecov -bga -ibam $file -g /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/$1/genome/$1.fa > coverage/$file.statscov.bg
  bedtools genomecov $file -g /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/$1/genome/$1.fa > coverage/$file.hist
done
