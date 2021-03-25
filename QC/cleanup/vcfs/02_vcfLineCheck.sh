#!/bin/bash

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/vcfs
# ./vcfLineCheck.sh spp_name

zgrep -v '^#' $1_combined.vcf.gz | wc -l > $1.log
zgrep -v '^#' ../../$1/vcf/*.vcf.gz | wc -l >> $1_split.log

## NB: not needed when refactoring to snakemake
