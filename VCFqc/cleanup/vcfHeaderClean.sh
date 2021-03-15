#!/bin/bash

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/vcfs

# ./vcfHeaderClean.sh spp_name

zcat $1_combined.vcf.gz | sed 's/##FILTER=<ID=LowQual,Description="Low quality">/#/g' | sed 's/##FILTER=<ID=GATK_default/#/g'  > $1_comboClean.vcf.gz

