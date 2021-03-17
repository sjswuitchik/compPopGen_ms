#!/bin/bash

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/vcfs

# ./vcfHeaderClean.sh spp_name

### NB: in progress

zcat $1_updatedFilter.vcf.gz | sed -e 's/##FILTER=<ID=LowQual,Description="Low quality">/#/g' | sed -e 's/##FILTER=<ID=GATK_default/#/g' | sed -e 's/##FILTER=<ID=FS_SOR_filter,Description="(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))">/##FS_SOR_filter,Description="filter SNPs for strand bias with Phred-scaled p-value for Fisher's exact test above 60 and Symmetric Odds Ratio above 3; or indels; or if mixed, a Phred-scaled p-value above 200 and a Symmetric Odds Ratio above 10"/g' | sed -e 's/> $1_final.clean.vcf.gz

