#!/bin/bash

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/vcfs

# ./vcfHeaderClean.sh spp_name

module load htslib/1.5-fasrc02 bcftools/1.5-fasrc02

zcat $1_updatedFilter.vcf.gz | \
sed -e 's/##FILTER=<ID=LowQual,Description="Low quality">/#/g' \
-e 's/##FILTER=<ID=GATK_default/#/g' \
-e $'s/##FILTER=<ID=FS_SOR_filter,Description="(vc.isSNP() && ((vc.hasAttribute(\'FS\') && FS > 60.0) || (vc.hasAttribute(\'SOR\') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute(\'FS\') && FS > 200.0) || (vc.hasAttribute(\'SOR\') &&  SOR > 10.0)))">/##FS_SOR_filter,Description="filter SNPs for strand bias with Phred-scaled p-value for Fishers exact test above 60 and Symmetric Odds Ratio above 3; or indels; or if mixed, a Phred-scaled p-value above 200 and a Symmetric Odds Ratio above 10"/g' \
-e $'s/##FILTER=<ID=MQ_filter,Description="(vc.isSNP() && ((vc.hasAttribute(\'MQ\') && MQ < 40.0) || (vc.hasAttribute(\'MQRankSum\') && MQRankSum < -12.5))">/##MQ_filter,Description="filter SNPs with RMS mapping quality less than 40 and Z-score for Wilcoxon rank sum test for read mapping quality less than -12.5"/g' \
-e $'s/##FILTER=<ID=RPRS_filter,Description="(vc.isSNP() && (vc.hasAttribute(\'ReadPosRankSum\') && ReadPosRankSum < -8.0))|| ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute(\'ReadPosRankSum\') && ReadPosRankSum < -20.0)) || (vc.hasAttribute(\'QD\') && QD < 2.0)">/##RPRS_filter,Description="filter SNPs with Z-Score for Wilcoxon rank sum test for read position bias less than -8; or indels; or if mixed, a Z-Score for Wilcoxon rank sum test less than -20"/g' | \
bgzip -c > $1_final.clean.vcf.gz

bcftools index -t $1_final.clean.vcf.gz

