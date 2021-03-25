#!/bin/bash
#SBATCH -J gatkUpdate
#SBATCH -o logs/out_%j
#SBATCH -e logs/err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=12000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/vcfs
# sbatch run_gatkUpdate.sh spp_name

set -o errexit

source activate gatk

bgzip -dc $1_combined.vcf.gz | sed -e 's/\<nan\>/NaN/g' | bgzip > $1_combined.rep.vcf.gz 

bcftools index -t $1_combined.rep.vcf.gz

picard SortVcf -Xmx8g --TMP_DIR=$PWD/tmp -I $1_combined.rep.vcf.gz -O $1_comboSorted.vcf.gz

gatk IndexFeatureFile -I $1_comboSorted.vcf.gz

gatk VariantFiltration -R /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/$1/genome/$1.fa -V $1_comboSorted.vcf.gz --output $1_updatedFilter.vcf.gz --filter-name "RPRS_filter" --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" --filter-name "FS_SOR_filter" --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" --filter-name "MQ_filter" --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))"


## NB: not needed when refactoring to snakemake (bc filter description should already be updated in main pipeline)
