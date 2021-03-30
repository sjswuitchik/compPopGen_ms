### VCF quality control pipeline ###  

This pipeline is in development to be added on to https://github.com/harvardinformatics/shortRead_mapping_variantCalling to provide some qc before downstream analyses are run.  

This pipeline currently relies on a conda environment called vcfqc. It can be built with the following command:  

```conda create -n vcfqc -c bioconda plink vcftools r-base r-tidyverse snakemake```  

NB: Once coverage and mappability is added, this env will be updated to be built by the env.yml
