# Reproducible workflow for comparative population genomics
Repository for manuscript

MKpipeline: snakemake workflow for running McDonald-Kreitman (MK) tests and SnIPRE.  

QC: currently standalone scripts to standardize output from two similar fastq2vcf workflows. Next step: snakemake workflow of basic coverage calculations and VCF QC to be run between fastq2vcf and MKpipeline.  

SRA: Entrez searchs and SRA run selector to create sample metadata.  

Mappability: snakemake workflow to compute mappability to be included in MKpipeline workflow.  
