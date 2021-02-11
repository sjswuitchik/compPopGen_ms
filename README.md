# MK pipeline for comparative population genomics collaboration 

Pipeline to filter VCFs, annotate variants with snpEff, produce MK tables for downstream analyses, and run MK tests, SnIPRE, and calculate various statistics.

Authors: 


Sara Wuitchik (Postdoctoral associate, Boston University & Harvard University; sjswuit@g.harvard.edu)  

Allison Shultz (Assistant Curator of Ornithology, LA Natural History Museum; ashultz@nhm.org)

Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)

## Configuration and set up

First, set up a conda environment that will allow access to Snakemake, Python/R packages, java, and required command line tools:

```conda create -n mk -c bioconda snakemake cyvcf2 tqdm bcftools vcftools htslib java-jdk bedtools r-base r-tidyverse r-rjags r-r2jags r-lme4 r-arm```  

Activate the environment so you have access to java in order to build the snpEff database:

```conda activate mk```

### SnpEff

We use SnpEff (http://snpeff.sourceforge.net/download.html) to build databases and annotate the variants in the VCFs. It should be downloaded in your project directory and set up prior to running the pipeline.

```wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip```

```unzip snpEff_latest_core.zip```

```rm snpEff_latest_core.zip``` 

```mkdir -p snpEff/data/ingroup_species_name/```  

Note: the ingroup_species_name can be any code or form of the species identifier you'd like, so long as it's consistent throughout.  

Ensure reference sequence (FASTA) and genome annotation (GFF3) are in the appropriate data directory, rename files to sequences.fa and genes.gff, then gzip.

#### Add genome information to config file

Add the following to the snpEff.config file, under the Databases & Genomes - Non-standard Genomes section:

\# Common name genome, Source and Version

ingroup_species_name.genome : genome_name

For example: 

\# Black-headed duck genome, NCBI version 1

hetAtr.genome : Heteronetta_atricapilla


#### Build a snpEff database

From the snpEff directory, run: 

```java -jar snpEff.jar build -gff3 -v ingroup_species_name```  

For example:  

```java -jar snpEff.jar build -gff3 hetAtr```

### In your working directory, you'll need: 

- single VCF for each of the ingroup and outgroup species  

- Missingness data for both ingroup and outgroup (\*_missing_data.txt)

- Coverage site data for both ingroup and outgroup (\*_coverage_sites_clean_merged.bed)

- genes.gff (same file that's in the snpEff data directory, uncompressed)
