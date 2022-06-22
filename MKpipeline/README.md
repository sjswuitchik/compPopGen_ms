# Automated pipeline for MK tests, SnIPRE, and comparative population genomics stats

Pipeline to filter VCFs, annotate variants with snpEff, produce MK tables for downstream analyses, then run MK tests, SnIPRE, and calculate various statistics.

Authors: 


Sara Smith Wuitchik (Postdoctoral fellow, University of Victoria & Simon Fraser University; sjswuit@g.harvard.edu)  

Allison Shultz (Assistant Curator of Ornithology, LA Natural History Museum; ashultz@nhm.org)

Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)

## Configuration and set up

Create a conda environment that will allow access to Snakemake:

```conda install -n base -c conda-forge mamba```

```conda activate base```

```mamba create -c conda-forge -c bionconda -n snakemake snakemake mamba```

Activate the environment:

```conda activate snakemake```  

The exact name of the environment (```snakemake```) is necessary for the run scripts. If you would like to change the name of the environment, make sure to edit the run*.sh files to reflect the name change.

### SnpEff

We use SnpEff (http://snpeff.sourceforge.net/download.html) to build databases and annotate the variants in the VCFs. The first Snakefile will download and organize the files required to build the snpEff annotation database, but you will need to manually edit the snpEff config before running the second Snakefile.

#### Add genome information to config file

Add the following to the snpEff.config file, under the Databases & Genomes - Non-standard Genomes section:

\# Common name genome, Source and Version

ingroup_species_name.genome : genome_name

For example: 

\# Black-headed duck genome, NCBI version 1

hetAtr.genome : Heteronetta_atricapilla  

### In the parent directory (MKpipeline), you'll need a directory called input_data that contains:

- single VCF for each of the ingroup and outgroup species and the associated index (\*.vcf.gz and \*.vcf.gz.tbi) e.g., hetAtr.vcf.gz and hetAtr.vcf.gz.tbi 

- a list of individuals to remove for each of the ingroup and outgroup species (\*.remove.indv) e.g., hetAtr.remove.indv

- a callable sites BED for each of the ingroup and outgroup (\*\.callable_sites.bed) e.g., hetAtr.callable_sites.bed
