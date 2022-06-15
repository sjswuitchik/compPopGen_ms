# Automated pipeline for MK tests, SnIPRE, and comparative population genomics stats

Pipeline to filter VCFs, annotate variants with snpEff, produce MK tables for downstream analyses, and run MK tests, SnIPRE, and calculate various statistics.

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

We use SnpEff (http://snpeff.sourceforge.net/download.html) to build databases and annotate the variants in the VCFs. It should be downloaded in your project directory and set up prior to running the pipeline.

```wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip```

```unzip snpEff_latest_core.zip```

```rm snpEff_latest_core.zip``` 

```mkdir -p snpEff/data/ingroup_species_name/```  

Note: the ingroup_species_name can be any code or form of the species identifier you'd like, so long as it's consistent throughout.  

The first Snakefile will download and organize the files required to build the snpEff annotation database, but you will need to manually edit the snpEff config before building the database.

#### Add genome information to config file

Add the following to the snpEff.config file, under the Databases & Genomes - Non-standard Genomes section:

\# Common name genome, Source and Version

ingroup_species_name.genome : genome_name

For example: 

\# Black-headed duck genome, NCBI version 1

hetAtr.genome : Heteronetta_atricapilla  

Note: you can also put the data directory into the config to avoid confusion between your working directory and your snpEff directory.  


#### Build a snpEff database

From the snpEff directory with the conda environment activated, run: 

```snpEff build -c snpEff.config -gff3 -v -noCheckCds -noCheckProtein ingroup_species_name```  

For example:  

```snpEff -Xmx8g build -c snpEff.config -gff3 -v -noCheckCds -noCheckProtein hetAtr```  

Note: the ```-Xmx8g``` flag may be required if you are getting Java heap space OOM errors  

### In your working directory, you'll need: 

- single VCF for each of the ingroup and outgroup species  

- Missingness data for both ingroup and outgroup (\*_missing_data.txt)

- Coverage site data for both ingroup and outgroup (\*_coverage_sites_clean_merged.bed)

- genes.gff (same file that's in the snpEff data directory, uncompressed)  

- snpEff
