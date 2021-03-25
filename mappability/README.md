
# Snakemake Pipeline for Computing Mappability
This small Snakemake pipeline makes use of Genmap (https://github.com/cpockrandt/genmap) to compute the mappability of an input genome fasta file.

# Setup and Usage
First, setup your conda environment:

``` conda create --name <env> --file <requirements> ```

And then activate it:

``` conda activate <env> ```

In your working directory you should have:
- `mappability.smk`
- `samples.txt`
- `Directories containing fasta file` 
