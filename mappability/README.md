
# Snakemake Pipeline for Computing Mappability
This small Snakemake pipeline makes use of Genmap (https://github.com/cpockrandt/genmap) to compute the mappability of an input genome fasta file.

# Setup and Usage  
#### note from Sara: Dependabot gave three high severity alerts regarding requirements.txt. For now, created env.yml and commented out content of requirements.txt
First, setup your conda environment:

``` conda create --name <env> --file <requirements> ```

And then activate it:

``` conda activate <env> ```

The `samples.txt`file should contain the names (1 per line) of the names of the directories containing the fasta (.fna) file you would like to process.

In your working directory you should have:
- `mappability.smk`
- `samples.txt`
- `Directories containing fasta file` 

Execute the pipeline using snakemake:

`snakemake -s mappability.smk --cores <n>`
