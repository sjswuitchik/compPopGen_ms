import csv

"""
Snakemake rules for downloading fastqs from SRA and gzipping them. Reads from tsv to determine SRRs to download. Change path to tsv file below. 

When executing snakemake, make sure to include the option '--use-conda'
"""

TSV_FILE = ".tsv"

accessions = []

with open(TSV_FILE, "r") as f:
    rd = csv.reader(f, delimiter="\t")
    for row in rd:
        accessions.append(row[0])

rule all:
    input: expand("fastqs/{accession}_{num}.fastq.gz", accession=accessions, num=[1,2])
rule get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "fastqs/{accession}_1.fastq",
        "fastqs/{accession}_2.fastq"
    params:
        # optional extra arguments
        extra=""
    threads: 6  # defaults to 6
    wrapper:
        "0.73.0/bio/sra-tools/fasterq-dump"
    
rule gzip_fastq:
    input:
        "fastqs/{accession}_1.fastq",
        "fastqs/{accession}_2.fastq"
    output:
        "fastqs/{accession}_1.fastq.gz",
        "fastqs/{accession}_2.fastq.gz"
    shell:
        "gzip {input}"



