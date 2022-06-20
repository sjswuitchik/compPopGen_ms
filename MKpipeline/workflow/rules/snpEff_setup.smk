import glob
import re
import os

import pandas as pd
from collections import defaultdict, deque
from snakemake.exceptions import WorkflowError

### INPUT FUNCTION ###
def get_ref(wildcards):
    if 'refPath' in samples.columns:
        _refs = samples.loc[(samples['refGenome'] == wildcards.refGenome)]['refPath'].dropna().unique().tolist()
        for ref in _refs:
            print(ref)
            if not os.path.exists(ref):
                raise WorkflowError(f"Reference genome {ref} does not exist")
            elif ref.rsplit(".", 1)[1] == '.gz':
                raise WorkflowError(f"Reference genome {ref} must be unzipped first.")
        return _refs
    else:
        return []
    
### RULES ### 

rule download_reference:
 """
 This rule downloads the NCBI dataset for the reference genome and annotation to build the snpEff database with
 """
    input:
        ref = get_ref
    output:
        ref = directory(config['refGenomeDir']) + "{refGenome}/{refGenome}.fna",
        gff = directory(config['refGenomeDir']) + "{refGenome}/genomic.gff"
    params:
        dataset = directory(config['refGenomeDir']) + "{refGenome}/{refGenome}_dataset.zip",
        outdir = directory(config['refGenomeDir']) + "{refGenome}"
    conda:
        "../envs/ncbi.yml"
    shell:
        "mkdir -p {params.outdir}\n"
        "datasets download genome accession {wildcards.refGenome} --exclude-protein --exclude-rna --filename {params.dataset}\n"
        "7z x {params.dataset} -aoa -o{params.outdir}\n"
        "cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}\n"
        "mv {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/genomic.gff {output.gff}"

rule reorganize:
  """
  This rule organizes & renames the data for the snpEff database creation, and copies the GFF to the correct location for MK and SnIPRE
  """
  input:
    seq = expand(directory(config['refGenomeDir']) + "{refGenome}/{refGenome}.fna", refGenome=REFGENOME),
    genes = expand(directory(config['refGenomeDir']) + "{refGenome}/genomic.gff", refGenome=REFGENOME)
  output:
    ref = "data/" + directory(config["ingroup"]) + "/sequences.fa",
    gff = "data/" + directory(config["ingroup"]) + "/genes.gff"
  params:
    ingroup = config['ingroup']
  shell:
    "mkdir -p snpEff/data/{params.ingroup}\n"
    "cp {input.seq} {output.ref}\n"
    "cp {input.genes} {output.gff}"
    "mkdir -p data/mk_tests/\n"
    "cp {output.gff} data/mk_tests/"
    
rule compress:
  """
  This rule gzips the reference genome and annotation for the snpEff database build
  """
  input:
    ref = "data/" + directory(config["ingroup"]) + "/sequences.fa",
    gff = "data/" + directory(config["ingroup"]) + "/genes.gff"
  output:
    ref = "data/" + directory(config["ingroup"]) + "/sequences.fa.gz",
    gff = "data/" + directory(config["ingroup"]) + "/genes.gff.gz"
  shell:
    "gzip {input.ref}\n"
    "gzip {input.gff}"
   
   
