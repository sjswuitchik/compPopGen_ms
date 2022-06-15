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
        
ruleorder: compress > reorganize > download_reference
    
### RULES ###

rule download_reference:
 """
 This rule downloads the NCBI dataset for the reference genome and annotation to build the snpEff database with
 """
    input:
        ref = get_ref
    output:
        ref = directory(config['refGenomeDir']) + "{refGenome}.fna",
        gff = directory(config['refGenomeDir']) + "{refGenome}.gff"
    params:
        dataset = directory(config['refGenomeDir']) + "{refGenome}_dataset.zip",
        outdir = directory(config['refGenomeDir']) + "{refGenome}"
    log:
        "logs/dl_reference/{refGenome}_snpEff.log"
    conda:
        "../envs/ncbi.yml"
    shell:
        "mkdir -p {params.outdir} \n"
        "datasets download genome accession --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} &> {log} \n"
        "&& 7z x {params.dataset} -aoa -o{params.outdir} \n"
        "&& cat {output.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}"
  
rule reorganize:
  """
  This rule organizes & renames the data for the snpEff database creation
  """
  input:
    seq = directory(config['refGenomeDir']) + "{refGenome}.fna",
    genes = directory(config['refGenomeDir']) + "{refGenome}.gff"
  output:
    ref = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/sequences.fa",
    gff = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/genes.gff"
  params:
    ingroup = config['ingroup']
  shell:
    """
    mkdir -p snpEff/{wildcards.refGenome}/data/{params.ingroup} \
    mv {input.seq} {output.ref}
    my {input.genes} {output.gff}
    """
    
rule compress:
  """
  This rule gzips the reference genome and annotation for the snpEff database build
  """
  input:
    ref = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/sequences.fa",
    gff = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/genes.gff"
  output:
    ref = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/sequences.fa.gz",
    gff = directory(config["snpEffDir"]) + "data/" + directory(config["ingroup"]) + "{refGenome}/genes.gff.gz"
  shell:
    "gzip {input.ref}\n"
    "gzip {input.gff}"
   
   
