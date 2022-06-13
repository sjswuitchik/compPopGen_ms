localrules: snpEffdb
  
rule download_reference:
 """
 This rule downloads the NCBI dataset for the reference genome and annotation to build the snpEff database with
 """
    output:
        outdir = directory(config["refGenomeDir"] + "{refGenome}"),
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    params:
        dataset = config["refGenomeDir"] + "{refGenome}_dataset.zip"
    log:
        "logs/dl_reference/{refGenome}_snpEff.log"
    conda:
        "../envs/ncbi.yml"
    shell:
        "datasets download genome accession --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} &> {log}"
        "&& 7z x {params.dataset} -aoa -o{output.outdir}"
        "&& cat {output.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}"
  
rule reorganize:
  """
  This rule organizes & renames the data for the snpEff database creation
  """
  input:
    ref = config["refGenomeDir"] + "{refGenome}.fna",
    gff = config["refGenomeDir"] + "{refGenome}.gff"
  output:
    ref = config["snpEffDir"] + "data/" + config["ingroup"] + "/sequences.fa",
    gff = config["snpEffDir"] + "data/" + config["ingroup"] + "/genes.gff"
  shell:
    "mv {input.ref} {output.ref}\n"
    "mv {input.gff} {output.gff}"
    
rule compress:
  """
  This rule gzips the reference genome and annotation for the snpEff database build
  """
  input:
    ref = config["snpEffDir"] + "data/" + config["ingroup"] + "/sequences.fa",
    gff = config["snpEffDir"] + "data/" + config["ingroup"] + "/genes.gff"
  output:
    ref = config["snpEffDir"] + "data/" + config["ingroup"] + "/sequences.fa.gz",
    gff = config["snpEffDir"] + "data/" + config["ingroup"] + "/genes.gff.gz"
  shell:
    "gzip {input.ref}\n"
    "gzip {input.gff}"
