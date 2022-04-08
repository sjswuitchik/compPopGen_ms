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
    ref = 
    gff = 
  output:
    ref = 
    gff = 
  shell:
    "mv {input.} {output.}\n"
    "mv {input.} {output.}"
