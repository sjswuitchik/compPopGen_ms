import glob
import os
directories = []
with open("samples.txt", "r") as f:
    for line in f:
        line = line.strip()
        directories.append(line)
samples = []
for d in directories:
    samp = os.path.basename(glob.glob(f"{d}/*.fna")[0])
    samples.append(samp.split(".")[0])
rule all:
    input: expand("{directory}/{sample}_genmap_output/{sample}.genmap.sorted.bedgraph", zip, directory=directories, sample=samples)

rule genmap_index:
    input: 
        fasta = "{directory}/{sample}.fna",
    log: "{directory}/{sample}_genmap_logs/index_log.txt"
    output: directory("{directory}/{sample}_genmap_index")
    shell: "genmap index -F {input.fasta} -I {output} &> {log}"

rule genmap_map:
    input: 
        index = "{directory}/{sample}_genmap_index",
    log: "{directory}/{sample}_genmap_logs/map_log.txt"
    params: outdir = directory("{directory}/{sample}_genmap_output") 
    output: bg = temp("{directory}/{sample}_genmap_output/{sample}.genmap.bedgraph")     
    shell: "genmap map -K 150 -E 0 -I {input.index} -O {params.outdir} -bg -T 2 -v  > {log}"

rule sort:
    input: "{directory}/{sample}_genmap_output/{sample}.genmap.bedgraph"
    output: "{directory}/{sample}_genmap_output/{sample}.genmap.sorted.bedgraph"
    shell: "sort -k1,1 -k2,2n {input} > {output}" 