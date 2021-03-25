genomes = []
with open("samples.txt", "r") as f:
    for line in f:
        line = line.strip()
        genomes.append(line)
samples = glob_wildcards(expand("{directory}/{sample}.fna", directory=genomes, allow_missing=True))
print(samples)
"""rule all:
    input: expand("{directory}/genmap_output/{sample}_genomic.genmap.sorted.bedgraph", directory=genomes, allow_missing=True)

rule genmap_index:
    input: 
        fasta = "{directory}/{sample}.fna",
    log: "{directory}/{sample}_genmap_logs/index_log.txt"
    output: directory("{directory}/{sample}_genmap_index")
    shell: "genmap index -F {input.fasta} -I {output} &> {log}"

rule genmap_map:
    input: 
        index = "{directory}/{sample}_genmap_index",
    #log: "{directory}/genmap_logs/map_log.txt"
    params: outdir = directory("{directory}/genmap_output") 
    output: bg = temp("{directory}/genmap_output/{sample}_genomic.genmap.bedgraph")     
    shell: "genmap map -K 150 -E 0 -I {input.index} -O {params.outdir} -bg -T 2 -v  > {log}"

rule sort:
    input: "{directory}/genmap_output/{sample}_genomic.genmap.bedgraph"
    output: "{directory}/genmap_output/{sample}_genomic.genmap.sorted.bedgraph"
    shell: "sort -k1,1 -k2,2n {input} > {output}" """