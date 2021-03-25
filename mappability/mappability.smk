genomes = []
with open("samples.txt", "r") as f:
    for line in f:
        line = line.strip()
        genomes.append(line)

rule all:
    input: expand("{sample}/genmap_output/{sample}_genomic.genmap.sorted.bedgraph", sample=genomes)

rule genmap_index:
    input: 
        fasta = "{sample}/{sample}_genomic.fna",
    log: "{sample}/genmap_logs/index_log.txt"
    output: directory("{sample}/genmap_index")
    shell: "genmap index -F {input.fasta} -I {output} &> {log}"

rule genmap_map:
    input: 
        index = "{sample}/genmap_index",
    log: "{sample}/genmap_logs/map_log.txt"
    params: outdir = directory("{sample}/genmap_output") 
    output: bg = temp("{sample}/genmap_output/{sample}_genomic.genmap.bedgraph")     
    shell: "genmap map -K 150 -E 0 -I {input.index} -O {params.outdir} -bg -T 2 -v  > {log}"

rule sort:
    input: "{sample}/genmap_output/{sample}_genomic.genmap.bedgraph"
    output: "{sample}/genmap_output/{sample}_genomic.genmap.sorted.bedgraph"
    shell: "sort -k1,1 -k2,2n {input} > {output}" 