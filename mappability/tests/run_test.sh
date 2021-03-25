#!/bin/bash

cp ../mappability.smk ./
rm -r testgenome/genmap_*
snakemake -s mappability.smk --cores 2 && \
cmp testgenome/genmap_output/testgenome_genomic.genmap.sorted.bedgraph testgenome/test_output/test_genome.genmap.sorted.bedgraph && echo "Test completed successfully"
