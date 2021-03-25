#!/bin/bash

cp ../mappability.smk ./
rm -r testgenome/*genmap_*
rm -r testgenome2/*genmap_*
snakemake -s mappability.smk --cores 2 && \
cmp testgenome/testgenome_genomic_genmap_output/testgenome_genomic.genmap.sorted.bedgraph testgenome/test_output/test_genome.genmap.sorted.bedgraph && cmp testgenome2/testgenome_genomic_genmap_output/testgenome_genomic.genmap.sorted.bedgraph testgenome2/test_output/test_genome.genmap.sorted.bedgraph && echo "Test completed successfully"
