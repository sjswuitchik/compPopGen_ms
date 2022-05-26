#!/bin/bash

## useage: ./gff_chr_repl_athCun.sh will overwrite the scaffold names in the GFF provided to the sed command based on the names in the translator file


#Loop through each chromosome and scaffold
cat athCun_ids_translator.tsv | while read LINE
do

GFF_CHR=$(echo $LINE | awk '{print $1}')
NCBI_ID=$(echo $LINE | awk '{print $2}')

#Sed to find and replace chromosome names
sed -i_translated "s/\b${GFF_CHR}\b/${NCBI_ID}/g"  genes.gff

done



