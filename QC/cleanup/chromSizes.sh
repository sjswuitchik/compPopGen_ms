#!/bin/bash

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir 
# ./chromSizes.sh spp_name

coverage/./faToTwoBit ../$1/genome/$1.fa coverage/$1.2bit
coverage/./twoBitInfo coverage/$1.2bit stdout | sort -k2rn > coverage/$1.chrom.sizes
