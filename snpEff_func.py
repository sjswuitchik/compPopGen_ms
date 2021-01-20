#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:49:24 2021

@author: sjswuitchik
"""
import os, argparse, subprocess

parser = argparse.ArgumentParser(description= 'snpEff commands & cleaning for MK pipeline')

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')
requiredParam.add_argument('-p', type = str, metavar = 'snpEff_path', required = True, help = 'Path for snpEff directory without terminal forward slash')

args = parser.parse_args()

print('VCF annotation with snpEff')

# annotate ingroup VCF
command = ('java -jar '+args.p+'/snpEff.jar '+args.i+' '+args.i+'.call.vcf > '+args.i+'ann.vcf')
print('Ingroup annotating command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# annotate outgroup VCF
command = ('java -jar '+args.p+'/snpEff.jar '+args.i+' '+args.o+'.call.vcf > '+args.o+'ann.vcf')
print('Outgroup annotating command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# parse variants from annotated VCF to output an annotated BED file 
## I think this can be done with subprocess.call(), need to check if there is a way to set parameters 
command = ('python3 annot_parser.py '+args.i+'.ann.vcf '+args.i+'.ann.bed -key missense_variant -key synonymous_variant')
print('Ingroup parsing command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('python3 annot_parser.py '+args.o+'.ann.vcf '+args.o+'.ann.bed -key missense_variant -key synonymous_variant')
print('Outgroup parsing command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# create onlyCDS.bed from GFF 
command = ('column -s, -t < genes.gff | awk \'$3 == "CDS"\' > onlyCDS.gff\n'
           'awk -f gff2bed.awk onlyCDS.gff > onlyCDS.bed\n'
           'cat onlyCDS.bed | python genenames.py > onlyCDS.genes.bed')
p = subprocess.Popen(command, shell = True)
sts = os.waitpid/(p.pid, 0)[1]

# associate genes with annotations in ingroup
command = ('bedtools intersect -a '+args.i+'.ann.bed -b onlyCDS.genes.bed -ab | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > '+args.i+'.final.bed')
print('Ingroup intersect command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# associate genes with annotations in outgroup
command = ('bedtools intersect -a '+args.o+'.ann.bed -b onlyCDS.genes.bed -ab | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > '+args.o+'.final.bed')
print('Outgroup intersect command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# associate callable sites with genes 
command = ('bedtools intersect -a callable.bed -b onlyCDS.genes.bed -wb | cut -f1,2,3,7 | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct > callable.cds.bed')
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

print('\nFinished!\n\n')

