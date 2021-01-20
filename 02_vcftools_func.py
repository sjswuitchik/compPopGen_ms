#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:03:25 2021

@author: sjswuitchik
"""
import os, argparse, subprocess

parser = argparse.ArgumentParser(description= 'vcftools commands for filtering in MK pipeline')

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

optionalParam = parser.add_argument_group('optional parameters')
optionalParam.add_argument('-mm', type = str, metavar = 'max_missing', default = '0.75', help = 'Maximum missing argument for vcftools')
optionalParam.add_argument('-maf', type = str, metavar = 'minor_allele_freq', help = 'Minor allele frequency argument for vcftools')
optionalParam.add_argument('-mac', type = str, metavar = 'minor_allele_count', help = 'Minor allele count argument for vcftools')

args = parser.parse_args()

print('VCF filtering with vcftools')

# run vcftools on ingroup
command = ('vcftools --gzvcf '+args.i+'.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac '+args.mac+' --max-missing '+args.mm+' --recode --recode-INFO-all --out '+args.i+'.filter')
print('Ingroup filtering command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('export ININDV=`cat ingroup.remove.indv | wc -l` \n'
           'if [ $ININDV -gt 0 ]\n'
           'then\n'
           'vcftools --gzvcf '+args.i+'.filter.vcf.gz --remove-indv ingroup.remove.indv --out '+args.i+'.clean.vcf.gz\n'
           'else\n'
           'echo "No individuals to remove from ingroup"\n'
           'fi')
print('Ingroup missingness command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# run vcftools on outgroup
command = ('vcftools --gzvcf '+args.o+'.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf '+args.maf+' --max-missing '+args.mm+' --recode --recode-INFO-all --out '+args.o+'.clean')
print('Outgroup filtering command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('export OUTINDV=`cat outgroup.remove.indv | wc -l` \n'
           'if [ $OUTINDV -gt 1 ]\n'
           'then\n'
           'vcftools --gzvcf '+args.o+'.filter.vcf.gz --remove-indv outgroup.remove.indv --out '+args.o+'.clean.vcf.gz\n'
           'else\n'
           'echo "No individuals to remove from outgroup"\n'
           'fi')
print('Outgroup missingness command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# create callable sites for ingroup and outgroup
parser = argparse.ArgumentParser(description= 'bedtools command for creating callable sites in MK pipeline')

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

args = parser.parse_args()

print('Callable sites filtering with bedtools')

command = ('bedtools intersect -a '+args.i+'_coverage_sites_clean_merged.bed -b '+args.o+'_coverage_sites_clean_merged.bed > callable.bed')
print('Callable sites command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('bedtools intersect -a '+args.i+'.clean.vcf.gz -b callable.bed -header > '+args.i+'call.vcf')
print('Ingroup callable sites intersect command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('bedtools intersect -a '+args.o+'.clean.vcf.gz -b callable.bed -header > '+args.o+'call.vcf')
print('Outgroup callable sites intersect command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]


print('\nFinished!\n\n')

