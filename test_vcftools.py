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
optionalParam.add_argument('-mm', type = str, metavar = 'max_missing', default = '0.5', help = 'Maximum missing argument for vcftools')
optionalParam.add_argument('-maf', type = str, metavar = 'minor_allele_freq', default = '1', help = 'Minor allele frequency argument for vcftools')
optionalParam.add_argument('-mac', type = str, metavar = 'minor_allele_count', default = '0', help = 'Minor allele count argument for vcftools')

args = parser.parse_args()

print('VCF filtering with vcftools')

# run vcftools on ingroup
command = ('vcftools --gzvcf '+args.i+'.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac '+args.mac+' --max-missing '+args.mm+' --remove ingroup.remove.indv --recode --recode-INFO-all --out '+args.i+'.filter')
print('Ingroup filtering command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

# run vcftools on outgroup
command = ('vcftools --gzvcf '+args.o+'.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf '+args.maf+' --max-missing '+args.mm+' --remove outgroup.remove.indv --recode --recode-INFO-all --out '+args.o+'.filter')
print('Outgroup filtering command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

print('\nFinished!\n\n')
