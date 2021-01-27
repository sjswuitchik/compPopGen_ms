#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:46:13 2021

@author: sjswuitchik
"""
import os, argparse, subprocess

 parser = argparse.ArgumentParser(description= 'bedtools commands to create callable sites for MK pipeline')
    
requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

args = parser.parse_args()

command = ('bedtools intersect -a '+args.i+'_coverage_sites_clean_merged.bed -b '+args.o+'_coverage_sites_clean_merged.bed > callable.bed')
print('Callable sites command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('bedtools intersect -a '+args.i+'.filter.recode.vcf -b callable.bed -header > '+args.i+'call.vcf')
print('Ingroup callable sites intersect command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

command = ('bedtools intersect -a '+args.o+'.filter.recode.vcf -b callable.bed -header > '+args.o+'call.vcf')
print('Outgroup callable sites intersect command :'+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]


print('\nFinished!\n\n')
