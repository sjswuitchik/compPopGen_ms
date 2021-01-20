#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:22:09 2021

@author: sjswuitchik
"""
import os, argparse, subprocess

parser = argparse.ArgumentParser(description= 'calculating missingness for MK pipeline')

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

args = parser.parse_args()

print('Identify individuals to remove in R based on missingness')

# identify ingroup and outgroup individuals to remove (if any)
command = ('Rscript --vanilla missingness.R '+args.i+'_missing_data.txt '+args.o+'_missing_data.txt')
print('Command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]

print('\nFinished!\n\n')

