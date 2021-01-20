#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:54:39 2021

@author: sjswuitchik
"""
import os, argparse, subprocess

parser = argparse.ArgumentParser(description= 'Prep for and running MK test, SnIPRE, DoS calculations')

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

args = parser.parse_args()

print('MK test, SnIPRE, DoS calculations in R')

command = ('Rscript --slave --vanilla mktest.R '+args.i+'.final.bed '+args.o+'.final.bed '+args.i+'.missingness.lmiss '+args.o+'.missingness.lmiss >std.Rout 2>std.Rerr')
print('Ingroup annotating command: '+command)
p = subprocess.Popen(command, shell = True)
sts = os.waitpid(p.pid, 0)[1]
