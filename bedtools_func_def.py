#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:58:54 2021

@author: sjswuitchik
"""
import os, argparse, subprocess, shlex #redo this with pybedtools in time

def create_callable_sites(ingroup_short, outgroup_short, output):
    '''Create a BED file of callables sites common to both species
    '''
    output_file = 'callable.bed'
    command = 'bedtools intersect -a '+ingroup_short+'_coverage_sites_clean_merged.bed -b '+outgroup_short+'_coverage_sites_clean_merged.bed > ' + output_file
    try:
        p = subprocess.run(shlex.split(command), capture_output = True, check = True, text = True)
    except subprocess.CalledProcessError as e:
        raise e
        
    return output_file


def filter_callable_sites(group):
    '''common function for filtering VCF for callable sites with bedtools
    '''
    output_file = group+'call.vcf'
    command = 'bedtools intersect -a '+group+'.clean.vcf -b callable.bed -header > '+group+'call.vcf'
    try:
        p = subprocess.run(shlex.split(command), capture_output = True, check = True, text = True)
    except subprocess.CalledProcessError as e:
        raise e
        
    return output_file


def filter_ingroup_call_sites(ingroup_short):
    '''Intersect the ingroup VCF with callable sites BED'''
    return filter_callable_sites(ingroup_short)

def filter_outgroup_call_sites(outgroup_short):
    '''Intersect the outgroup VCF with callable sites BED'''
    return filter_callable_sites(outgroup_short)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description= 'bedtools commands to create callable sites for MK pipeline')
    
    requiredParam = parser.add_argument_group('required parameters')
    requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
    requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

    args = parser.parse_args()
