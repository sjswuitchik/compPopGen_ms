#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:03:25 2021
@author: sjswuitchik
"""
import os, argparse, subprocess, shlex

def filter_VCF(group, mac, maf, mm):
    '''common function for filtering VCFs with vcftools
    group - base name for species
    mac - minor allele count
    maf - minor allele frequency
    mm - max missing
    
    '''
    output_file = group+'clean.vcf'
    command = 'vcftools --gzvcf '+group+'.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --remove ingroup.remove.indv --max-alleles 2 --mac '+mac+' --max-missing '+mm+' --out '+group+'.clean'
    try:
        p = subprocess.run(shlex.split(command), capture_output = True, check = True, text = True)
    except subprocess.CalledProcessError as e:
        raise e
        
    return output_file


def filter_ingroup_VCF(ingroup_short, mac, mm):
    '''Filter the ingroup VCF'''
    return filter_VCF(ingroup_short, mac, mm, ingroup_short)

def filter_outgroup_VCF(outgroup_short, maf, mm):
    '''Filter the outgroup VCF'''
    return filter_VCF(outgroup_short, maf, mm, outgroup_short)



if __name__=='__main__':
    parser = argparse.ArgumentParser(description= 'vcftools & bedtools commands for filtering in MK pipeline')
    
    requiredParam = parser.add_argument_group('required parameters')
    requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
    requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')

    optionalParam = parser.add_argument_group('optional parameters')
    optionalParam.add_argument('-mm', type = str, metavar = 'max_missing', default = '0.5', help = 'Maximum missing argument for vcftools')
    optionalParam.add_argument('-maf', type = str, metavar = 'minor_allele_freq', default = '1', help = 'Minor allele frequency argument for vcftools')
    optionalParam.add_argument('-mac', type = str, metavar = 'minor_allele_count', default = '0', help = 'Minor allele count argument for vcftools')

    args = parser.parse_args()
