#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:49:24 2021

@author: sjswuitchik
"""
import os, argparse, subprocess, shlex, annot_parser

def annotate_group_VCF(snpEff_path, group_1, group_2, group_out):
    ''' common function for calling snpEff.jar and annotating groups
        
        snpEff_path - directory of snpEff.jar
        group_1 - 1st arg to snpEff.jar
        group_2 - 2nd arg to snpEff.jar:   group_2.call.vcf
        group_out - output file prefix:    group_out.ann.vcf
    
    '''
    output_file = group_out+'ann.vcf'
    command = 'java -jar '+snpEff_path+'/snpEff.jar ' + group_1 + ' ' + group_2 + '.call.vcf > '+ output_file
    try:
        # Run the subprocess, capture its output, throw an exception if there's an error,
        # and return the output as text not binary strings.
        p = subprocess.run(shlex.split(command), capture_output=True, check=True, text=True) 
    except subprocess.CalledProcessError as e:
        # error in subprocess!
        raise e 
 
    return output_file     

def annotate_ingroup_VCF(snpEff_path,ingroup_short):
    ''' Annotate the ingroup VCF with the ingroup genome as reference'''
    return annotate_group_VCF(snpEff_path, ingroup_short, ingroup_short, ingroup_short)   

def annotate_outgroup_VCF(snpEff_path,ingroup_short, outgroup_short):
    '''Annotate the outgroup VCF with the ingroup genome as reference'''
    return annotate_group_VCF(snpEff_path, ingroup_short, outgroup_short, outgroup_short)   

def parse_variants(outfile,infile,keys):
    ''' Description goes here
        keys should be a list of keys '''
    return annot_parser.proc_file(outfile,infile,keys)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description= 'snpEff commands & cleaning for MK pipeline')

    # Note the addition of the key for the annot_parser call
    requiredParam = parser.add_argument_group('required parameters')
    requiredParam.add_argument('-i', type = str, metavar = 'ingroup_short', required = True, help = 'Base name for ingroup species')
    requiredParam.add_argument('-o', type = str, metavar = 'outgroup_short', required = True, help = 'Base name for outgroup species')
    requiredParam.add_argument('-p', type = str, metavar = 'snpEff_path', required = True, help = 'Path for snpEff directory without terminal forward slash')
    requiredParam.add_argument('-key',help='Effects to extract. Can be repeated, must be mutually exclusive', action='append')

    args = parser.parse_args()
