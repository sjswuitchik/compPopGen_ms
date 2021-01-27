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

print('\nFinished!\n\n')
