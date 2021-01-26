#### First function: Checks  
- check that ingroup and outgroup files exist, number of scaffolds in both VCFs are equal, etc.; throw exception if not  
- def main  
      add arg in.vcf
      add arg out.vcf
      parse args
      check vcf func
        print 'checking vcfs'
        return True
      else
        raise Exception  

#### Second function: Missingness  
- call once for ingroup, once for outgroup
- takes missing_data_per_indv as input, returns a list of indv to remove as output (maybe use pandas?) 
- threshold for missingness can be either a fixed or scaled parameter set by the user  e.g., thresh_scaled would be scaled by the median and thresh_fixed would be a hard threshold, would have to be mutually exclusive  
- remove indv; check if list is empty, if not empty, pass to vcftools filtering step  

#### Third function: Filtering with vcftools  
- takes vcftools params but allows user to change mac, maf, and max_missing  
- add in callable sites BED filtering and remove indv (instead of remove indv in missingness?) 
- output is "final" clean VCF that is ready for analyses
- need to write in a check to make sure the VCF has at least one indv  

#### Fourth function: snpEff  
- assume user has correctly set up snpEff and database before running the pipeline  
- user will need to specific path to snpEff.jar 
- use annot_parser.py to extract synonymous and nonsynonymous sites 
- with callable sites and onlyCDS BED files, output the MK table  

### Script breakdown  
#### Script 1: VCF cleaning, suggestions for best practices  
- filter for missing, callable sites, remove sites, remove indv (functions 1 - 3)
- outputs analysis-ready VCF  

#### Script 2: snpEff + MK  
- input = clean VCF, output = MK table (function 4)  

#### Script 3: SnIPRE  
- input = clean VCF, genome FASTA, genome annotation GFF  
- output = table of gene, dn, ds, pn, ps, npop, nout, Tsil, Trepl
- run SnIPRE & MK tests (from mktest.R) and output results
