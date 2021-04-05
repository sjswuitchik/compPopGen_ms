localrules: covCalcs
  
rule prep_genome:
      """
      Describe rule
      """
      input: 
              genome = config['spp] + ".fa"             
      output:
              twobit = bamDir + config['spp'] + ".2bit",
              chrom = bamDir + config['spp'] + ".chrom.sizes",
              bed = bamDir + config['spp'] + ".genome.bed" 
      conda:
              "../envs/coverage.yml"
      shell:
              "faToTwoBit -long {input.genome} {output.twobit}\n"
              "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}\n"
              """awk -f helper_scripts/sizes2genome.awk {output.chrom} > {output.bed}"""

rule make_bigWigs:
      """
      Describe rule
      """
      input:
              bam = bamDir + "{sample}" + bam_suffix,
              chrom = bamDir + config['spp'] + ".chrom.sizes" 
      output:
              bg = bamDir + "{sample}" + ".bg",
              bgsort = bamDir + "{sample}" + ".sorted.bg",
              bw = bamDir + "{sample}" + ".bw"                 
      shell:
              "bedtools genomecov -bga -ibam {input.bam} -g {input.chrom} > {output.bg}\n"
              "sort -k1,1 -k2,2n {output.bg} > {output.bgsort}\n"
              "bedGraphToBigWig {output.bgsort} {input.chrom} {output.bw}"
                              
rule compress_bedgraphs:
      """
      Describe rule
      """
      input:
              bg = bamDir + "{sample}" + ".bg",
              bgsort = bamDir + "{sample}" + ".sorted.bg"                   
      output:                        
              bgzip = bamDir + "{sample}" + ".bg.gz",
              bgsortzip = bamDir + "{sample}" + ".sorted.bg.gz"  
      conda:
              "../envs/coverage.yml"
      shell:                       
              "gzip {input.bg}\n"                
              "gzip {input.bgsort}"                
                              
rule bigWig_summary:
      """
      Describe rule
      """
      input:                        
              bws = bamDir + "{sample}" + ".bw",
              chrom = bamDir + config['spp'] + ".chrom.sizes",
              bed = bamDir + config['spp'] + ".genome.bed"                
      output:                        
              sample_list = bamDir + "list",
              bgmerge = bamDir + config['spp'] + ".merge.bg",
              bwmerge = bamDir + config['spp'] + ".merge.bw",                
              summ = bamDir + config['spp'] + ".summary.tab" 
      conda:
              "../envs/coverage.yml"
      shell:                        
              "ls {input.bws} > {output.sample_list}\n"
              "bigWigMerge -inList {output.sample_list} {output.bgmerge}\n"
              "bedGraphToBigWig {output.bgmerge} {input.chrom} {output.bwmerge}\n"
              "bigWigAverageOverBed {output.bwmerge} {input.bed} {output.summ}\n"                
              "gzip {output.bgmerge}"                
                              
rule coverage_beds: ## will need to add bedtools merge to this eventually?  
      """
      Describe rule
      """
      input:                        
              summ = bamDir + config['spp'] + ".summary.tab",                  
              bgmerge = bamDir + config['spp'] + ".merge.bg.gz"                
      output:
              clean = bamDir + config['spp'] + "_coverage_sites_clean.bed",
              low = bamDir + config['spp'] + "_coverage_sites_low.bed", 
              high = bamDir + config['spp'] + "_coverage_sites_high.bed"
      conda:
              "../envs/coverage.yml"
      params:
              spp = config['spp]                 
      shell:                        
              ## talk to Tim about best way to snakemake this awkward awk
                           
                         
                             
                              
                              
