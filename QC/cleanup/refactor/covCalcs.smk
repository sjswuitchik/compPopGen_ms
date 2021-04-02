localrules: covCalcs
  
rule prep_genome:
      """
      Describe rule
      """
      input: 
              genome = config['spp] + ".fa"             
      output:
              twobit = bamsDir + config['spp'] + ".2bit",
              chrom = bamsDir + config['spp'] + ".chrom.sizes",
              bed = bamsDir + config['spp'] + ".genome.bed"                        
      shell:
              "faToTwoBit -long {input.genome} {output.twobit}\n"
              "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}\n"
              """awk -f helper_scripts/sizes2genome.awk {output.chrom} > {output.bed}"""

rule make_bigWigs:
      """
      Describe rule
      """
      input:
              bam = bamsDir + "{sample}" + bam_suffix,
              chrom = bamsDir + config['spp'] + ".chrom.sizes" 
      output:
              bg = bamsDir + "{sample}" + ".bg",
              bgsort = bamsDir + "{sample}" + ".sorted.bg",
              bw = bamsDir + "{sample}" + ".bw"                 
      shell:
              "bedtools genomecov -bga -ibam {input.bam} -g {input.chrom} > {output.bg}\n"
              "sort -k1,1 -k2,2n {output.bg} > {output.bgsort}\n"
              "bedGraphToBigWig {output.bgsort} {input.chrom} {output.bw}"
                              
rule compress_bedgraphs:
      """
      Describe rule
      """
      input:
              bg = bamsDir + "{sample}" + ".bg",
              bgsort = bamsDir + "{sample}" + ".sorted.bg"                   
      output:                        
              bgzip = bamsDir + "{sample}" + ".bg.gz",
              bgsortzip = bamsDir + "{sample}" + ".sorted.bg.gz"                
      shell:                       
              "gzip {input.bg}\n"                
              "gzip {input.bgsort}"                
                              
rule bigWig_summary:
      """
      Describe rule
      """
      input:                        
              bws = bamsDir + "{sample}" + ".bw",
              chrom = bamsDir + config['spp'] + ".chrom.sizes",
              bed = bamsDir + config['spp'] + ".genome.bed"                
      output:                        
              sample_list = bamsDir + "list",
              bgmerge = bamsDir + config['spp'] + ".merge.bg",
              bwmerge = bamsDir + config['spp'] + ".merge.bw",                
              summ = bamsDir + config['spp'] + ".summary.tab"                
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
              summ = bamsDir + config['spp'] + ".summary.tab",                  
              bgmerge = bamsDir + config['spp'] + ".merge.bg.gz"                
      output:
              clean = bamsDir + config['spp'] + "_coverage_sites_clean.bed",
              low = bamsDir + config['spp'] + "_coverage_sites_low.bed", 
              high = bamsDir + config['spp'] + "_coverage_sites_high.bed"                 
      params:
              spp = config['spp]                 
      shell:                        
              ## talk to Tim about best way to snakemake this awkward awk
                           
                         
                             
                              
                              
