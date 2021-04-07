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
  
