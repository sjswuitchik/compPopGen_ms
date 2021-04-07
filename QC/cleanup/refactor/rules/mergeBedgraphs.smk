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
