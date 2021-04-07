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
              """awk -f scripts/sizes2genome.awk {output.chrom} > {output.bed}"""
