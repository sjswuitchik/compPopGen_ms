localrules: covCalcs
  
rule prep_genome:
      """
      Describe rule
      """
      input: 
              genome = config['spp] + ".fa"
      output:
              twobit = config['spp'] + ".2bit",
              chrom = config['spp'] + ".chrom.sizes",
              bed = config['spp'] + ".genome.bed"
      params:
              spp = config['spp']
      shell:
              "faToTwoBit -long {input.genome} {output.twobit}\n"
              "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}\n"
              """awk -f helper_scripts/sizes2genome.awk {output.chrom} > {output.bed}"           
