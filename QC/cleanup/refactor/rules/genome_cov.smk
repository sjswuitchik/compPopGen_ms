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
      conda:
              "../envs/coverage.yml"
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
