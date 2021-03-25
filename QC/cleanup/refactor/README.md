### Conda env build notes

`conda create -n qc -c bioconda ucsc-fatotwobit ucsc-twobitinfo bedtools ucsc-bedgraphtobigwig ucsc-bigwigmerge ucsc-bigwigaverageoverbed`  

For snakemake env build  


channels:  
  \- bioconda  
  \- defaults  
dependencies:  
  \- ucsc-fatotwobit==377  
  \- ucsc-twobitinfo==377  
  \- bedtools==2.30.0  
  \- ucsc-bedgraphtobigwig==377  
  \- ucsc-bigwigmerge==377  
  \- ucsc-bigwigaverageoverbed==377  
