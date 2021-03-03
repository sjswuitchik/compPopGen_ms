localrules: vcfqc

rule missingness:
  """
  Calculate relative missingness per individual to output both a table and a dotplot
  """
  input:
          script = "helper_scripts/missingness.R"
  output:
          table = "relative_missing_per_ind.txt",
          plot = "relative_missing.pdf"
  shell:
          "Rscript {input.script}"

rule plink_pruning:
  """
  (Aggressively) prune for LD to process for PCA
  """
  input:
          vcf = "Combined_hardFiltered.vcf"
  output:
          prune = config['ingroup'] + ".prune.in",
          ld = config['ingroup'] + ".ld_pruned"
  params:
          ingroup = config['ingroup']
  shell:
          "plink --vcf {input.vcf} --make-bed --out {params.ingroup} --allow-extra-chr\n"
          "plink --bfile {params.ingroup} --indep-pairwise 500 50 0.1 --out {params.ingroup} --allow-extra-chr\n"
          "plink --bfile {params.ingroup} --make-bed --extract {output.prune} --out {output.ld} --allow-extra-chr --geno 0.95"

rule pca_raw:
  """
  Plot PCA of raw data
  """
  input:
          script = "helper_scripts/PCA.R",
          val = config['ingroup'] + ".eigenval",
          vec = config['ingroup'] + ".eigencev"
  output:
          plot = "PCA.pdf"
  shell:
          "Rscript {input.script} {input.val} {input.vec}"

rule relatedness:
  """
  Output relatedness statistic from Yang et al. (2010) (doi:10.1038/ng.608)
  """
  input:
          vcf = "Combined_hardFiltered.vcf"
  output:
          rel = config['ingroup'] + ".relatedness"
  params:
          ingroup = config['ingroup']
  shell:
          "vcftools --vcf {input.vcf} --out {params.ingroup} --relatedness

rule matrix:
  """
  Creating a 012 matrix for PCA input
  """
  input:
        vcf = "Combined_hardFiltered.vcf"
  output:
        pos = config['ingroup'] + ".012.pos",
        indv = config['ingroup'] + ".012.indv",
        matrix = config['ingroup'] + ".012"
  params:
        ingroup = config['ingroup']
  shell:
        "vcftools --vcf {input.vcf} --out {params.ingroup} --012"

rule pca_clean:
  """
  Rule description
  """
  input:
        script = "helper_scripts/SOMETHING.R"
        pos = config['ingroup'] + ".012.pos"
        indv = config['ingroup'] + ".012.indv"
        matrix = config['ingroup'] + ".012"
  output:
        SOMETHING
  shell:
        "Rscript {input.script} ... probably commandArgs 







