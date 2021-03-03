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
          "RScript {input.script}"

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

rule pca:
  """
  Plot PCA
  """
  input:
          script = "helper_scripts/PCA.R",
          val = config['ingroup'] + ".eigenval",
          vec = config['ingroup'] + ".eigencev"
  output:
          plot = "PCA.pdf"
  shell:
          "RScript {input.script} {input.val} {input.vec}"

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





