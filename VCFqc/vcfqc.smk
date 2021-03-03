localrules: vcfqc

rule missingness:
  """
  Rule description
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
  Rule description
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
  Rule description
  """
  input:
          script = "PCA.R",
          ld = config['ingroup'] + ".ld_pruned"
  output:
          plot = SOMETHING
  shell:
          "RScript {input.script}
