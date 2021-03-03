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

rule rule2:
  """
  Rule description
  """
  input:
          
