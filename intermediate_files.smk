localrules: intermediate_files

rule calc_missingness:
	input:
		ingroup = "{ingroup}_missing_data.txt"
		outgroup = "{outgroup}_missing_data.txt"
	output: 
		"{ingroup}.remove.indv"
		"{outgroup}.remove.indv"
	shell:
		"Rscript --vanilla missingness.R {input.ingroup} {input.outgroup}"

rule callable_sites:
	input:
		ingroup = "{ingroup}_coverage_sites_clean_merged.bed"
		outgroup = "{outgroup}_coverage_sites_clean_merged.bed"
	output:
		call = "callable.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.outgroup} > {output.call}"

rule cds_genes:
	input:
		genes = "genes.gff"
	output:
		cds = "onlyCDS.genes.bed"
	shell:
		"column -s, -t < {input.genes} | awk '$3 == "CDS"' > onlyCDS.gff"
		"awk -f gff2bed.awk onlyCDS.gff > onlyCDS.bed"
		"cat onlyCDS.bed | python3 genenames.py > {output.cds}"
