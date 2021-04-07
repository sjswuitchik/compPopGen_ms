rule calc_missingness:
	"""
	This rule calculates the proportion of missing data and outputs a list of individuals to be removed in the vcf_filter rule
	"""
	input:
		script = "scripts/missingness.R",
		ingroup = config['ingroup'] + "_missing_data.txt",
		outgroup = config['outgroup'] + "_missing_data.txt"
	output: 
		"ingroup.remove.indv",
		"outgroup.remove.indv"
	shell:
		"Rscript {input.script} {input.ingroup} {input.outgroup}"

rule callable_sites:
	"""
	This rule takes the clean coverage sites from each species and intersects to output a set of callable sites common between both species to be used in the vcf_filter rule
	"""
	input:
		ingroup = config['ingroup'] + "_coverage_sites_clean_merged.bed",
		outgroup = config['outgroup'] + "_coverage_sites_clean_merged.bed",
		map = "WHATEVER THE MAPPABILITY BED IS CALLED"
	output:
		clean = "clean.callable.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.outgroup} > callable.bed\n"
		"bedtools intersect -a {output.call} -b {PUT MAPPABILITY BED HERE} > {output.clean}"

rule cds:
	"""
	This rule pulls out the CDS regions from the GFF to be used in the cds_genes rule
	"""
	input:
		genes = "genes.gff"
	output:
		cdsGFF = "onlyCDS.gff",
		cdsBED = "onlyCDS.bed"
	shell:
		"""awk -f helper_scripts/cds.awk {input.genes} > {output.cdsGFF}\n"""
		"""awk -f helper_scripts/gff2bed.awk {output.cdsGFF} > {output.cdsBED}"""

rule cds_genes:
	"""
	This rule associates the CDS regions from the GFF with the gene names to be used in the gene_annot rule
	"""
	input:
		bed = "onlyCDS.bed"
	output:
		bed = "onlyCDS.genes.bed"
	shell:
		"""awk -F '["\t ]' -v OFS='\t' '$(NF - 1) {print $1, $2, $3, $(NF-1)}' {input.bed} > {output.bed}"""
