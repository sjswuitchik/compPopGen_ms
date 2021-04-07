rule miss_snipre:
	"""
	This rule outputs the missingness on a per-site basis for use in the prep_snipre rule 
	"""
	input: 
		ingroup = config['ingroup'] + ".ann.vcf",
		outgroup = config['outgroup'] + ".ann.vcf"
	output:
		config['ingroup'] + ".lmiss",
		config['outgroup'] + ".lmiss" 
	params:
		ingroup = config['ingroup'],
		outgroup = config['outgroup']	
	shell:
		"vcftools --vcf {input.ingroup} --missing-site --out {params.ingroup}\n"
		"vcftools --vcf {input.outgroup} --missing-site --out {params.outgroup}"
		
rule prep_snipre:
	"""
	This rule calculates the missingness on a per-site basis and, with the final BED files, outputs an MK table that is ready for use in the mk_snipre_stats rule
	"""
	input:
		script = "helper_scripts/prep_snipre.R",
		ingroupBED = config['ingroup'] + ".final.bed",
		outgroupBED = config['outgroup'] + ".final.bed",
		ingroupM = config['ingroup'] + ".lmiss",
		outgroupM = config['outgroup'] + ".lmiss",
		call = "callable.cds.bed",
		cds = "onlyCDS.genes.bed"
	output:
		"snipre_data.tsv"
	shell:
		"Rscript --slave --vanilla {input.script} {input.ingroupBED} {input.outgroupBED} {input.ingroupM} {input.outgroupM}"

rule mk_snipre_stats:
	"""
	This rule takes the MK table formatted for SnIPRE and runs an MK test, SnIPRE, and calculates a number of statistics 
	"""
	input:
		script = "helper_scripts/run_snipre.R",
		data = "snipre_data.tsv"
	output:
		"mk_output.tsv",
		"snipre_output.tsv"
	shell:
		"Rscript --slave --vanilla {input.script}"
