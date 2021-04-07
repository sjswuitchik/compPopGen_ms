rule vcf_filter:
	"""
	This rule filters both VCFs for sites and individuals to produce filtered VCFs to be used in the vcf_call rule
	"""
	input:
		ingroup = config['ingroup'] + ".vcf.gz",
		outgroup = config['outgroup'] + ".vcf.gz",
		ingroupR = "ingroup.remove.indv",
		outgroupR = "outgroup.remove.indv"
	output:
		ingroup = config['ingroup'] + ".filter.recode.vcf",
		outgroup = config['outgroup'] + ".filter.recode.vcf"
	params:
		mac = config['mac'],
		maf = config['maf'],
		mm = config['mm'],
		ingroup = config['ingroup'],
		outgroup = config['outgroup']
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove {input.ingroupR} --max-missing {params.mm} --recode --recode-INFO-all --out {params.ingroup}.filter\n"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf {params.maf} --remove {input.outgroupR} --max-missing {params.mm} --recode --recode-INFO-all --out {params.outgroup}.filter"
		
rule vcf_call: 
	"""
	This rule filters both VCFs for callable sites common between the two species to be used in the vcf_annotate rule
	"""
	input:
		ingroup = config['ingroup'] + ".filter.recode.vcf",
		outgroup = config['outgroup'] + ".filter.recode.vcf",
		call = "callable.bed"
	output:
		ingroup = config['ingroup'] + ".clean.vcf",
		outgroup = config['outgroup'] + ".clean.vcf"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.call} -header > {output.ingroup}\n"
		"bedtools intersect -a {input.outgroup} -b {input.call} -header > {output.outgroup}"
