rule vcf_annotate:
	"""
	This rule annotates clean VCFs with snpEff then parses out the missense and synonymous variants to a BED file to be used in the gene_annot rule
	"""
	input:
		ingroup = config['ingroup'] + ".clean.vcf",
		outgroup = config['outgroup'] + ".clean.vcf"
	output:
		ingroup = config['ingroup'] + ".ann.vcf",
		outgroup = config['outgroup'] + ".ann.vcf"
	params:
		snpEffGenome = config['ingroup']
	shell:
		"snpEff ann -i vcf -o vcf -c snpEff/snpEff.config {params.snpEffGenome} {input.ingroup} > {output.ingroup}\n"
		"snpEff ann -i vcf -o vcf -c snpEff/snpEff.config {params.snpEffGenome} {input.outgroup} > {output.outgroup}"
		
rule vcf_parse:
	"""
	This rule parses the variant effects of interest from the annotated VCF and ouputs a BED file for use in the gene_annot rule
	"""
	input:
		script = "scripts/annot_parser.py",
		ingroup = config['ingroup'] + ".ann.vcf",
		outgroup = config['outgroup'] + ".ann.vcf"
	output:
		ingroup = config['ingroup'] + ".ann.bed",
		outgroup = config['outgroup'] + ".ann.bed"
	shell:
		"python3 {input.script} {input.ingroup} {output.ingroup} -key missense_variant -key synonymous_variant\n"
		"python3 {input.script} {input.outgroup} {output.outgroup} -key missense_variant -key synonymous_variant"

rule gene_annot:
	"""
	This rule intersects the annotated BED files with the CDS gene names BED to create the final BED files for use in the prep_snipre rule
	"""
	input:
		ingroup = config['ingroup'] + ".ann.bed",
		outgroup = config['outgroup'] + ".ann.bed"
	output:
		ingroup = config['ingroup'] + ".final.bed",
		outgroup = config['outgroup'] + ".final.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.ingroup}\n"
		"bedtools intersect -a {input.outgroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.outgroup}"
