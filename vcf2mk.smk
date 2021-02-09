localrules: vcf2mk

rule calc_missingness:
	"""
	This rule calculates the proportion of missing data and outputs a list of individuals to be removed in the vcf_filter rule
	"""
	input:
		script = "helper_scripts/missingness.R",
		ingroup = config['ingroup'] + "_missing_data.txt",
		outgroup = config['outgroup'] + "_missing_data.txt"
	output: 
		config['ingroup'] + ".remove.indv",
		config['outgroup'] + ".remove.indv"
	shell:
		"Rscript {input.script} {input.ingroup} {input.outgroup}"

rule callable_sites:
	"""
	This rule takes the clean coverage sites from each species and intersects to output a set of callable sites common between both species to be used in the vcf_filter rule
	"""
	input:
		ingroup = config['ingroup'] + "_coverage_sites_clean_merged.bed",
		outgroup = config['outgroup'] + "_coverage_sites_clean_merged.bed"
	output:
		call = "callable.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.outgroup} > {output.call}"

rule cds:
	"""
	This rule pulls out the CDS regions from the GFF to be used in the cds_genes rule
	"""
	input:
		genes = "genes.gff"
	output:
		cds = "onlyCDS.genes.bed"
	shell:
		"""awk "\$3 == "CDS"" {input.genes} | awk -f helper_scripts/gff2bed.awk > {output.cds}"""

rule cds_genes:
	"""
	This rule associates the CDS regions from the GFF with the gene names to be used in the gene_annot rule
	"""
	input:
		"onlyCDS.bed"
	output:
		"onlyCDS.genes.bed"
	script:
		"helper_scripts/genenames.py"

rule vcf_filter:
	"""
	This rule filters both VCFs for sites and individuals to produce clean VCFs to be used in the vcf_annotate rule
	"""
	input:
		ingroup = config['ingroup'] + ".vcf.gz",
		outgroup = config['outgroup'] + ".vcf.gz",
		ingroupR = config['ingroup'] + ".remove.indv",
		outgroupR = config['outgroup'] + ".remove.indv"
	output:
		ingroup = config['ingroup'] + ".clean.vcf",
		outgroup = config['outgroup'] + ".clean.vcf"
	params:
		mac = config['mac'],
		maf = config['maf'],
		mm = config['mm']
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove {input.ingroupR} --max-missing {params.mm} --recode --recode-INFO-all --out ingroup.filter\n"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min--alleles 2 --max-alleles 2 --maf {params.maf} --remove {input.outgroupR} --max-missing {params.mm} --recore --recode-INFO-all --out outgroup.filter\n"
		"bedtools intersect -a ingroup.filter.recode.vcf -b callable.bed -header > {output.ingroup}\n"
		"bedtools intersect -a outgroup.filter.recode.vcf -b callable.bed -header > {output.outgroup}"

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
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.ingroup} > {output.ingroup}\n"
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.outgroup} > {output.outgroup}"
		
rule vcf_parse:
	"""
	This rule ...
	"""
	input:
		ingroup = config['ingroup'] + ".ann.vcf",
		outgroup = config['outgroup'] + ".ann.vcf"
	output:
		ingroup = config['ingroup'] + ".ann.bed",
		outgroup = config['outgroup'] + ".ann.bed"
	script:
		"helper_scripts/annot_parser.py {input.ingroup} {output.ingroup} -key missense_variant -key synonymous_variant\n"
		"helper_scripts/annot_parser.py {input.outgroup} {output.outgroup} -key missense_variant -key synonymous_variant"

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

rule miss_snipre:
	"""
	This rule ... 
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
		ingroupBED = config['ingroup'] + ".final.bed",
		outgroupBED = config['outgroup'] + ".final.bed",
		ingroupM = config['ingroup'] + ".lmiss",
		outgroupM = config['outgroup'] + ".lmiss"
	output:
		"snipre_data.tsv"
	shell:
		"Rscript --slave --vanilla helper_scripts/prep_snipre.R {input.ingroupBED} {input.outgroupBED} {input.ingroupM} {input.outgroupM} > prep_std.Rout"

rule mk_snipre_stats:
	"""
	This rule takes the MK table formatted for SnIPRE and runs an MK test, SnIPRE, and calculates a number of statistics 
	"""
	input:
		"snipre_data.tsv"
	output:
		"mk_output.tsv",
		"snipre_output.tsv"
	shell:
		"Rscript --slave --vanilla helper_scripts/run_snipre.R > mk_std.Rout"
