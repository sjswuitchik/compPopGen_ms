localrules: vcf2mk

rule calc_missingness:
	"""
	This rule calculates the proportion of missing data and outputs a list of individuals to be removed in the vcf_filter rule
	"""
	input:
		ingroup = config['ingroup'] + "_missing_data.txt",
		outgroup = config['outgroup'] + "_missing_data.txt"
	output: 
		config['ingroup'] + ".remove.indv",
		config['outgroup'] + ".remove.indv"
	shell:
		"Rscript --vanilla helper_scripts/missingness.R {input.ingroup} {input.outgroup}"

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

rule cds_genes:
	"""
	This rule pulls out the CDS regions from the GFF and parses the gene names to be used in the gene_annot rule
	"""
	input:
		genes = "genes.gff"
	output:
		cds = "onlyCDS.genes.bed"
	shell:
		"column -s, -t < {input.genes} | awk '$3 == "CDS"' > onlyCDS.gff\n"
		"awk -f helper_scripts/gff2bed.awk onlyCDS.gff > onlyCDS.bed\n"
		"cat onlyCDS.bed | python3 helper_scripts/genenames.py > {output.cds}"

rule vcf_filter:
	"""
	This rule filters both VCFs for sites and individuals to produce clean VCFs to be used in the vcf_annotate rule
	"""
	input:
		ingroup = config['ingroup'] + ".vcf.gz",
		outgroup = config['outgroup'] + ".vcf.gz"
	output:
		ingroup = config['ingroup'] + ".clean.vcf",
		outgroup = config['outgroup'] + ".clean.vcf"
	params:
		mac = config['mac'],
		maf = config['maf'],
		mm = config['mm']
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove ingroup.remove.indv --max-missing {params.mm} --recode --recode-INFO-all --out ingroup.filter\n"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min--alleles 2 --max-alleles 2 --maf {params.maf} --remove outgroup.remove.indv --max-missing {params.mm} --recore --recode-INFO-all --out outgroup.filter\n"
		"bedtools intersect -a ingroup.filter.recode.vcf -b callable.bed -header > {output.ingroup}\n"
		"bedtools intersect -a outgroup.filter.recode.vcf -b callable.bed -header > {output.outgroup}"

rule vcf_annotate:
	"""
	This rule annotates clean VCFs with snpEff then parses out the missense and synonymous variants to a BED file to be used in the gene_annot rule
	"""
	input:
		ingroup = config['ingroup'] + ".clean.vcf"
		outgroup = config['outgroup'] + ".clean.vcf"
	output:
		ingroup = config['ingroup'] + ".ann.bed"
		outgroup = config['outgroup'] + ".ann.bed"
	params:
		snpEffGenome = config['ingroup']
		ingroup = config['ingroup']
		outgroup = config['outgroup']
	shell:
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.ingroup} > {params.ingroup}.ann.vcf\n"
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.outgroup} > {params.outgroup}.ann.vcf\n"
		"python3 helper_scripts/annot_parser.py {params.ingroup}.ann.vcf {output.ingroup} -key missense_variant -key synonymous_variant\n"
		"python3 helper_scripts/annot_parser.py {params.outgroup}.ann.vcf {output.outgroup} -key missense_variant -key synonymous_variant"

rule gene_annot:
	"""
	This rule intersects the annotated BED files with the CDS gene names BED to create the final BED files for use in the prep_snipre rule
	"""
	input:
		ingroup = config['ingroup'] + ".ann.bed"
		outgroup = config['outgroup'] + ".ann.bed"
	output:
		ingroup = config['ingroup'] + ".final.bed"
		outgroup = config['outgroup'] + ".final.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.ingroup}\n"
		"bedtools intersect -a {input.outgroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.outgroup}"

rule prep_snipre:
	"""
	This rule calculates the missingness on a per-site basis and, with the final BED files, outputs an MK table that is ready for use in the mk_snipre_stats rule
	"""
	input:
		ingroupVCF = config['ingroup'] + ".ann.vcf"
		outgroupVCF = config['outgroup'] + ".ann.vcf"
		ingroupBED = config['ingroup'] + ".final.bed"
		outgroupBED = config['outgroup'] + ".final.bed"
	output:
		snipre = "snipre_data.tsv"
	params:
		ingroup = config['ingroup']
		outgroup = config['outgroup']
	shell:
		"vcftools --vcf {input.ingroupVCF} --missing-site --out {params.ingroup}\n"
		"vcftools --vcf {input.outgroupVCF} --missing-site --out {params.outgroup}\n"
		"Rscript --slave --vanilla helper_scripts/prep_snipre.R {input.ingroupBED} {input.outgroupBED} {params.ingroup}.lmiss {params.outgroup}.lmiss > prep_std.Rout"

rule mk_snipre_stats:
	"""
	This rule takes the MK table formatted for SnIPRE and runs an MK test, SnIPRE, and calculates a number of statistics 
	"""
	input:
		"snipre_data.tsv"
	output:
		"mk_output.tsv"
		"snipre_output.tsv"
	shell:
		"Rscript --slave --vanilla helper_scripts/run_snipre.R > mk_std.Rout"
