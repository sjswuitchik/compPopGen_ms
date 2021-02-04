localrules: vcf2mk

rule calc_missingness:
	input:
		ingroup = config['ingroup'] + "_missing_data.txt",
		outgroup = config['outgroup'] + "_missing_data.txt"
	output: 
		config['ingroup'] + ".remove.indv",
		config['outgroup'] + ".remove.indv"
	shell:
		"Rscript --vanilla missingness.R {input.ingroup} {input.outgroup}"

rule callable_sites:
	input:
		ingroup = config['ingroup'] + "_coverage_sites_clean_merged.bed",
		outgroup = config['outgroup'] + "_coverage_sites_clean_merged.bed"
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

rule vcf_filter:
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
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove ingroup.remove.indv --max-missing {params.mm} --recode --recode-INFO-all --out ingroup.filter"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min--alleles 2 --max-alleles 2 --maf {params.maf} --remove outgroup.remove.indv --max-missing {params.mm} --recore --recode-INFO-all --out outgroup.filter"
		"bedtools intersect -a ingroup.filter.recode.vcf -b callable.bed -header > {output.ingroup}"
		"bedtools intersect -a outgroup.filter.recode.vcf -b callable.bed -header > {output.outgroup}"

rule vcf_annotate:
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
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.ingroup} > {params.ingroup}.ann.vcf"
		"java -jar snpEff/snpEff.jar {params.snpEffGenome} {input.outgroup} > {params.outgroup}.ann.vcf"
		"python3 annot_parser.py {params.ingroup}.ann.vcf {output.ingroup} -key missense_variant -key synonymous_variant"
		"python3 annot_parser.py {params.outgroup}.ann.vcf {output.outgroup} -key missense_variant -key synonymous_variant"

rule gene_annot:
	input:
		ingroup = config['ingroup'] + ".ann.bed"
		outgroup = config['outgroup'] + ".ann.bed"
	output:
		ingroup = config['ingroup'] + ".final.bed"
		outgroup = config['outgroup'] + ".final.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.ingroup}"
		"bedtools intersect -a {input.outgroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.outgroup}"

rule prep_snipre:
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
		"vcftools --vcf {input.ingroupVCF} --missing-site --out {params.ingroup}"
		"vcftools --vcf {input.outgroupVCF} --missing-site --out {params.outgroup}"
		"Rscript --slave --vanilla prep_snipre.R {input.ingroupBED} {input.outgroupBED} {params.ingroup}.lmiss {params.outgroup}.lmiss > prep_std.Rout"

rule mk_snipre_stats:
	input:
		"snipre_data.tsv"
	output:
		"mk_output.tsv"
		"snipre_output.tsv"
	shell:
		"Rscript --slave --vanilla run_snipre.R > mk_std.Rout"
