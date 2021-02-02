localrules: vcf2mk

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

rule vcf_filter:
	input:
		ingroup = "{ingroup}.vcf.gz"
		outgroup = "{outgroup}.vcf.gz"
	output:
		ingroup = "{ingroup}.clean.vcf"
		outgroup = "{outgroup}.clean.vcf"
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --remove ingroup.remove.indv --max-missing 0.5 --recode --recode-INFO-all --out ingroup.filter"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min--alleles 2 --max-alleles 2 --maf 0 --remove outgroup.remove.indv --max-missing 0.5 --recore --recode-INFO-all --out outgroup.filter"
		"bedtools intersect -a ingroup.filter.recode.vcf -b callable.bed -header > {output.ingroup}"
		"bedtools intersect -a outgroup.filter.recode.vcf -b callable.bed -header > {output.outgroup}"

rule vcf_annotate:
	input:
		ingroup = "{ingroup}.clean.vcf"
		outgroup = "{outgroup}.clean.vcf"
	output:
		ingroup = "{ingroup}.ann.bed"
		outgroup = "{outgroup}.ann.bed"
	shell:
		"java -jar snpEff/snpEff.jar {ingroup} {input.ingroup} > {ingroup}.ann.vcf"
		"java -jar snpEff/snpEff.jar {ingroup} {input.outgroup} > {outgroup}.ann.vcf"
		"python3 annot_parser.py {ingroup}.ann.vcf {output.ingroup} -key missense_variant -kay synonymous_variant"
		"python3 annot_parser.py {outgroup}.ann.vcf {output.outgroup} -key missense_variant -kay synonymous_variant"

rule gene_annot:
	input:
		ingroup = "{ingroup}.ann.bed"
		outgroup = "{outgroup}.ann.bed"
	output:
		ingroup = "{ingroup}.final.bed"
		outgroup = "{outgroup}.final.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.ingroup}"
		"bedtools intersect -a {input.outgroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.outgroup}"

rule prep_snipre:
	input:
		ingroupVCF = "{ingroup}.ann.vcf"
		outgroupVCF = "{outgroup}.ann.vcf"
		ingroupBED = "{ingroup}.final.bed"
		outgroupBED = "{outgroup}.final.bed"
	output:
		snipre = "snipre_data.tsv"
	shell:
		"vcftools --vcf {input.ingroupVCF} --missing-site --out {ingroup}"
		"vcftools --vcf {input.outgroupVCF} --missing-site --out {outgroup}"
		"Rscript --slave --vanilla snipre_prep.R {input.ingroupBED} {input.outgroupBED} {ingroup}.lmiss {outgroup}.lmiss > prep_std.Rout"

rule mk_snipre_stats:
	input:
		"snipre_data.tsv"
	output:
		"mk_output.tsv"
		"snipre_output.tsv"
	shell:
		"Rscript --slave --vanilla run_snipre.R > mk_std.Rout"
