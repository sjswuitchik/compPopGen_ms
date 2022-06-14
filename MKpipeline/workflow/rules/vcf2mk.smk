localrules: vcf2mk
	
rule db_build: 
	"""
	This rule builds the snpEff annotation database 
	"""
	conda:
		"../envs/vcfSnpEff.yml"
	log:
		"logs/snpEff/{Organism}_dbBuild.txt"
	params:
		ref = config['ingroup'],
		config = config['snpEffDir'] + "snpEff.config"
	shell:
		"snpEff -Xmx8g build -c {params.config} -gff3 -v -noCheckCds -noCheckProtein {params.ref}"
	
rule cds:
	"""
	This rule pulls out the CDS regions from the GFF to be used in the cds_genes rule
	"""
	input:
		genes = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "genes.gff" 
	output:
		cdsGFF = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "onlyCDS.gff",
		cdsBED = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "onlyCDS.bed"
	script:
		cds = "../scripts/cds.awk",
		bed = "../scripts/gff2bed.awk"
	shell:
		"""awk -f {script.cds} {input.genes} > {output.cdsGFF}\n"""
		"""awk -f {script.bed} {output.cdsGFF} > {output.cdsBED}"""	

rule callable_sites:
	"""
	This rule takes the clean coverage sites from each species and intersects to output a set of callable sites common between both species to be used in the vcf_filter rule
	"""
	input:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".callable_sites.bed",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".callable_sites.bed"
	output:
		call = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "callable.bed"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.outgroup} > {output.call}\n"

rule cds_genes:
	"""
	This rule associates the CDS regions from the GFF with the gene names to be used in the gene_annot rule
	"""
	input:
		bed = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "onlyCDS.bed"
	output:
		bed = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "onlyCDS.genes.bed"
	shell:
		"""awk -v OFS='\t' 'match($0, /gene=[^;]+/) {print $1, $2, $3, substr($0, RSTART+5, RLENGTH-5)}' {input.bed} > {output.bed}"""
		
rule callable_cds:
	"""
	This rule intersects the callable site and the gene names for CDS regions for use in the prep_snipre rule
	"""
	input:
		call = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "callable.bed",
		cds = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "onlyCDS.genes.bed"
	output:
		call = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "callable.cds.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"bedtools intersect -a {input.call} -b {input.cds} -wb | cut -f1,2,3,7 | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct > {output.call}"

rule vcf_filter:
	"""
	This rule filters both VCFs for sites and individuals to produce filtered VCFs to be used in the vcf_call rule
	"""
	input:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".vcf.gz",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".vcf.gz",
		ingroupR = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".remove.indv",
		outgroupR = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".remove.indv"
	output:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".filter.recode.vcf",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".filter.recode.vcf"
	params:
		mac = config['mac'],
		maf = config['maf'],
		mm = config['mm'],
		ingroup = config['ingroup'],
		outgroup = config['outgroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove {input.ingroupR} --max-missing {params.mm} --recode --recode-INFO-all --out {params.ingroup}.filter\n"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf {params.maf} --remove {input.outgroupR} --max-missing {params.mm} --recode --recode-INFO-all --out {params.outgroup}.filter"
		
rule vcf_call: 
	"""
	This rule filters both VCFs for callable sites common between the two species to be used in the vcf_annotate rule
	"""
	input:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".filter.recode.vcf",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".filter.recode.vcf",
		call = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + "callable.bed"
	output:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".clean.vcf",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".clean.vcf"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.call} -header > {output.ingroup}\n"
		"bedtools intersect -a {input.outgroup} -b {input.call} -header > {output.outgroup}"

rule vcf_annotate:
	"""
	This rule annotates clean VCFs with snpEff then parses out the missense and synonymous variants to a BED file to be used in the gene_annot rule
	"""
	input:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".clean.vcf",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".clean.vcf"
	output:
		ingroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['ingroup'] + ".ann.vcf",
		outgroup = config['output'] + "{Organism}/{refGenome}/" + config['mkDir'] + config['outgroup'] + ".ann.vcf"
	params:
		ref = config['ingroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"snpEff ann -Xmx8g -i vcf -o vcf -c snpEff/snpEff.config {params.ref} {input.ingroup} > {output.ingroup}\n"
		"snpEff ann -Xmx8g -i vcf -o vcf -c snpEff/snpEff.config {params.ref} {input.outgroup} > {output.outgroup}"
		
rule vcf_parse: #### start here after call
	"""
	This rule parses the variant effects of interest from the annotated VCF and ouputs a BED file for use in the gene_annot rule
	"""
	input:
		ingroup = config['ingroup'] + ".ann.vcf",
		outgroup = config['outgroup'] + ".ann.vcf"
	output:
		ingroup = config['ingroup'] + ".ann.bed",
		outgroup = config['outgroup'] + ".ann.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	script:
		"../scripts/annot_parser.py"
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
