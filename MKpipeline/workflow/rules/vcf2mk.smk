rule db_build: 
	"""
	This rule builds the snpEff annotation database 
	"""
	output: 
		"data/" + config['ingroup'] + "/snpEffectPredictor.bin"
	params:
		ref = config['ingroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"snpEff -Xmx8g build -c snpEff.config -gff3 -v -noCheckCds -noCheckProtein {params.ref}"
	
rule cds:
	"""
	This rule pulls out the CDS regions from the GFF to be used in the cds_genes rule
	"""
	input:
		genes = "data/" + config['mkDir'] + "genes.gff" 
	output:
		cdsGFF = "data/" + config['mkDir'] + "onlyCDS.gff",
		cdsBED = "data/" + config['mkDir'] + "onlyCDS.bed"
	params:
		ingroup = config['ingroup']
	shell:
		"cp data/{params.ingroup}/genes.gff.gz data/mk_tests\n"
		"gunzip data/mk_tests/genes.gff.gz\n"
		"""awk -f workflow/scripts/cds.awk {input.genes} > {output.cdsGFF}\n"""
		"""awk -f workflow/scripts/gff2bed.awk {output.cdsGFF} > {output.cdsBED}"""

rule callable_sites:
	"""
	This rule takes the clean coverage sites from each species and intersects to output a set of callable sites common between both species to be used in the vcf_filter rule
	"""
	input:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".callable_sites.bed",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".callable_sites.bed"
	output:
		call = "data/" + config['mkDir'] + "callable.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"bedtools intersect -a {input.ingroup} -b {input.outgroup} > {output.call}\n"

rule cds_genes:
	"""
	This rule associates the CDS regions from the GFF with the gene names to be used in the gene_annot rule
	"""
	input:
		bed = "data/" + config['mkDir'] + "onlyCDS.bed"
	output:
		bed = "data/" + config['mkDir'] + "onlyCDS.genes.bed"
	shell:
		"""awk -v OFS="\t" "match($0, /gene=[^;]+/) {{print $1, $2, $3, substr($0, RSTART+5, RLENGTH-5)}}" {input.bed} > {output.bed}"""
		
rule callable_cds:
	"""
	This rule intersects the callable site and the gene names for CDS regions for use in the prep_snipre rule
	"""
	input:
		call = "data/" + config['mkDir'] + "callable.bed",
		cds = "data/" + config['mkDir'] + "onlyCDS.genes.bed"
	output:
		call = "data/" + config['mkDir'] + "callable.cds.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"bedtools intersect -a {input.call} -b {input.cds} -wb | cut -f1,2,3,7 | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct > {output.call}"

rule vcf_filter:
	"""
	This rule filters both VCFs for sites and individuals to produce filtered VCFs to be used in the vcf_call rule
	"""
	input:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".vcf.gz",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".vcf.gz",
		ingroupR = "data/" + config['mkDir'] + config['ingroup'] + ".remove.indv",
		outgroupR = "data/" + config['mkDir'] + config['outgroup'] + ".remove.indv"
	output:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".filter.recode.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".filter.recode.vcf"
	params:
		mac = config['mac'],
		maf = config['maf'],
		mm = config['mm'],
		ingroup = config['ingroup'],
		outgroup = config['outgroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"vcftools --gzvcf {input.ingroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac {params.mac} --remove {input.ingroupR} --max-missing {params.mm} --recode --recode-INFO-all --out data/mk_tests/{params.ingroup}/{params.ingroup}.filter\n"
		"vcftools --gzvcf {input.outgroup} --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf {params.maf} --remove {input.outgroupR} --max-missing {params.mm} --recode --recode-INFO-all --out data/mk_tests/{params.ingroup}/{params.outgroup}.filter"
		
rule vcf_call: 
	"""
	This rule filters both VCFs for callable sites common between the two species to be used in the vcf_annotate rule
	"""
	input:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".filter.recode.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".filter.recode.vcf",
		call = "data/" + config['mkDir'] + "callable.bed"
	output:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".clean.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".clean.vcf"
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
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".clean.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".clean.vcf"
	output:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".ann.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".ann.vcf"
	params:
		ref = config['ingroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"snpEff ann -Xmx8g -i vcf -o vcf -c snpEff.config {params.ref} {input.ingroup} > {output.ingroup}\n"
		"snpEff ann -Xmx8g -i vcf -o vcf -c snpEff.config {params.ref} {input.outgroup} > {output.outgroup}"
		
rule vcf_parse: 
	"""
	This rule parses the variant effects of interest from the annotated VCF and ouputs a BED file for use in the gene_annot rule
	"""
	input:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".ann.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".ann.vcf"
	output:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".ann.bed",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".ann.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"python3 ../scripts/annot_parser.py {input.ingroup} {output.ingroup} -key missense_variant -key synonymous_variant\n"
		"python3 ../scripts/annot_parser.py {input.outgroup} {output.outgroup} -key missense_variant -key synonymous_variant"

rule gene_annot:
	"""
	This rule intersects the annotated BED files with the CDS gene names BED to create the final BED files for use in the prep_snipre rule
	"""
	input:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".ann.bed",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".ann.bed"
	output:
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".final.bed",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".final.bed"
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"bedtools intersect -a {input.ingroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.ingroup}\n"
		"bedtools intersect -a {input.outgroup} -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > {output.outgroup}"

rule miss_snipre:
	"""
	This rule outputs the missingness on a per-site basis for use in the prep_snipre rule 
	"""
	input: 
		ingroup = "data/" + config['mkDir'] + config['ingroup'] + ".ann.vcf",
		outgroup = "data/" + config['mkDir'] + config['outgroup'] + ".ann.vcf"
	output:
		"data/" + config['mkDir'] + config['ingroup'] + ".lmiss",
		"data/" + config['mkDir'] + config['outgroup'] + ".lmiss" 
	params:
		ingroup = config['ingroup'],
		outgroup = config['outgroup']
	conda:
		"../envs/vcfSnpEff.yml"
	shell:
		"vcftools --vcf {input.ingroup} --missing-site --out data/mk_tests/{params.ingroup}/{params.ingroup}\n"
		"vcftools --vcf {input.outgroup} --missing-site --out data/mk_tests/{params.ingroup}/{params.outgroup}"
		
rule prep_snipre:
	"""
	This rule calculates the missingness on a per-site basis and, with the final BED files, outputs an MK table that is ready for use in the mk_snipre_stats rule
	"""
	input:
		ingroupBED = "data/" + config['mkDir'] + config['ingroup'] + ".final.bed",
		outgroupBED = "data/" + config['mkDir'] + config['outgroup'] + ".final.bed",
		ingroupM = "data/" + config['mkDir'] + config['ingroup'] + ".lmiss",
		outgroupM = "data/" + config['mkDir'] + config['outgroup'] + ".lmiss",
		call = "data/" + config['mkDir'] + "callable.cds.bed",
		cds = "data/" + config['mkDir'] + "onlyCDS.genes.bed"
	output:
		"data/" + config['mkDir'] + "snipre_data.tsv"
	conda:
		"../envs/r.yml"
	shell:
		"Rscript --slave --vanilla ../scripts/prep_snipre.R {input.ingroupBED} {input.outgroupBED} {input.ingroupM} {input.outgroupM}"

rule mk_snipre_stats:
	"""
	This rule takes the MK table formatted for SnIPRE and runs an MK test, SnIPRE, and calculates a number of statistics 
	"""
	input:
		"data/" + config['mkDir'] + "snipre_data.tsv"
	output:
		"data/" + config['mkDir'] + "mk_output.tsv",
		"data/" + config['mkDir'] + "snipre_output.tsv"
	conda:
		"../envs/r.yml"
	shell:
		"Rscript --slave --vanilla ../scripts/run_snipre.R"
