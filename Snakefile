configfile: "config.yaml"

rule search:
	input:
		fasta=config["search"]["fasta"],
		gene_family_info=config["search"]["gene_family_info"]
	output:
		"results/{prefix}.{gene_family}.domains.fasta"
	params:
		gene_family_name=config["search"]["gene_family_name"]
	log:
		"logs/search/test.log"
	benchmark:
		"benchmarks/search/test.benchmark.txt"
	conda:
	  "envs/search.yaml"
	shell:
		"python helper/s01_search.py -f {input.fasta} -g {input.gene_family_info} {params.gene_family_name}"




