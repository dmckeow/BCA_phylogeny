configfile: "config/search.yaml"

def get_fasta_inputs(wildcards):
	return config["samples"][wildcards.sample]

rule all:
	input:
		"results_annotation/searches/{wildcards.sample}_search_done.txt"

rule search:
	input:
		get_fasta_inputs
	output:
		"results_annotation/searches/{wildcards.sample}_search_done.txt" # Dummy output
	params:
		gene_family_info=config["search"]["gene_family_info"],
		gene_family_name=config["search"]["gene_family_name"],
		hmm_dir=config["search"]["hmm_dir"]
	conda:
		"../envs/search_cluster.yaml"
	shell:
		"""
		python workflow/scripts/s01_search.py -f {input}.fa -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir}
		"""
