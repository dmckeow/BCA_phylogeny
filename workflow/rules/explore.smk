import os
import glob

configfile: "config/explore.yaml"

# Check for illegal config settings

if config["workflow_control"]["broccoli"] != "TRUE" and config["workflow_control"]["orthofinder"] != "TRUE":
	print("Neither 'broccoli' nor 'orthofinder' is enabled in the workflow control. You must eat your vegetables or your orthofinder. Exiting.")
	sys.exit(1)

# FUNCTIONS

# this version aggregates:
#def get_fasta_files(wildcards):
#	ckpt_output = checkpoints.orthofinder.get().output[0]  # Retrieve the checkpoint output
#	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
#	return expand("results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.fa", sample=fasta_files)

def get_fasta_samples_orthofinder():
	ckpt_output = checkpoints.orthofinder.get().output[0]  # Get the checkpoint output directory
	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
	return fasta_files

def get_fasta_samples_broccoli():
	ckpt_output = checkpoints.broccoli.get().output[0]  # Get the checkpoint output directory
	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
	return fasta_files

# RULES BEGIN

rule all:
	input:
		"results_annotation/broccoli/dir_step3/orthologous_groups.txt",
		"results_annotation/broccoli/dir_step4/orthologous_pairs.txt",
		"results_annotation/OrthoFinder/Orthogroup_Sequences",
		"results_annotation/broccoli/orthologous_groups_fastas",
		#expand("results_annotation/OrthoFinder/test/{sample}.fa", sample=get_fasta_files) # to aggregate
		expand("results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.txt", sample=lambda wildcards: get_fasta_samples_orthofinder()),
		expand("results_annotation/broccoli/orthologous_groups_fastas/{sample}.txt", sample=lambda wildcards: get_fasta_samples_broccoli())

# INITIAL ORTHOLOGY

if config["workflow_control"]["broccoli"] == "TRUE":
	print("Running INITIAL ORTHOLOGY with Broccoli")
	checkpoint broccoli:
		input:
			config["input_dir"]
		output:
			directory("results_annotation/broccoli/orthologous_groups_fastas"),
			"results_annotation/broccoli/dir_step3/orthologous_groups.txt",
			"results_annotation/broccoli/dir_step4/orthologous_pairs.txt"
		params:
			tool_params=config["broccoli"]["broccoli_params"]
		conda:
			"../envs/initial_ortho.yaml"
		shell:
			"""
			rm -fr results_annotation/broccoli &&
			mkdir -p results_annotation/broccoli &&
			python Broccoli/broccoli.py -dir {input} -threads {threads} {params.tool_params} &&
			mv dir_step1 results_annotation/broccoli &&
			mv dir_step2 results_annotation/broccoli &&
			mv dir_step3 results_annotation/broccoli &&
			mv dir_step4 results_annotation/broccoli &&
			python workflow/scripts/parse_fastas_broccoli.py \
				-b results_annotation/broccoli/dir_step3/orthologous_groups.txt \
				-f {input} \
				-o results_annotation/broccoli/orthologous_groups_fastas
			"""

if config["workflow_control"]["orthofinder"] == "TRUE":
	print("Running INITIAL ORTHOLOGY with OrthoFinder")
	checkpoint orthofinder:
		input:
			config["input_dir"]
		output:
			directory("results_annotation/OrthoFinder/Orthogroup_Sequences")
		params:
			tool_params=config["orthofinder"]["orthofinder_params"]
		conda:
			"../envs/initial_ortho.yaml"
		shell:
			"""
			rm -fr results_annotation/OrthoFinder_tmp &&
			rm -fr results_annotation/OrthoFinder &&
			orthofinder -n OrthoFinder -f {input} -o results_annotation/OrthoFinder_tmp -t {threads} {params.tool_params} &&
			mv results_annotation/OrthoFinder_tmp/Results_OrthoFinder results_annotation/OrthoFinder &&
			rm -fr results_annotation/OrthoFinder_tmp
			"""



# example:
#rule intermediate:
#	input:
		#get_fasta_files # to aggregate
#		"results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.fa"
#	output:
#		"results_annotation/OrthoFinder/test/{sample}.fa"
#	shell:
#		"cp {input} {output}"


# SEARCH

## HMM search for orthofinder
if config["workflow_control"]["search"] == "TRUE" and config["workflow_control"]["orthofinder"] == "TRUE":
	print("Running SEARCH with hmmsearch against your hmms of interest for OrthoFinder results")
	rule search_orthofinder:
		input:
			"results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.fa"
		output:
			"results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.txt"
		params:
			gene_family_info=config["search"]["gene_family_info"],
			gene_family_name=config["search"]["gene_family_name"],
			hmm_dir=config["search"]["hmm_dir"],
			input_source="OrthoFinder"
		conda:
			"../envs/search_cluster.yaml"
		shell:
			"""
			python workflow/scripts/s01_search.py -f {input} -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir} -i {params.input_source} &&
			touch {output}
			"""
else: 
	print("hmm search step SKIPPED for OrthoFinder")

## HMM search for broccoli

if config["workflow_control"]["search"] == "TRUE" and config["workflow_control"]["broccoli"] == "TRUE":
	print("Running SEARCH with hmmsearch against your hmms of interest for broccoli results")
	rule search_broccoli:
		input:
			"results_annotation/broccoli/orthologous_groups_fastas/{sample}.fa"
		output:
			"results_annotation/broccoli/orthologous_groups_fastas/{sample}.txt"
		params:
			gene_family_info=config["search"]["gene_family_info"],
			gene_family_name=config["search"]["gene_family_name"],
			hmm_dir=config["search"]["hmm_dir"],
			input_source="broccoli"
		conda:
			"../envs/search_cluster.yaml"
		shell:
			"""
			python workflow/scripts/s01_search.py -f {input} -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir} -i {params.input_source} &&
			touch {output}
			"""
else: 
	print("hmm search step SKIPPED for broccoli")

# CLUSTER

if config["workflow_control"]["diamond_mcl"] == "TRUE":
	print("Running clustering using diamond and MCL")
	if config["workflow_control"]["search"] == "TRUE":
		if config["workflow_control"]["broccoli"] == "TRUE":
			rule cluster_diamond_mcl_broccoli:
				input:
					"results_annotation/searches/broccoli/{sample}.domains.fasta"
				output:
					"results_annotation/clusters/broccoli/{sample}.domains.abc"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"]
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input}  -q {input} -o {input}.csv --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {input}.csv > {input}.csv.tmp &&
					mv {input}.csv.tmp {input}.csv &&
					mcl {input}.csv {params.mcl_params} -I {mcl_inflation} -o {output}
					"""
		if config["workflow_control"]["orthofinder"] == "TRUE":
			rule cluster_diamond_mcl_orthofinder:
				input:
					"results_annotation/searches/OrthoFinder/{sample}.domains.fasta"
				output:
					"results_annotation/clusters/OrthoFinder/{sample}.domains.abc"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"]
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input}  -q {input} -o {input}.csv --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {input}.csv > {input}.csv.tmp &&
					mv {input}.csv.tmp {input}.csv &&
					mcl {input}.csv {params.mcl_params} -I {mcl_inflation} -o {output}
					"""

