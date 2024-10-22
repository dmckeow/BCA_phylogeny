import os
import glob

configfile: "config/explore.yaml"

# this version aggregates:
#def get_fasta_files(wildcards):
#	ckpt_output = checkpoints.orthofinder.get().output[0]  # Retrieve the checkpoint output
#	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
#	return expand("results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.fa", sample=fasta_files)

def get_fasta_samples():
	ckpt_output = checkpoints.orthofinder.get().output[0]  # Get the checkpoint output directory
	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
	return fasta_files


rule all:
	input:
		"results_annotation/broccoli/dir_step3/orthologous_groups.txt",
		"results_annotation/broccoli/dir_step4/orthologous_pairs.txt",
		"results_annotation/OrthoFinder/Orthogroup_Sequences",
		#expand("results_annotation/OrthoFinder/test/{sample}.fa", sample=get_fasta_files) # to aggregate
		expand("results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.txt", sample=lambda wildcards: get_fasta_samples())

# INITIAL ORTHOLOGY ASSIGNMENT

rule broccoli:
	input:
		config["broccoli"]["input_dir"]
	output:
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
		python resources/Broccoli/broccoli.py -dir {input} -threads {threads} {params.tool_params} &&
		mv dir_step1 results_annotation/broccoli &&
		mv dir_step2 results_annotation/broccoli &&
		mv dir_step3 results_annotation/broccoli &&
		mv dir_step4 results_annotation/broccoli
		"""

checkpoint orthofinder:
	input:
		config["orthofinder"]["input_dir"]
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




rule search:
	input:
		"results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.fa"
	output:
		temp("results_annotation/OrthoFinder/Orthogroup_Sequences/{sample}.txt")
	params:
		gene_family_info=config["search"]["gene_family_info"],
		gene_family_name=config["search"]["gene_family_name"],
		hmm_dir=config["search"]["hmm_dir"]
	conda:
		"../envs/search_cluster.yaml"
	shell:
		"""
		python workflow/scripts/s01_search.py -f {input} -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir} &&
		touch {output}
		"""
