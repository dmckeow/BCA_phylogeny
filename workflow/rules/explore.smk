
configfile: "config/explore.yaml"

rule all:
	input:
		"results_annotation/broccoli/dir_step3/orthologous_groups.txt",
		"results_annotation/broccoli/dir_step4/orthologous_pairs.txt",
		"results_annotation/OrthoFinder"

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

rule orthofinder:
	input:
		config["orthofinder"]["input_dir"]
	output:
		directory("results_annotation/OrthoFinder")
	params:
		tool_params=config["orthofinder"]["orthofinder_params"]
	conda:
		"../envs/initial_ortho.yaml"
	shell:
		"""
		orthofinder -f {input} -o {output} -t {threads} {params.tool_params}
		"""
