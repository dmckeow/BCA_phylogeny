import os
import glob

configfile: "config/explore.yaml"

# Check for illegal config settings

if config["ww_cl"]["broccoli"] != "TRUE" and config["ww_cl"]["orthofinder"] != "TRUE":
	print("Neither 'broccoli' nor 'orthofinder' is enabled in the workflow control. You must eat your vegetables or your orthofinder. Exiting.")
	sys.exit(1)

# FUNCTIONS

# this version aggregates:
#def get_fasta_files(wildcards):
#	ckpt_output = checkpoints.orthofinder.get().output[0]  # Retrieve the checkpoint output
#	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
#	return expand("{config['output_dir']}/OrthoFinder/Orthogroup_Sequences/{sample}.fa", sample=fasta_files)

def get_fasta_samples_orthofinder():
	ckpt_output = checkpoints.orthofinder.get().output[0]  # Get the checkpoint output directory
	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
	return fasta_files

def get_fasta_samples_broccoli():
	ckpt_output = checkpoints.broccoli.get().output[0]  # Get the checkpoint output directory
	fasta_files = glob_wildcards(os.path.join(ckpt_output, "{sample}.fa")).sample
	return fasta_files

def get_orthofinder_domains_files():
	# Retrieve the actual output from the checkpoint
	ckpt_output = checkpoints.search_orthofinder_ckpt.get().output[0]
	# Extract all files with the specific naming pattern
	domain_files = glob_wildcards(os.path.join(ckpt_output, "{unique_sample_name}.domains.fasta")).unique_sample_name
	return domain_files

def get_broccoli_domains_files():
	# Retrieve the actual output from the checkpoint
	ckpt_output = checkpoints.search_broccoli_ckpt.get().output[0]
	# Extract all files with the specific naming pattern
	domain_files = glob_wildcards(os.path.join(ckpt_output, "{unique_sample_name}.domains.fasta")).unique_sample_name
	return domain_files

# RULES BEGIN

rule all:
	input:
		f"{config["output_dir"]}/broccoli/dir_step3/orthologous_groups.txt",
		f"{config["output_dir"]}/broccoli/dir_step4/orthologous_pairs.txt",
		f"{config["output_dir"]}/OrthoFinder/Orthogroup_Sequences",
		f"{config["output_dir"]}/broccoli/Orthogroup_Sequences",
		#expand("{config["output_dir"]}/OrthoFinder/test/{sample}.fa", sample=get_fasta_files) # to aggregate
		expand(f"{config["output_dir"]}/OrthoFinder/Orthogroup_Sequences/{{sample}}.txt", sample=lambda wildcards: get_fasta_samples_orthofinder()),
		expand(f"{config["output_dir"]}/broccoli/Orthogroup_Sequences/{{sample}}.txt", sample=lambda wildcards: get_fasta_samples_broccoli()),

		# CLUSTER OUTPUTS:
			# diamond mcl
		expand(f"{config["output_dir"]}/searches/OrthoFinder/{{unique_sample_name}}.domains.fasta", unique_sample_name=lambda wildcards: get_orthofinder_domains_files()),
		expand(f"{config["output_dir"]}/clusters/searches/OrthoFinder/dmnd_mcl/{{unique_sample_name}}.domains.dmnd.csv", unique_sample_name=lambda wildcards: get_orthofinder_domains_files()),
		
		expand(f"{config["output_dir"]}/searches/broccoli/{{unique_sample_name}}.domains.fasta", unique_sample_name=lambda wildcards: get_broccoli_domains_files()),
		expand(f"{config["output_dir"]}/clusters/searches/broccoli/dmnd_mcl/{{unique_sample_name}}.domains.dmnd.csv", unique_sample_name=lambda wildcards: get_broccoli_domains_files()),

		expand(f"{config["output_dir"]}/clusters/all/OrthoFinder/dmnd_mcl/{{sample}}.dmnd.csv", sample=lambda wildcards: get_fasta_samples_orthofinder()),
		expand(f"{config["output_dir"]}/clusters/all/broccoli/dmnd_mcl/{{sample}}.dmnd.csv", sample=lambda wildcards: get_fasta_samples_broccoli())

# INITIAL ORTHOLOGY

if config["ww_cl"]["broccoli"] == "TRUE":
	print("Running INITIAL ORTHOLOGY with Broccoli")
	checkpoint broccoli:
		input:
			config["input_dir"]
		output:
			directory(f"{config["output_dir"]}/broccoli/Orthogroup_Sequences"),
			f"{config["output_dir"]}/broccoli/dir_step3/orthologous_groups.txt",
			f"{config["output_dir"]}/broccoli/dir_step4/orthologous_pairs.txt"
		params:
			tool_params=config["broccoli"]["broccoli_params"],
			output_dir=config["output_dir"]
		conda:
			"../envs/initial_ortho.yaml"
		shell:
			"""
			rm -fr {params.output_dir}/broccoli &&
			mkdir -p {params.output_dir}/broccoli &&
			python Broccoli/broccoli.py -dir {input} -threads {threads} {params.tool_params} &&
			mv dir_step1 {params.output_dir}/broccoli &&
			mv dir_step2 {params.output_dir}/broccoli &&
			mv dir_step3 {params.output_dir}/broccoli &&
			mv dir_step4 {params.output_dir}/broccoli &&
			python workflow/scripts/parse_fastas_broccoli.py \
				-b {params.output_dir}/broccoli/dir_step3/orthologous_groups.txt \
				-f {input} \
				-o {params.output_dir}/broccoli/Orthogroup_Sequences
			"""

if config["ww_cl"]["orthofinder"] == "TRUE":
	print("Running INITIAL ORTHOLOGY with OrthoFinder")
	checkpoint orthofinder:
		input:
			config["input_dir"]
		output:
			directory(f"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences")
		params:
			tool_params=config["orthofinder"]["orthofinder_params"],
			output_dir=config["output_dir"]
		conda:
			"../envs/initial_ortho.yaml"
		shell:
			"""
			rm -fr {params.output_dir}/OrthoFinder_tmp &&
			rm -fr {params.output_dir}/OrthoFinder &&
			orthofinder -n OrthoFinder -f {input} -o {params.output_dir}/OrthoFinder_tmp -t {threads} {params.tool_params} &&
			mv {params.output_dir}/OrthoFinder_tmp/Results_OrthoFinder {params.output_dir}/OrthoFinder &&
			rm -fr {params.output_dir}/OrthoFinder_tmp
			"""



# example:
#rule intermediate:
#	input:
		#get_fasta_files # to aggregate
#		"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences/{sample}.fa"
#	output:
#		"{config['output_dir']}/OrthoFinder/test/{sample}.fa"
#	shell:
#		"cp {input} {output}"


# SEARCH

## HMM search for orthofinder
if config["ww_cl"]["search"] == "TRUE" and config["ww_cl"]["orthofinder"] == "TRUE":
	print("Running SEARCH with hmmsearch against your hmms of interest for OrthoFinder results")
	rule search_orthofinder:
		input:
			f"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences/{{sample}}.fa"
		output:
			f"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences/{{sample}}.txt"
		params:
			gene_family_info=config["search"]["gene_family_info"],
			gene_family_name=config["search"]["gene_family_name"],
			hmm_dir=config["search"]["hmm_dir"],
			input_source="OrthoFinder",
			output_dir=config["output_dir"]
		conda:
			"../envs/search_cluster.yaml"
		shell:
			"""
			python workflow/scripts/s01_search.py -f {input} -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir} -i {params.input_source} -o {params.output_dir}/searches &&
			touch {output}
			"""
	
	checkpoint search_orthofinder_ckpt:
		input:
			f"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences"
		output:
			directory(f"{config['output_dir']}/searches/OrthoFinder")

else: 
	print("hmm search step SKIPPED for OrthoFinder")

## HMM search for broccoli

if config["ww_cl"]["search"] == "TRUE" and config["ww_cl"]["broccoli"] == "TRUE":
	print("Running SEARCH with hmmsearch against your hmms of interest for broccoli results")
	rule search_broccoli:
		input:
			f"{config['output_dir']}/broccoli/Orthogroup_Sequences/{{sample}}.fa"
		output:
			f"{config['output_dir']}/broccoli/Orthogroup_Sequences/{{sample}}.txt"
		params:
			gene_family_info=config["search"]["gene_family_info"],
			gene_family_name=config["search"]["gene_family_name"],
			hmm_dir=config["search"]["hmm_dir"],
			input_source="broccoli",
			output_dir=config["output_dir"]
		conda:
			"../envs/search_cluster.yaml"
		shell:
			"""
			python workflow/scripts/s01_search.py -f {input} -g {params.gene_family_info} -t {threads} {params.gene_family_name} -H {params.hmm_dir} -i {params.input_source} -o {params.output_dir}/searches &&
			touch {output}
			"""

	checkpoint search_broccoli_ckpt:
		input:
			f"{config['output_dir']}/broccoli/Orthogroup_Sequences"
		output:
			directory(f"{config['output_dir']}/searches/broccoli")

else: 
	print("hmm search step SKIPPED for broccoli")

# CLUSTER

if config["ww_cl"]["diamond_mcl"] == "TRUE":
	print("Running CLUSTER using diamond and MCL")
	if config["ww_cl"]["search"] == "TRUE": # when SEARCH is enabled
		print("Running CLUSTER on outputs of SEARCH")
		if config["ww_cl"]["broccoli"] == "TRUE":
			rule cluster_diamond_mcl_broccoli_search:
				input:
					f"{config['output_dir']}/searches/broccoli/{{unique_sample_name}}.domains.fasta"
				output:
					f"{config['output_dir']}/clusters/searches/broccoli/dmnd_mcl/{{unique_sample_name}}.domains.dmnd.csv"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"],

					mcl_fasta_outdir=f"{config['output_dir']}/clusters/searches/broccoli/dmnd_mcl"
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input} -q {input} -o {output} --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {output} > {output}.tmp &&
					mv {output}.tmp {output} &&
					mcl {output} {params.mcl_params} -I {params.mcl_inflation} -o {output}.abc &&

					python workflow/scripts/parse_fastas_mcl.py \
						-m {output}.abc \
						-f {input} \
						-o {params.mcl_fasta_outdir}
					"""
		if config["ww_cl"]["orthofinder"] == "TRUE":
			rule cluster_diamond_mcl_orthofinder_search:
				input:
					f"{config['output_dir']}/searches/OrthoFinder/{{unique_sample_name}}.domains.fasta"
				output:
					f"{config['output_dir']}/clusters/searches/OrthoFinder/dmnd_mcl/{{unique_sample_name}}.domains.dmnd.csv"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"],

					mcl_fasta_outdir=f"{config['output_dir']}/clusters/searches/OrthoFinder/dmnd_mcl"
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input} -q {input} -o {output} --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {output} > {output}.tmp &&
					mv {output}.tmp {output} &&
					mcl {output} {params.mcl_params} -I {params.mcl_inflation} -o {output}.abc &&

					python workflow/scripts/parse_fastas_mcl.py \
						-m {output}.abc \
						-f {input} \
						-o {params.mcl_fasta_outdir}
					"""
	if config["ww_cl"]["downstream_all"] == "TRUE": # if true, then also run cluster, etc on all the orthogroups
		print("Running CLUSTER on all orthogroups identified")
		if config["ww_cl"]["broccoli"] == "TRUE":
			rule cluster_diamond_mcl_broccoli_all:
				input:
					f"{config['output_dir']}/broccoli/Orthogroup_Sequences/{{sample}}.fa"
				output:
					f"{config['output_dir']}/clusters/all/broccoli/dmnd_mcl/{{sample}}.dmnd.csv"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"],

					mcl_fasta_outdir=f"{config['output_dir']}/clusters/all/broccoli/dmnd_mcl"
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input} -q {input} -o {output} --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {output} > {output}.tmp &&
					mv {output}.tmp {output} &&
					mcl {output} {params.mcl_params} -I {params.mcl_inflation} -o {output}.abc &&

					python workflow/scripts/parse_fastas_mcl.py \
						-m {output}.abc \
						-f {input} \
						-o {params.mcl_fasta_outdir}
					"""
		if config["ww_cl"]["orthofinder"] == "TRUE":
			rule cluster_diamond_mcl_orthofinder_all:
				input:
					f"{config['output_dir']}/OrthoFinder/Orthogroup_Sequences/{{sample}}.fa"
				output:
					f"{config['output_dir']}/clusters/all/OrthoFinder/dmnd_mcl/{{sample}}.dmnd.csv"
				params:
					diamond_params=config["cluster_diamond_mcl"]["diamond_params"],
					mcl_params=config["cluster_diamond_mcl"]["mcl_params"],
					mcl_inflation=config["cluster_diamond_mcl"]["mcl_inflation"],

					mcl_fasta_outdir=f"{config['output_dir']}/clusters/all/OrthoFinder/dmnd_mcl"
				conda:
					"../envs/search_cluster.yaml"
				shell:
					"""
					diamond blastp {params.diamond_params} -d {input} -q {input} -o {output} --threads {threads} &&
					awk '{{ print $1,$2,$12 }}' {output} > {output}.tmp &&
					mv {output}.tmp {output} &&
					mcl {output} {params.mcl_params} -I {params.mcl_inflation} -o {output}.abc &&

					python workflow/scripts/parse_fastas_mcl.py \
						-m {output}.abc \
						-f {input} \
						-o {params.mcl_fasta_outdir}
					"""




