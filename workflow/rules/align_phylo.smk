configfile: "config/align_phylo.yaml"

def get_input_fastas(wildcards):
	return config["samples"][wildcards.sample]

rule all:
	input:
		expand("results_annotation/alignments/{sample}.aln", sample=config["samples"].keys()),
		expand("results_annotation/alignments/{sample}.aln.clipkit", sample=config["samples"].keys())

rule align:
	input:
		get_input_fastas
	output:
		"results_annotation/alignments/{sample}.aln"
	params:
		tool_params=config["align"]["mafft_params"]
	conda:
		"../envs/align_phylo.yaml"
	shell:
		"""
		mafft --reorder --thread {threads} {params.tool_params} {input} > {output}
		"""

rule trim:
	input:
		"results_annotation/alignments/{sample}.aln"
	output:
		"results_annotation/alignments/{sample}.aln.clipkit"
	params:
		tool_params=config["trim"]["clipkit_params"]
	conda:
		"../envs/align_phylo.yaml"
	shell:
		"""
		clipkit {input} {params.tool_params} -o {output}
		"""


