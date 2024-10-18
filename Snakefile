## Main Snakefile - not currently in use
#
## configfile: "config.yaml"
#
## Define the final output for the main workflow
#rule all:
#	input:
#		expand("results_annotation/alignments/{sample}.aln", sample=config["samples"].keys()),
#		expand("results_annotation/alignments/{sample}.aln.clipkit", sample=config["samples"].keys()),
#		expand("results_annotation/searches/{gene_family_name}.domains.fasta", gene_family_name=config["gene_families"])
#
## Subworkflows
#subworkflow search:
#	snakefile: "rules/search.smk"
#	configfile: "config.yaml"
#
#subworkflow align_phylo:
#	snakefile: "rules/align_phylo.smk"
#	configfile: "config/align_phylo.yaml"
#
#
#
#