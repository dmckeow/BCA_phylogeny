configfile: "config.yaml"

rule s02_cluster:
	input:
		fasta=config["s02_cluster"]["fasta_file"]
	output:
		clusters=config["s02_cluster"]["output_file"]
	params:
		inflation=config["s02_cluster"]["inflation"]
	log:
		"logs/s02_cluster/test.log"
	benchmark:
		"benchmarks/s02_cluster/test.benchmark.txt"
	conda:
	  "envs/s02_cluster.yaml"
	shell:
		"""
		python helper/s02_cluster.py {input.fasta} {output.clusters} {params.inflation}
		"""





