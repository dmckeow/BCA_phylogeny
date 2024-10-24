configfile: "config.yaml"

rule all:
	input:
		"plots/quals.svg"

def get_bwa_map_input_fastqs(wildcards):
	return config["samples"][wildcards.sample]

rule bwa_map:
	input:
		"data/genome.fa",
		get_bwa_map_input_fastqs
	output:
		temp("mapped_reads/{sample}.bam")
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}"
	log:
		"logs/bwa_map/{sample}.log"
	benchmark:
		"benchmarks/{sample}.bwa_map.benchmark.txt"
	conda:
	  "envs/snakemake_tutorial.yaml"
	shell:
		"(bwa mem -R '{params.rg}' -t {threads} {input} | "
		"samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		protected("sorted_reads/{sample}.bam")
	conda:
	  "envs/snakemake_tutorial.yaml"
	shell:
		"samtools sort -T sorted_reads/{wildcards.sample} "
		"-O bam {input} > {output}"

rule samtools_index:
  input:
	  "sorted_reads/{sample}.bam"
  output:
	  "sorted_reads/{sample}.bam.bai"
  conda:
	  "envs/snakemake_tutorial.yaml"
  shell:
	  "samtools index {input}"


# aggregate outputs from all mappings then jointly run variant calling on them
rule bcftools_call:
	input:
		fa="data/genome.fa",
		bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
		bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
	output:
		"calls/all.vcf"
	params:
		mr=config["mutation_rate"]
	log:
		"logs/bcftools_call/bcftools_call.log"
	conda:
	  "envs/snakemake_tutorial.yaml"
	shell:
		"bcftools mpileup -f {input.fa} {input.bam} | "
		"(bcftools call -mv -P '{params.mr}' - > {output}) 2> {log}"

# use a custom script
rule plot_quals:
	input:
		"calls/all.vcf"
	output:
		"plots/quals.svg"
	log:
		"logs/bcftools_call/plot_quals.log"
	conda:
	  "envs/snakemake_tutorial.yaml"
	script:
		"scripts/plot-quals.py"

