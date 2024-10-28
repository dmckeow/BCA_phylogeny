import sys
import os
import subprocess
import csv
import argparse
import logging

def do_hmmsearch(hmm, out, fasta, threads, threshold, verbose = False):
	logging.info("Running HMMSEARCH")
	logging.info(f"hmm: {hmm}")
	logging.info(f"out: {out}")
	logging.info(f"fasta: {fasta}")
	logging.info(f"threads: {threads}")
	logging.info(f"threshold: {threshold}")

	if not os.path.exists(hmm):
		logging.error(f"HMM model file {hmm} not found!\nCheck the config.yaml file!")
		sys.exit(1)

	if threshold == "GA":
		cmd = (
			f"hmmsearch --domtblout {out}.domtable --cut_ga "
			f"--cpu {threads} {hmm} {fasta} 1> /dev/null"
		)
	else:
		cmd = (
			f"hmmsearch --domtblout {out}.domtable --domE {threshold} "
			f"--cpu {threads} {hmm} {fasta} 1> /dev/null"
		)
	
	subprocess.run(cmd, shell=True, check=True)
	hmmsearch_outfile = f"{out}.domtable.csv.tmp"
	with open(f"{out}.domtable", 'r') as infile, open(hmmsearch_outfile, 'w') as outfile:
		for line in infile:
			if not line.startswith('#'):
				columns = line.split()
				outfile.write('\t'.join([columns[0], columns[17], columns[18], columns[3],columns[4],columns[11]]) + '\n')
	if verbose:
		print(f"Hmmsearch created: {hmmsearch_outfile}")

def merge_results(prefix,gene_family_name, searches_dir, fasta_file, verbose = False):
	#prefix = gene_family_name
	fam_id = gene_family_name
	# Check for any hits
	cmd = f"cat {searches_dir}/{prefix}.{fam_id}.hmmsearch.*.domtable.csv.tmp | grep -v '#' | cut -f 1 | sort -u > {searches_dir}/{prefix}.{fam_id}.genes.list"
	subprocess.run(cmd, shell=True, check=True)
	num_hit_genes_cmd = f"cat {searches_dir}/{prefix}.{fam_id}.genes.list | wc -l"
	num_hit_genes = int(subprocess.check_output(num_hit_genes_cmd, shell=True).strip().decode('utf-8'))
	logging.info(f"# {fasta_file}: {fam_id} | # genes found = {num_hit_genes}")

	if num_hit_genes == 0:
		logging.info(f"# {fasta_file}: {fam_id} | Omit downstream analyses")
	else:
		# Find most-inclusive region that includes all domain hits in the protein
		if verbose: 
			print("-----")
		cmd = f"cat {searches_dir}/{prefix}.{fam_id}.hmmsearch.*.domtable.csv.tmp > {searches_dir}/{prefix}.{fam_id}.domtable.csv.tmp"
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
		cmd = (
			f"bedtools merge -i <(sort -k1,1 -k2,2n {searches_dir}/{prefix}.{fam_id}.domtable.csv.tmp) "
			f"-c 4 -o collapse -d 100000 > {searches_dir}/{prefix}.{fam_id}.domains.csv.tmp2 "
		)
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
		# Extract complete sequences
		cmd = f"samtools faidx {fasta_file} -r {searches_dir}/{prefix}.{fam_id}.genes.list > {searches_dir}/{prefix}.{fam_id}.seqs.fasta"
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, check=True)
		# Dict sequence lengths
		cmd = f"samtools faidx {searches_dir}/{prefix}.{fam_id}.seqs.fasta"
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, check=True)
		
		# Expand domain region by a fixed amount of aa
		cmd = (
			f"bedtools slop -i {searches_dir}/{prefix}.{fam_id}.domains.csv.tmp2 "
			f"-g {searches_dir}/{prefix}.{fam_id}.seqs.fasta.fai -b 50 "
			f"| awk 'BEGIN{{OFS=@\\t@}}{{ print $1, $2+1, $3, $4 }}' > {searches_dir}/{prefix}.{fam_id}.domains.csv "
		)
		cmd = cmd.replace("@",'"')
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
		# Extract domain region
		cmd = (
			f"bedtools getfasta -fi {searches_dir}/{prefix}.{fam_id}.seqs.fasta -bed "
			f"<(awk 'BEGIN{{OFS=@\\t@}}{{ print $1, $2, $3, $1 }}' {searches_dir}/{prefix}.{fam_id}.domains.csv) "
			f"> {searches_dir}/{prefix}.{fam_id}.domains.fasta"
		)
		cmd = cmd.replace("@",'"')
		if verbose:
			print(cmd)
		subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
		# Report
		num_unique_domains_cmd = f"cut -f1 {searches_dir}/{prefix}.{fam_id}.domains.csv | wc -l"
		num_unique_domains = subprocess.check_output(num_unique_domains_cmd, shell=True).strip().decode('utf-8')
		logging.info(f"# {fasta_file}: {fam_id} | # unique domains = {num_unique_domains}")


def parse_gene_family_info(gene_family_info):
	gene_families = {}
	with open(gene_family_info, 'r') as file:
		reader = csv.DictReader(file, delimiter='\t', fieldnames=["class_name", "hmms", "inflation", "min_seqs", "threshold", "pref1", "pref2"])
		for row in reader:
			class_name = row["class_name"]
			gene_families[class_name] = {
				"hmms": row["hmms"].split(','),
				"inflation": row["inflation"],
				"min_seqs": row["min_seqs"],
				"threshold": row["threshold"],
				"pref1": row["pref1"],
				"pref2": row["pref2"]
			}
	return gene_families

def search(fasta_file, gene_family_info, gene_family_name, hmm_dir, threads, output_dir):
	logging.info(f"# {fasta_file}: {gene_family_name} | HMM search")
	gene_families = parse_gene_family_info(gene_family_info)
   # if verbose: 
		#print(gene_families)
	if gene_family_name not in gene_families:
		logging.error(f"Gene family {gene_family_name} not found in gene family info file.")
		return

	gene_family = gene_families[gene_family_name]

	base_name = os.path.basename(os.path.splitext(fasta_file)[0])

    # If you want to handle multiple extensions
	for ext in [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn"]:
		if base_name.endswith(ext):
			base_name = base_name[:-len(ext)]
			break

	prefix = f"{base_name}.{gene_family['pref2']}"
	#prefix = gene_family["pref2"]
	
	hmms = gene_family["hmms"]
	threshold = gene_family["threshold"]

	searches_dir = os.path.join(f"{args.output_dir}/", f"{args.input_source}/")
	os.makedirs(searches_dir, exist_ok=True)
	
	tmp_file = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.domains.csv.tmp")
	with open(tmp_file, 'w') as outfile:
		for hmm in hmms:
			logging.info(f"Running HMM search for {gene_family_name} with HMM: {hmm}")
			output_pref = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.hmmsearch.domain_{hmm}")
			do_hmmsearch(os.path.join(hmm_dir, f"{hmm}.hmm"), output_pref, fasta_file, threads, threshold)
			# Concatenate results
			output_file = f'{output_pref}.domtable'
			with open(output_file, 'r') as infile:
				for line in infile:
					outfile.write(line)
	# Check for any hits
	genes_list_file = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.genes.list")
	cmd = f"grep -v '#' {tmp_file} | cut -f 1 -d ' ' | sort -u > {genes_list_file}"
	subprocess.run(cmd, shell=True, check=True)
	num_hit_genes_cmd = f"cat {genes_list_file} | wc -l"
	num_hit_genes = int(subprocess.check_output(num_hit_genes_cmd, shell=True).strip().decode('utf-8'))
	logging.info(f"# {fasta_file}: {gene_family_name} | # genes found = {num_hit_genes}")

	if num_hit_genes == 0:
		logging.info(f"# {fasta_file}: {gene_family_name} | Omit downstream analyses")
	else:
		merge_results(prefix,gene_family_name, searches_dir, fasta_file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""
	Python wrapper around some useful commands
	""")

	# Search 
	parser.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
	parser.add_argument('-g', '--gene_family_info', required=True, help='Path to the gene family info file specifying HMMs and parameters')
	parser.add_argument('-H', '--hmm', required=True, help='Path to the directory containing hmm profiles against which the search will be performed')
	parser.add_argument('-t', '--threads', required=True, help='Num cpus for hmmsearch')
	parser.add_argument('-i', '--input_source', required=True, help='A meaningful name indicating where your input fastas came from - will be used as subdirectory name in results_annotations/searches')
	parser.add_argument('-o', '--output_dir', required=True, help='Directory for outputs in which a directory <input_source> will be created')
	parser.add_argument('gene_family_name', help='Name of the gene family to search')

	args = parser.parse_args()

	verbose = True

	search(args.fasta, args.gene_family_info, args.gene_family_name, args.hmm, args.threads, args.output_dir)
