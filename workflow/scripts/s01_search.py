import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml

# Ensure PyYAML is installed
try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is not installed. Install it using 'pip install pyyaml'.")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config():
    config_file = "config.yaml"
    tool_directory = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(tool_directory, config_file)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file {config_path} not found.")
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def check(tool_name):
    try:
        # Try to run the tool with a harmless argument like --help
        cmd = [tool_name, '--help']
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # If return code is 0, it means the tool executed successfully
        if result.returncode == 0:
            print(f"{tool_name} is available and can be launched.")
        else:
            print(f"Error: {tool_name} is not functioning properly.")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: {tool_name} is not available on your system.")
        sys.exit(1)



def activate_conda_environment(env_name):
    logging.info(f"Activating conda environment: {env_name}")
    subprocess.run(f"conda activate {env_name}", shell=True, check=True)

def do_hmmsearch(hmm, out, fasta, cpu, threshold, verbose = False):
    logging.info("Running HMMSEARCH")
    logging.info(f"hmm: {hmm}")
    logging.info(f"out: {out}")
    logging.info(f"fasta: {fasta}")
    logging.info(f"cpu: {cpu}")
    logging.info(f"threshold: {threshold}")

    if not os.path.exists(hmm):
        logging.error(f"HMM model file {hmm} not found!\nCheck the config.yaml file!")
        sys.exit(1)

    if threshold == "GA":
        cmd = (
            f"hmmsearch --domtblout {out}.domtable --cut_ga "
            f"--cpu {cpu} {hmm} {fasta} 1> /dev/null"
        )
    else:
        cmd = (
            f"hmmsearch --domtblout {out}.domtable --domE {threshold} "
            f"--cpu {cpu} {hmm} {fasta} 1> /dev/null"
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

def search(fasta_file, gene_family_info, gene_family_name, config, verbose):
    logging.info(f"# {fasta_file}: {gene_family_name} | HMM search")
    gene_families = parse_gene_family_info(gene_family_info)
   # if verbose: 
        #print(gene_families)
    if gene_family_name not in gene_families:
        logging.error(f"Gene family {gene_family_name} not found in gene family info file.")
        return

    gene_family = gene_families[gene_family_name]
    prefix = gene_family["pref2"]
    hmms = gene_family["hmms"]
    threshold = gene_family["threshold"]
    cpu = config.get("cpu", 4)
    hmm_dir = config["hmm_dir"]
    searches_dir = 'results_annotation/searches/'
    os.makedirs(searches_dir, exist_ok=True)
    
    tmp_file = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.domains.csv.tmp")
    with open(tmp_file, 'w') as outfile:
        for hmm in hmms:
            logging.info(f"Running HMM search for {gene_family_name} with HMM: {hmm}")
            output_pref = os.path.join(searches_dir, f"{prefix}.{gene_family_name}.hmmsearch.domain_{hmm}")
            do_hmmsearch(os.path.join(hmm_dir, f"{hmm}.hmm"), output_pref, fasta_file, cpu, threshold)
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
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # Search 
    parser_search = subparsers.add_parser('search', help='Search for a family using HMMER')
    parser_search.add_argument('-f','--fasta', required=True, help='Path to the input fasta file')
    parser_search.add_argument('-g', '--gene_family_info', required=True, help='Path to the gene family info file specifying HMMs and parameters')
    parser_search.add_argument('gene_family_name', help='Name of the gene family to search')

    args = parser.parse_args()

    verbose = True

    config = load_config()

    if args.command == 'search':
        logging.info("Command: Search")
        search(args.fasta, args.gene_family_info, args.gene_family_name, config, verbose)
    else:
        parser.print_help()
        print("main.py now only used for search or clustering")
        sys.exit(1)