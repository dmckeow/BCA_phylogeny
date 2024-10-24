import sys
import os
import subprocess
import csv
import argparse
import logging
import yaml

def cluster(fasta_file, output_file, inflation, ncpu = 5, max_target_seqs = 100, clean = True):
    if not os.path.exists(fasta_file):
        logging.error(f"{fasta_file} doesn't exist")
        sys.exit(1)

    logging.info(f"Clustering sequences {fasta_file} with diamond & MCL {inflation}")

    cmd = (f"diamond makedb --in {fasta_file} -d {fasta_file} --quiet")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    output_csv = output_file + '.csv'
    cmd = (f"diamond blastp  --more-sensitive --max-target-seqs {max_target_seqs}  -d {fasta_file}  -q {fasta_file} -o {output_csv} --quiet --threads {ncpu}")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    output_abc = output_file + '.abc'
    cmd = (f"awk '{{ print $1,$2,$12 }}' {output_csv} > {output_abc}")
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # run the mcl clustering
    cmd = f"mcl {output_abc} --abc -I {inflation} -o {output_file} 2> /dev/null"
    logging.info(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # remove intermediate files
    if clean:
        cmd = f"rm {output_abc}; rm {output_csv}"
        logging.info(cmd)
        subprocess.run(cmd, shell=True, check=True)
    logging.info(f"Clustering done: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster sequences using DIAMOND and MCL.")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("outfile", help="Output file prefix")
    parser.add_argument("inflation", type=float, help="Inflation parameter for MCL")
    parser.add_argument("--ncpu", type=int, default=5, help="Number of CPU threads")
    
    args = parser.parse_args()
    cluster(args.fasta, args.outfile, args.inflation, args.ncpu)
    
