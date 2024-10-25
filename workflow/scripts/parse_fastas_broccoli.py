import os
import argparse
from Bio import SeqIO


def extract_fasta_sequences(og_file, fasta_dir, output_dir):
	# Read the OG to protein names mapping from the input file
	with open(og_file, "r") as og:
		for line in og.readlines()[1:]:  # Skip the header
			og_name, proteins = line.strip().split("\t")
			protein_list = proteins.split(" ")  # Get list of protein names

			# Create a new fasta file for the orthogroup
			output_file = os.path.join(output_dir, f"{og_name}.fa")
			
			with open(output_file, "w") as out_fasta:
				# Loop through each fasta file in the fasta_dir
				for fasta_file in os.listdir(fasta_dir):
					if fasta_file.endswith(".fa") or fasta_file.endswith(".fasta"):
						fasta_path = os.path.join(fasta_dir, fasta_file)
						# search for matching fastas between orthogroup list and input fastas in dir
						for record in SeqIO.parse(fasta_path, "fasta"):
							if record.id in protein_list:
								SeqIO.write(record, out_fasta, "fasta")

	print("parse_fastas_broccoli.py: Fasta extraction complete.")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-f','--fasta', required=True, help='Path to the input fasta file directory containing .fa or .fasta files')
	parser.add_argument('-b', '--broccoli_ogs', required=True, help='Path to the Broccoli output dir_step3/orthologous_groups.txt')
	parser.add_argument('-o', '--outdir', required=True, help='Path to the directory where broccoli orthogroup fastas will be put')
	
	args = parser.parse_args()

	os.makedirs(args.outdir, exist_ok=True)

	# Just extract using Broccoli's dir_step3/orthologous_groups.txt to make fastas for OGs which Broccoli should ALREADY FUCKING DO
	extract_fasta_sequences(args.broccoli_ogs, args.fasta, args.outdir)
