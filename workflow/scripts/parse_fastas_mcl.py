import os
import argparse
from Bio import SeqIO

def mcl_input(input_file, output_file):
	# Extract the base file name without extension
	base_filename = input_file.split('/')[-1].replace('.dmnd.csv.abc', '')
	
	# Initialize an empty list to store the modified lines
	modified_lines = []
	
	# Add the header line
	modified_lines.append("#placeholder\theader")  # Modify headers as needed
	
	# Read the input file and modify the lines
	with open(input_file, 'r') as file:
		for index, line in enumerate(file):
			# Replace all tabs with spaces
			line = line.replace('\t', ' ')
			# Create the new identifier
			new_id = f"{base_filename}.mcl.{index + 1}"
			# Combine the new identifier with the modified line
			modified_line = f"{new_id}\t{line.strip()}"
			modified_lines.append(modified_line)
	
	# Write the modified lines to the output file
	with open(output_file, 'w') as file:
		for modified_line in modified_lines:
			file.write(modified_line + '\n')


def extract_fasta_sequences(og_file, fasta_file, output_dir):
	# Read the OG to protein names mapping from the input file
	with open(og_file, "r") as og:
		for line in og.readlines()[1:]:  # Skip the header
			og_name, proteins = line.strip().split("\t")
			protein_list = proteins.split(" ")  # Get list of protein names

			# Create a new fasta file for the orthogroup
			output_file = os.path.join(output_dir, f"{og_name}.fa")
			
			with open(output_file, "w") as out_fasta:
				# Parse the single fasta file for matching records
				for record in SeqIO.parse(fasta_file, "fasta"):
					if record.id in protein_list:
						SeqIO.write(record, out_fasta, "fasta")

	print("parse_fastas_mcl.py: Fasta extraction complete.")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-f','--fasta', required=True, help='Path to the input fasta file containing sequences')
	parser.add_argument('-m', '--mcl_ogs', required=True, help='Path to the .abc cluster file from MCL')
	parser.add_argument('-o', '--outdir', required=True, help='Path to the directory where mcl orthogroup fastas will be put')
	
	args = parser.parse_args()

	os.makedirs(args.outdir, exist_ok=True)

	# use MCL .abc as input to make fasta files for lower granularity clusters within an OG
	output = f"{args.mcl_ogs}.mcl.tmp"
	mcl_input(args.mcl_ogs, output)
	extract_fasta_sequences(output, args.fasta, args.outdir)
