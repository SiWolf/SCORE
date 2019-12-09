# ----------------------------------
# Title: preprocess_transcriptome.py
# Author: Silver A. Wolf
# Last Modified: Mo, 09.12.2019
# Version: 0.0.2
# ----------------------------------

# Imports
import argparse
import os

def update_headers(input_file, id):
	found_id = False
	output = open(input_file + ".tmp", "w")
	seq = ""

	with open(input_file) as file:
		for line in file:
			if id in line:
				name = "> " + line.split(id + "=")[1].split("]")[0] + "\n"
				output.write(name)
				found_id = True
			elif found_id == True:
				seq = line.strip() + "\n"
				output.write(seq)
			else:
				break

	output.close()

	if found_id == False:
		os.system("cp " + input_file + " " + input_file + ".tmp")
		os.system("sed -i 's/ /|/g' " + input_file + ".tmp")

# Main
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-f", "--transcriptome_fasta", type = str, default = "references/REF.ffn", required = False, help = "Transcriptome fasta to be updated")
	parser.add_argument("-i", "--transcript_identifier", type = str, default = "locus_tag", required = False, help = "Identifier for individual transcripts")
	args = parser.parse_args()

	update_headers(args.transcriptome_fasta, args.transcript_identifier)
