# --------------------------------
# Title: generate_transcriptome.py
# Author: Silver A. Wolf
# Last Modified: Wed, 18.12.2019
# Version: 0.0.2
# --------------------------------

import argparse
import csv
import fasta_wrapper
import gffutils
import pyfaidx

def generate_transcriptome(ANNOT, fasta, gff, ID, ffn):
	errors = FALSE
	ID_list = []
	# First checks for duplicate entries with ANNOT/ID pair
	with open(gff) as gff_file:
		for line in csv.reader(gff_file, delimiter = "\t"):
			if line[0][0] != "#":
				if line[2] == ANNOT:
					current_id = line[8].split(ID + "=")[1].split(";")[0].strip()
					if current_id in ID_list:
						errors = TRUE
						print("Error - Duplicate ID - " + current_id)
					else:
						ID_list.append(current_id)

	db = gffutils.create_db(gff, "tmp.db", merge_strategy = "create_unique")
	fasta = pyfaidx.Fasta(fasta)
	results = open(ffn, "w")

	for cds in db.features_of_type(ANNOT, order_by = "start"):
		gff_id = ''.join(cds[ID]).strip()
		fasta_sequence = cds.sequence(fasta)
		clean_sequence = fasta_wrapper.fasta_wrapper(fasta_sequence, 60)
		results.write(">" + gff_id + "\n" + clean_sequence + "\n")
	
	results.close()

# Main
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-a", "--annotation_feature", type = str, default = "gene", required = False, help = "Feature in gff to consider for the analysis")
	parser.add_argument("-f", "--fasta_genome", type = str, default = "references/REF.fasta", required = False, help = "Reference genome in fasta format")
	parser.add_argument("-g", "--gff_annotation_file", type = str, default = "references/REF.gff", required = False, help = "Reference annotation file in gff format")
	parser.add_argument("-i", "--identifier", type = str, default = "locus_tag", required = False, help = "Identifier for individual transcripts")
	parser.add_argument("-t", "--fasta_transcriptome", type = str, default = "references/REF.ffn", required = False, help = "Reference transcriptome in fasta format")
	args = parser.parse_args()

	generate_transcriptome(args.annotation_feature, args.fasta_genome, args.gff_annotation_file, args.identifier, args.fasta_transcriptome)