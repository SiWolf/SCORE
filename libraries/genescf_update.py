# ------------------------------
# Title: genescf_update.py
# Author: Silver A. Wolf
# Last Modified: Thue, 05.08.2019
# Version: 0.0.3
# ------------------------------

# Imports
import argparse
import csv

def rewrite_genescf_results(database, obo_file):
	file = "deg/pathway_analysis_GO/deg_gene_symbols_GO_all_" + database + "_functional_classification.tsv"
	first_line = True
	output = open("deg/pathway_analysis_GO/" + file.split(".tsv")[0] + "_updated.tsv", "w")
	with open(file) as original_file:
		for line in csv.reader(original_file, delimiter = "\t"):
			if first_line == True:
				output.write("category\tprocess ID\tdescription\tgenes\tcorrected_pval\n")
				first_line = False
			else:
				genes = line[0]
				process_description = line[1].split("~")[1]
				process_id = line[1].split("~")[0]
				corrected_p_val = line [6]
				current_process_id = "id: " + line[1].split("~")[0]
				category = ""
				upcoming_namespace = False
				with open(obo_file) as go_file:
					for line in go_file:
						if current_process_id in line:
							upcoming_namespace = True
						else:
							if upcoming_namespace == True:
								if "namespace" in line:
									category = line.split(" ")[1].strip()
									break
				output.write(category + "\t" + process_id + "\t" + process_description + "\t" + genes + "\t" + corrected_p_val + "\n")
	output.close()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-d", "--genescf_go_database", type = str, default = "ecocyc", required = False, help = "GO Database used for GeneSCF")
	parser.add_argument("-o", "--go_reference_file", type = str, default = "references/go.obo", required = False, help = "GO Reference File")
	args = parser.parse_args()
	
	rewrite_genescf_results(args.genescf_go_database, args.go_reference_file)