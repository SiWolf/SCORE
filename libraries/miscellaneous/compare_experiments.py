# ------------------------------
# Title: compare_experiments.py
# Author: Silver A. Wolf
# Last Modified: Wed, 31.07.2019
# Version: 0.0.2
# ------------------------------

# Imports
import argparse
import csv
import os

def perform_analysis(e1, e2, path, genes, org_go, org_kegg):
	new_experiment = e1 + "_" + e2
	if not os.path.exists(new_experiment):
		os.makedirs(new_experiment)
	gene_list_e1 = e1 + "/deg_gene_symbols.csv"
	gene_list_e2 = e2 + "/deg_gene_symbols.csv"
	gene_list_filtered = new_experiment + "/deg_gene_symbols.csv"
	first_line = True
	output_list = open(gene_list_filtered, "w")
	output_list.write("gene name\n")
	with open(gene_list_e2) as e2_file:
		for line in csv.reader(e2_file, delimiter = ","):
			if first_line == True:
				first_line = False
			else:
				current_gene_e2 = line[0]
				background_gene = False
				with open(gene_list_e1) as e1_file:
					for line in csv.reader(e1_file, delimiter = ","):
						current_gene_e1 = line[0]
						if current_gene_e1 == current_gene_e2:
							background_gene = True
							break
				if background_gene == False:
					output_list.write(current_gene_e2 + "\n")
	
	gene_scf_folder_go = new_experiment + "/geneSCF_GO/"
	if not os.path.exists(gene_scf_folder_go):
		os.makedirs(gene_scf_folder_go)
	command_go = path + " -m=update -i=" + gene_list_filtered + " -t=sym -o=" + gene_scf_folder_go + " -db=GO_all -p=yes -bg=" + genes + " -org=" + org_go
	os.system(command_go)
	
	gene_scf_folder_kegg = new_experiment + "/geneSCF_KEGG/"
	if not os.path.exists(gene_scf_folder_kegg):
		os.makedirs(gene_scf_folder_kegg)
	command_go = path + " -m=update -i=" + gene_list_filtered + " -t=sym -o=" + gene_scf_folder_kegg + " -db=KEGG -p=yes -bg=" + genes + " -org=" + org_kegg
	os.system(command_go)
	
	output_list.close()

	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-a", "--experiment_1", type = str, default = "V3", required = False, help = "Name of experiment 1")
	parser.add_argument("-b", "--experiment_2", type = str, default = "V4", required = False, help = "Name of experiment 2")
	parser.add_argument("-g", "--gene_scf", type = str, default = "./libraries/geneSCF/geneSCF", required = False, help = "Path to GeneSCF")
	parser.add_argument("-t", "--total_genes", type = str, default = "5000", required = False, help = "Total amount of genes")
	parser.add_argument("-dbg", "--organism_db_go", type = str, default = "ecocyc", required = False, help = "Organism DB (GO)")
	parser.add_argument("-dbk", "--organism_db_kegg", type = str, default = "ecg", required = False, help = "Organism DB (KEGG)")
	args = parser.parse_args()
	
	perform_analysis(args.experiment_1, args.experiment_2, args.gene_scf, args.total_genes, args.organism_db_go, args.organism_db_kegg)