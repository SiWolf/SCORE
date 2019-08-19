# -------------------------------
# Title: merge_tigrfam_results.py
# Author: Silver A. Wolf
# Last Modified: Mo, 19.08.2019
# Version: 0.0.3
# -------------------------------

# Imports
import argparse
import csv
		
def refine_summary_file(tigrfams_links_file, tigrfams_roles_file):
	summary_file_temp = open("deg/summary_extended.tsv", "w")
	with open("deg/summary.tsv") as summary_file:
		for summary_line in csv.reader(summary_file, delimiter = "\t"):
			id = summary_line[0]
			gene = summary_line[1]
			product = summary_line[2]
			fold_change = summary_line[3]
			deg = summary_line[4]
			p_value_bayseq = summary_line[5]
			p_value_deseq2 = summary_line[6]
			p_value_edger = summary_line[7]
			p_value_limma = summary_line[8]
			p_value_noiseq = summary_line[9]
			p_value_sleuth = summary_line[10]
			group_01 = summary_line[11]
			group_02 = summary_line[12]
			nucleotide_sequence = summary_line[13]
			aa_sequence = summary_line[14]
			if id == "ID":
				summary_file_temp.write(id + "\t" + gene + "\t" + product + "\t" + fold_change + "\t" + deg + "\t" + p_value_bayseq + "\t" + p_value_deseq2 + "\t" + p_value_edger + "\t" + p_value_limma + "\t" + p_value_noiseq + "\t" + p_value_sleuth + "\t" + group_01 + "\t" + group_02 + "\t" + "TIGRFAM main role" + "\t" + "TIGRFAM sub role" + "\t" + "TIGRFAM description" + "\t" + "TIGRFAM ID" + "\t" + nucleotide_sequence + "\t" + aa_sequence + "\n")
			else:
				tigrfam_main_role = ""
				tigrfam_sub_role = ""
				tigrfam_description = ""
				tigrfam_id = ""
				with open("deg/hmmer_output.txt") as hmmer_file:
					for hmmer_line in hmmer_file:
						if hmmer_line[0] != "#":
							if id in hmmer_line:
								tigrfam_description = hmmer_line.split(": ")[1]
								tigrfam_id = hmmer_line.split(" ")[0]
								break
				with open(tigrfams_links_file) as link_file:
					for link_line in link_file:
						if tigrfam_id in link_line:
							role_id = link_line.split("\t")[1]
							break
				with open(tigrfams_roles_file) as roles_file:
					for roles_line in roles_file:
						current_role = roles_line.split("\t")[1]
						if role_id == current_role:
							role_type = roles_line.split("\t")[2]
							if role_type == "mainrole:":
								tigrfam_main_role = roles_line.split("\t")[3]
							else:
								tigrfam_sub_role = roles_line.split("\t")[3]
				summary_file_temp.write(id + "\t" + gene + "\t" + product + "\t" + fold_change + "\t" + deg + "\t" + p_value_bayseq + "\t" + p_value_deseq2 + "\t" + p_value_edger + "\t" + p_value_limma + "\t" + p_value_noiseq + "\t" + p_value_sleuth + "\t" + group_01 + "\t" + group_02 + "\t" + tigrfam_main_role + "\t" + tigrfam_sub_role + "\t" + tigrfam_description + "\t" + tigrfam_id + "\t" + nucleotide_sequence + "\t" + aa_sequence + "\n")

	summary_file_temp.close()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-l", "--links_file", type = str, default = "references/tigrfam/TIGRFAMS_ROLE_LINK", required = False, help = "Links file for TIGRFAM")
	parser.add_argument("-r", "--roles_file", type = str, default = "references/tigrfam/TIGRFAMS_ROLE_NAMES", required = False, help = "Roles file for TIGRFAM")
	args = parser.parse_args()
	
	refine_summary_file(args.links_file, args.roles_file)