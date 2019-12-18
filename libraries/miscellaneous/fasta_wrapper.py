# -------------------------------
# Title: fasta_wrapper.py
# Author: Silver A. Wolf
# Last Modified: Wed, 18.12.2019
# Version: 0.0.1
# -------------------------------

# Library for sequence formatting
# Inputs: sequence as string, width as int
# Will break sequence at every width using a newline

def fasta_wrapper(sequence, width):
	sequence = sequence.strip()
	sequence_count = 0
	sequence_length = len(sequence)
	new_sequence = ""
	new_sequence_count = 0
	width = width - 1
	for base in sequence:
		sequence_count += 1
		if (new_sequence_count == width) and (sequence_length - sequence_count > 0):
			new_sequence = new_sequence + base + "\n"
			new_sequence_count = 0
		else:
			new_sequence = new_sequence + base
			new_sequence_count += 1
	return(new_sequence)