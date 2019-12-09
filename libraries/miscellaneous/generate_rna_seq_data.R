# --------------------------------------------
# Title: generate_rna_seq_data.R
# Author: Silver A. Wolf
# Last Modified: Fr, 06.12.2019
# Version: 0.0.2
# --------------------------------------------

# Script to simulate RNA-Seq data using Polyester
# Used for benchmarking SCORE

# Installers
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("polyester")

# Imports
library("Biostrings")
library("polyester")

# Main
args <- commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]
argument_3 = args[3]
argument_4 = args[4]
argument_5 = args[5]
argument_6 = args[6]

# Special case if this script is executed manually without any given parameters
# Example: RStudio
if (is.na(argument_1)){
  argument_1 = "references/BENCHMARK_TRANSCRIPTOME.ffn"
  argument_2 = 20
  argument_3 = 4
  argument_4 = 50
  argument_5 = TRUE
  argument_6 = 100
}

reference_path = argument_1
coverage_val = as.numeric(argument_2)
deg_fold_change = as.numeric(argument_3)
overrepresented_genes = as.numeric(argument_4)
random_overrepresented_genes = as.logical(argument_5)
read_length = as.numeric(argument_6)

# Read fasta file
fasta = readDNAStringSet(reference_path)

# Value for coverage
readspertx = round(coverage_val * width(fasta) / read_length)

# Matrix of fold changes to be simulated
if (random_overrepresented_genes == FALSE){
  # Default: First 50 genes overexpressed in Group 1, Second 50 genes overexpressed in Group 2
  gene_position_first_group = length(fasta) - overrepresented_genes
  gene_position_second_group = length(fasta) - overrepresented_genes - overrepresented_genes
  fold_changes = matrix(c(rep(deg_fold_change, overrepresented_genes), rep(1, gene_position_first_group), rep(1, overrepresented_genes), rep(deg_fold_change, overrepresented_genes), rep(1, gene_position_second_group)), nrow = length(fasta))
} else {
  # Randomize transcript selection, fixed amount of 50 per group
  random_lists <- sample(1:length(fasta), 2*overrepresented_genes, replace = F)
  random_group_1 <- random_lists[0:overrepresented_genes]
  random_group_2 <- random_lists[(overrepresented_genes+1):(2*overrepresented_genes)]
  group_1 = ""
  group_2 = ""
  for (i in 0:length(fasta)){
    if (i %in% random_group_1){
      group_1[i] = deg_fold_change
    } else {
      group_1[i] = 1
    }
    if (i %in% random_group_2){
      group_2[i] = deg_fold_change
    } else {
      group_2[i] = 1
    }
  }
  fold_changes = matrix(c(as.numeric(group_1), as.numeric(group_2)), nrow = length(fasta))
}

# Run simulation
replicates = c(3, 3)
results_folder = "libraries/miscellaneous/simulation_data/"
simulate_experiment(reference_path, reads_per_transcript = readspertx, num_reps = replicates, fold_changes = fold_changes, outdir = results_folder, readlen = read_length, paired = TRUE)
system("gzip libraries/miscellaneous/simulation_data/*.fasta")