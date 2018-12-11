# Script to simulate RNA-Seq data using polyester
# Will be used for Benchmarking SCORE

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("polyester")

library("Biostrings")
library("polyester")

# Read fasta file
reference_path = "simulation_data/PROKKA_07132018.ffn"
fasta = readDNAStringSet(reference_path)

# Value for coverage
readspertx = round (20 * width(fasta) / 100)

# Matrix of fold changes to be simulated
# Default: First 50 genes overexpressed in Group 1, Second 50 genes overexpressed in Group 2
# TO-DO: Randomize transcript selection
gene_set = 50
gene_position_first_group = length(fasta) - gene_set
gene_position_second_group = length(fasta) - gene_set - gene_set
fold_changes = matrix(c(rep(4, gene_set), rep(1, gene_position_first_group), rep(1, gene_set), rep(4, gene_set), rep(1, gene_position_second_group)), nrow = length(fasta))

# Run simulation
replicates = c(3, 3)
results_folder = "simulation_data/"
simulate_experiment(reference_path, reads_per_transcript = readspertx, num_reps = replicates, fold_changes = fold_changes, outdir = results_folder, readlen = 100, paired = TRUE)