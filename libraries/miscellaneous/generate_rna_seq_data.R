# Script to simulate RNA-Seq data using Polyester
# Useful for benchmarking SCORE

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

overrepresented_genes = 50
random_overrepresented_genes = FALSE

# Matrix of fold changes to be simulated
if (random_overrepresented_genes == FALSE){
  # Default: First 50 genes overexpressed in Group 1, Second 50 genes overexpressed in Group 2
  gene_position_first_group = length(fasta) - overrepresented_genes
  gene_position_second_group = length(fasta) - overrepresented_genes - overrepresented_genes
  fold_changes = matrix(c(rep(4, overrepresented_genes), rep(1, gene_position_first_group), rep(1, overrepresented_genes), rep(4, overrepresented_genes), rep(1, gene_position_second_group)), nrow = length(fasta))
} else {
  # Randomize transcript selection, fixed amount of 50 per group
  random_lists <- sample(1:length(fasta), 2*overrepresented_genes, replace = F)
  random_group_1 <- random_lists[0:50]
  random_group_2 <- random_lists[51:100]
  group_1 = ""
  group_2 = ""
  for (i in 0:length(fasta)){
    if (i %in% random_group_1){
      group_1[i] = 4
    } else {
      group_1[i] = 1
    }
    if (i %in% random_group_2){
      group_2[i] = 4
    } else {
      group_2[i] = 1
    }
  }
  fold_changes = matrix(c(group_1, group_2), nrow = length(fasta))
}

# Run simulation
replicates = c(3, 3)
results_folder = "simulation_data/"
simulate_experiment(reference_path, reads_per_transcript = readspertx, num_reps = replicates, fold_changes = fold_changes, outdir = results_folder, readlen = 100, paired = TRUE)