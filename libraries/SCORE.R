# --------------------------------------------
# Title: SCORE.R
# Author: Silver A. Wolf
# Last Modified: Thur, 25.07.2019
# Version: 0.6.1
# --------------------------------------------

# Installers
#source("https://bioconductor.org/biocLite.R")
#biocLite("baySeq")
#biocLite("DESeq2")
#biocLite("edgeR")  
#biocLite("limma")
#biocLite("NOISeq")
#biocLite("rhdf5")
#install.packages("devtools")
#install.packages("UpSetR")
#devtools::install_github("pachterlab/sleuth")

# Imports
library("baySeq")
library("DESeq2")
library("edgeR")
library("limma")
library("NOISeq")
library("rhdf5")
library("sleuth")
library("UpSetR")

# Functions

# Calculates FN, FP, TN and TP ratios
# Uses these to compute ACC, TNR, TPR, FDR, FPR, FNR and PRE values
calculate_statistics <- function(binaries, consensus){
  total_binaries = as.numeric(length(rownames(binaries)))
  total_tools = as.numeric(length(colnames(binaries)))
  to = 1
  fn_vector = c()
  fp_vector = c()
  tn_vector = c()
  tp_vector = c()
  acc_vector = c()
  tnr_vector = c()
  tpr_vector = c()
  fdr_vector = c()
  fpr_vector = c()
  fnr_vector = c()
  pre_vector = c()
  for (tool in colnames(binaries)){
    fn_entry = 0
    fp_entry = 0
    tn_entry = 0
    tp_entry = 0
    tr = 1
    if (to == 7){
      for (transcript in rownames(binaries)){
        current_truth <- binaries[tr, total_tools]
        if (current_truth == TRUE){
          if (transcript %in% rownames(consensus)){
            tp_entry = tp_entry + 1
          } else {
            fn_entry = fn_entry + 1
          }
        } else {
          if (transcript %in% rownames(consensus)){
            fp_entry = fp_entry + 1
          } else {
            tn_entry = tn_entry + 1
          }
        }
        tr = tr + 1
      }
    } else {
      for (transcript in rownames(binaries)){
        current_frame <- binaries[tr, to]
        current_truth <- binaries[tr, total_tools]
        if (as.logical(current_frame) == current_truth){
          if (as.logical(current_frame) == FALSE){
            tn_entry = tn_entry + 1
          } else {
            tp_entry = tp_entry + 1
          }
        }
        if (as.logical(current_frame) != current_truth){
          if (current_truth == FALSE){
            fp_entry = fp_entry + 1
          } else {
            fn_entry = fn_entry + 1
          }
        }
        tr = tr + 1
      }
    }
    fn_vector[to] <- fn_entry
    fp_vector[to] <- fp_entry
    tn_vector[to] <- tn_entry
    tp_vector[to] <- tp_entry
    acc_vector[to] <- (tp_entry + tn_entry)/(tp_entry + tn_entry + fp_entry + fn_entry)
    tnr_vector[to] <- tn_entry/(fp_entry + tn_entry)
    tpr_vector[to] <- tp_entry/(tp_entry + fn_entry)
    fdr_vector[to] <- fp_entry/(fp_entry + tp_entry)
    fpr_vector[to] <- fp_entry/(fp_entry + tn_entry)
    fnr_vector[to] <- fn_entry/(fn_entry + tp_entry)
    pre_vector[to] <- tp_entry/(tp_entry + fp_entry)
    to = to + 1
  }
  statistics <- data.frame(FN = fn_vector, FP = fp_vector, TN = tn_vector, TP = tp_vector, TPR = tpr_vector, TNR = tnr_vector, ACC = acc_vector, FDR = fdr_vector, FPR = fpr_vector, FNR = fnr_vector, PRE = pre_vector)
  return(statistics)
}

# Estimates TPM values from raw counts
calculate_tpm <- function(transcript_counts){
  transcript_lengths <- read.table(file = "transcript_lengths.csv", sep = ",", header = TRUE)
  c = 0
  for (transcript in row.names(transcript_counts)){
    l = transcript_lengths[transcript_lengths$Transcript.ID == transcript, 2] / 1000
    transcript_counts[c,] <- transcript_counts[c,] / l
    c = c + 1
  }
  d = 1
  for (sample in colnames(transcript_counts)){
    sum_column = sum(transcript_counts[,d])
    tpm_scaling_factor = sum_column / 1000000
    transcript_counts[,d] <- transcript_counts[,d] / tpm_scaling_factor
    d = d + 1
  }
  return(transcript_counts)
}

# Reads the first counts-file in order to fetch a list of all gene symbols
create_gene_list <- function(sample){
  sample_path <- paste("mapped/bowtie2/featureCounts/", sample, sep = "")
  setwd(sample_path)
  gene_list = read.csv(paste("counts", sample, sep = "_"), sep = "", head = T, skip = 1)[,c("Geneid")]
  return(gene_list)
}

# Merges the individual count-files into a matrix
# Applies a low-expression cutoff in order to reduce computational time
create_count_matrix <- function(sample_list, gene_list, low_expression_cutoff){
  sample_nr = 0
  # Let's receive the rest of the counts...
  for (sample in sample_list){
    sample_nr = sample_nr + 1
    path <- paste("../", sample, sep = "")
    setwd(path)
    counts = read.csv(paste("counts", sample, sep = "_"), sep = "", head = T, skip = 1, row.names = 1)
    # Remove the first six columns (genesymbol, chr, start, end, strand, length)
    counts <- counts[ ,6:ncol(counts)]
    if (sample_nr == 1){
      count_matrix <- counts
    } else {
      count_matrix <- cbind(count_matrix, counts)  
    }
  }
  # Convert to a matrix
  count_matrix <- as.matrix(count_matrix)
  # Rename the columns
  colnames(count_matrix) <- paste(sample_list)
  rownames(count_matrix) <- paste(gene_list)
  # unfiltered_count_matrix <- count_matrix
  
  # Filter lowly expressed genes
  count_matrix <- filter_matrix(count_matrix, low_expression_cutoff)
  
  # Initial visualization
  # Histogram
  cpm_log <- cpm(count_matrix, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  hist(median_log2_cpm)
  abline(v = low_expression_cutoff, col = "red", lwd = 3)
  sum(median_log2_cpm > low_expression_cutoff)
  # Principal component analysis (PCA)
  cpm_log_filtered <- cpm(count_matrix, log = TRUE)
  pca <- prcomp(t(cpm_log_filtered), scale. = TRUE)
  plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
  text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log_filtered))
  summary(pca)
  
  return(count_matrix)
}

# Removes genes with low counts
# According to a cutoff for low-expressed genes
# Filter out all genes which do not have rowsum >= cutoff
# TO-DO: Remove ubiquitous genes?
filter_matrix <- function(input_matrix, cutoff){
  output_matrix <- input_matrix[rowSums(input_matrix) >= as.numeric(cutoff), ]
  return(output_matrix)
}

# Merges the results of all individual tools into a CSV file
export_results <- function(bayseq_result, deseq2_result, edgeR_result, limma_result, noiseq_result, sleuth_result, gene_names_list){
  # Write individual results to files
  write.csv(bayseq_result, file = "diffexpr_results_bayseq.csv")
  write.csv(deseq2_result, file = "diffexpr_results_deseq2.csv")
  write.csv(edgeR_result, file = "diffexpr_results_edger.csv")
  write.csv(limma_result, file = "diffexpr_results_limma.csv")
  write.csv(noiseq_result, file = "diffexpr_results_noiseq.csv")
  # Filtering sleuth -> gene subset
  sleuth_result <- subset(sleuth_result, target_id %in% gene_names_list)
  write.csv(sleuth_result, file = "diffexpr_results_sleuth.csv")
  
  # Fetch probabilities/p-values of differential expression
  bayseq_reordered <- bayseq_result[order(match(bayseq_result[[1]], gene_names_list)), ]
  bayseq_FDR <- bayseq_reordered$FDR.DE
  limma_reordered <- limma_result[order(match(rownames(limma_result), gene_names_list)), ]
  limma_pvalues <- limma_reordered$adj.P.Val
  edger_FDR <- unlist(edgeR_result[, "FDR"])[1:length(gene_names_list)]
  
  # Adding missing genes as NA values to NOISeq and sleuth
  for (gene in gene_names_list){
    if ((gene %in% rownames(noiseq_result)) == FALSE){
      new_row_noiseq <- matrix(c(NA, NA, NA, NA, NA, NA), nrow = 1, dimnames = list(gene))
      noiseq_result[nrow(noiseq_result) + 1, ] = as.data.frame(new_row_noiseq)
    }
    if ((gene %in% sleuth_result$target_id) == FALSE){
      new_row_sleuth <- matrix(c(gene, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), nrow = 1)
      sleuth_result[nrow(sleuth_result) + 1, ] = as.data.frame(new_row_sleuth)
    }
  }

  deseq2_pvalues <- deseq2_result[, "padj"]
  noiseq_reordered <- noiseq_result[order(match(rownames(noiseq_result), gene_names_list)), ]
  noiseq_probabilities <- noiseq_reordered$prob
  
  sleuth_reordered <- sleuth_result[order(match(sleuth_result$target_id, gene_names_list)), ]
  sleuth_qvalues <- sleuth_reordered$qval
    
  final_results <- structure(list(baySeq = bayseq_FDR, DESeq2 = deseq2_pvalues, edgeR = edger_FDR, limma = limma_pvalues, NOISeq = noiseq_probabilities, sleuth = sleuth_qvalues), row.names = gene_names_list, class = "data.frame")
  
  # Write summary of all results
  write.csv(final_results, file = "diffexpr_results_all.csv")
  return(final_results)
}

# Transforms the individual predictions into a binary table
# Based on p_value cutoff
probabilities_to_binaries <- function(cutoff_general){
  # Reads the results file and sets the first column as the rownames
  raw_binary_results <- read.csv(file = "diffexpr_results_all.csv", header = TRUE, sep = ",")
  binary_results <- raw_binary_results[,-1]
  rownames(binary_results) <- raw_binary_results[,1]
  
  # Special case for NOISeq since it uses probabilities not p_values
  # In addition, NOISeq was already filtered according to the significance level
  # This means all entrys which are not equal to NA are true DEGs
  noiseq_column <- binary_results$NOISeq
  noiseq_column[is.na(noiseq_column)] <- 0
  noiseq_column[noiseq_column > 0] <- 1
  
  # baySeq, DESeq, edgeR, limma and sleuth results are simply transformed
  # Compare p_values to cutoff
  # This includes the NOISeq column which will be replaced later
  binary_results[is.na(binary_results)] <- 100
  binary_results[binary_results > cutoff_general] <- 100
  binary_results[binary_results <= cutoff_general] <- 0
  binary_results[binary_results == 0] <- 1
  binary_results[binary_results == 100] <- 0
  
  # Overwrite the noiseq results with the precomputed ones
  binary_results$NOISeq <- noiseq_column
  
  write.csv(binary_results, file = "diffexpr_results_all_binary.csv") 
  
  return(binary_results)
}

# Function to call baySeq
run_bayseq <- function(gene_list, gene_counts, raw_replicates_list, total_genes_background){
  DE <- as.numeric(raw_replicates_list == unique(raw_replicates_list)[2])
  groups <- list(NDE = rep(1, length(raw_replicates_list)), DE = DE + 1)
  CD <- new("countData", data = gene_counts, replicates = raw_replicates_list, groups = groups)
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- as.data.frame(gene_list)
  cl <- NULL
  CDP.NBML <- getPriors.NB(CD, samplesize = as.numeric(total_genes_background), estimation = "QL", cl = cl)
  CDPost.NBML <- getLikelihoods(CDP.NBML, pET = "BIC", cl = cl)
  #CDPost.NBML@estProps
  #topCounts(CDPost.NBML, group = 2)
  #NBML.TPs <- getTPs(CDPost.NBML, group = 2, TPs= 1:100)
  # Fetching (all) top hits
  bayseq_de = topCounts(CDPost.NBML, group = 2, number = length(gene_list))
  return(bayseq_de)
}

# Function to call DESeq2
run_deseq2 <- function(list_of_gene_names, sample_counts, sample_conditions){
  # First create a coldata frame and instantiate the DESeqDataSet
  coldata <- data.frame(row.names = colnames(sample_counts), sample_conditions)
  dds <- DESeqDataSetFromMatrix(countData = sample_counts, colData = coldata, design =~ sample_conditions)
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  # Retrieve differential expression results
  res <- results(dds)
  table(res$padj<0.05)
  # Merge with normalized count data and gene symbols
  resdata <- merge(as.matrix(res), as.matrix(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
  new_resdata <- merge(as.matrix(resdata), as.matrix(list_of_gene_names), by = "row.names", sort = FALSE)
  new_resdata <- new_resdata[3:13]
  names(new_resdata)[11] <- "Gene"
  
  # Plots
  # Normalized counts across groups for most significant gene
  plotCounts(dds, gene = which.min(res$padj), intgroup = "sample_conditions")
  # Histogram of adjusted p-values
  hist(res$padj, breaks = 50, col = "grey", main = "Histogram of adjusted p-values", xlab = "p_adjust")
  # Principal component analysis (PCA)
  # This will not work with low-count genes and should be silenced in these cases
  vsd <- vst(dds, blind = FALSE)
  plotPCA(vsd, intgroup = c("sample_conditions"))
  
  return(new_resdata)
}

# Function to call edgeR
# TO-DO: Test different sequencing depth normalization (total count normalization vs. TMM)
# TO-DO: TMM recommended?
run_edger <- function(read_counts, metadata_labels){
  dgList <- DGEList(counts = read_counts, group = metadata_labels)
  # edgeR Normalization
  dgList <- calcNormFactors(dgList)
  dgList <- estimateDisp(dgList)
  # Biological coefficient of variation
  sqrt(dgList$common.dispersion)
  plotBCV(dgList)
  # Test for differential expression between the two classes
  et <- exactTest(dgList)
  edgeR_results <- topTags(et, n = nrow(read_counts), sort.by = "none")
  
  # Plots
  # How many genes are differentially expressed at an FDR of 5%?
  sum(edgeR_results$table$FDR < .05)
  plotSmear(et, de.tags = rownames(edgeR_results)[edgeR_results$table$FDR < .05])
  abline(h = c(-2, 2), col = "blue")
  
  return(edgeR_results)
}

# Function to call voom and limma
run_limma <- function(counts, groups){
  DE <- as.matrix(as.numeric(groups == unique(groups)[2]) + 1)
  DE <- as.data.frame(DE)
  vec <- as.data.frame(c(1))
  DE <- cbind(vec, DE)
  colnames(DE) <- c("Intercept", "Groups")
  DE <- as.matrix(DE)
  dge <- DGEList(counts = counts)
  v <- voom(counts, design = DE, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, DE)
  fit <- eBayes(fit)
  limma_results <- topTable(fit, coef = ncol(DE) - 0, number = length(counts))
  return(limma_results)
}

# Function to call NOIseq
run_noiseq <- function(names_noiseq, counts_noiseq, groups_noiseq, threshold, use_biological_samples){
  list_of_lengths <- read.table(file = "transcript_lengths.csv", sep = ",", header = TRUE)
  internal_threshold = 1 - as.numeric(threshold)
  lengths_DF <- as.data.frame(list_of_lengths$Length, levels(list_of_lengths$Transcript.ID))
  colnames(lengths_DF) <- c("Lengths")
  lengths_DF_new <- lengths_DF[rownames(lengths_DF) %in% names_noiseq, ]
  lengths_DF_new <- as.data.frame(lengths_DF_new, names_noiseq)
  #lengths_DF_turned <- as.data.frame(lengths_DF_new$lengths_DF_new,c("ID", "Length"))
  DE_noiseq <- as.data.frame(as.numeric(groups_noiseq == unique(groups_noiseq)[2]) + 1)
  colnames(DE_noiseq) <- c("Group")
  mydata <- readData(data = counts_noiseq, length = lengths_DF_new, factors = DE_noiseq)
  
  if(use_biological_samples == TRUE){
    # NOISeq biological data mode
    mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "tmm", factor = "Group", lc = 0, r = 50, adj = 1.5, plot = TRUE, a0per = internal_threshold, filter = 1)
  } else {
    # NOISeq non-biological data mode
    mynoiseqbio = noiseq(mydata, k = 0.5, norm = "tmm", factor = "Group", lc = 0, replicates = "biological")
  }

  # TO-DO: Implement degenes instead of filtering genes by myself
  noiseq_results = degenes(mynoiseqbio, q = internal_threshold, M = NULL)
  #noiseq_results <- mynoiseqbio@results[[1]]
  return(noiseq_results)
}

# Function to call sleuth
run_sleuth <- function(metadata_sleuth, benchmarking){
  colnames(metadata_sleuth) <- c("sample", "condition")
  path_list = ""
  if(benchmarking == FALSE){
    kallisto_path = "../mapped/kallisto/"
  } else {
    kallisto_path = "../libraries/miscellaneous/simulation_data/kallisto/"
  }
  for (sample in metadata_sleuth$sample){
    path_raw <- paste(kallisto_path, sample, sep = "")
    path_list <- cbind(path_list, path_raw)
  }
  path_list = path_list[2:(length(metadata_sleuth$sample) + 1)]
  metadata_sleuth_updated <- dplyr::mutate(metadata_sleuth, path = path_list)
  so <- sleuth_prep(metadata_sleuth_updated, num_cores = 1)
  so <- sleuth_fit(so, ~condition, "full")
  so <- sleuth_fit(so, ~1, "reduced")
  so <- sleuth_lrt(so, "reduced", "full")
  sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE)
  return(sleuth_table)
}

# Creates a consensus list of DEGs
# DEGs predicted using a weighted average (majority vote) of individual predictions
# Strict mode requires more than half of the tools to predict a DEG
smart_consensus <- function(b, w, s, t){
  if(s == TRUE){
    consensus_degs <- subset(b, rowWeightedMeans(as.matrix(b), w = w) > t) 
  } else{
    consensus_degs <- subset(b, rowWeightedMeans(as.matrix(b), w = w) >= t)
  }
  return(consensus_degs)
}

# Visualization of the DEGs as a Venn diagram and using UpSetR
# Option to save important diagrams as .png files instead of .pdf
# If 2 tools (e.g. limma and sleuth) both consist of vectors of 0, they will be merged in the upsetR diagram
# TO-DO: Prioritize overlapping 49 genes with all tools
# TO-DO: Then rank rest of genes according to overlaps
# TO-DO: What is the difference between the 88 and the 18 groups of genes detected by single tools?
visualization <- function(binary_table, merge_separate_images){
  number_of_sets = 6
  if (merge_separate_images == FALSE){
    # Venn diagrams
    png(filename = "deg_analysis_venn_diagram_01.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
    v1 <- vennCounts(binary_table[1:3])
    vennDiagram(v1, circle.col = c("blue", "red", "green"))
    dev.off()
    
    png(filename = "deg_analysis_venn_diagram_02.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
    v2 <- vennCounts(binary_table[4:6])
    vennDiagram(v2, circle.col = c("grey", "orange", "cyan"))
    dev.off()
    
    png(filename = "deg_analysis_venn_diagram_03.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
    v3 <- vennCounts(binary_table[1:4])
    vennDiagram(v3, circle.col = c("blue", "red", "green", "grey"))
    dev.off()
    
    # UpSetR images
    png(filename = "deg_analysis_upsetR_diagram.png", width = 20, height = 20, units = "cm", res = 600)
    upset(binary_table, nsets = number_of_sets, mainbar.y.label = "DEG Intersections", sets.x.label = "DEGs Per Tool", order.by = "freq")
    dev.off()

  } else {
    # Venn diagrams
    v1 <- vennCounts(binary_table[1:3])
    vennDiagram(v1, circle.col = c("blue", "red", "green"))
    
    v2 <- vennCounts(binary_table[4:6])
    vennDiagram(v2, circle.col = c("grey", "orange", "cyan")) 
    
    v3 <- vennCounts(binary_table[1:4])
    vennDiagram(v3, circle.col = c("blue", "red", "green", "grey"))
    
    # UpSetR images
    upset(binary_table, nsets = number_of_sets, mainbar.y.label = "DEG Intersections", sets.x.label = "DEGs Per Tool", order.by = "freq")
  }
}

# Main
args <- commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]
argument_3 = args[3]
argument_4 = args[4]
argument_5 = args[5]
argument_6 = args[6]
argument_7 = args[7]
argument_8 = args[8]
argument_9 = args[9]
argument_10 = args[10]
argument_11 = args[11]
argument_12 = args[12]
argument_13 = args[13]
argument_14 = args[14]
argument_15 = args[15]

# Special case if this script is executed manually without any given parameters
# Example: RStudio
if (is.na(argument_1)){
  argument_1 = "Metadata_C1.tsv"
  argument_2 = 5000
  argument_3 = FALSE
  argument_4 = 0.05
  argument_5 = 5
  argument_6 = 0.5
  argument_7 = 1.0
  argument_8 = 0.75
  argument_9 = 0.75
  argument_10 = 2.5
  argument_11 = 0.5
  argument_12 = FALSE
  argument_13 = TRUE
  argument_14 = TRUE
  argument_15 = 0.5
  setwd("../")
}

# Formatting options
# Deactivates e-values
# This is better for the in-script comparisons but makes output less readable
options(scipen = 999)

# Percentage of expected DEGs
# Might need to be adjusted per experiment
benchmark_mode = as.logical(argument_12)
genes_background = argument_2
merge_images = as.logical(argument_3)
noiseq_biological_mode = as.logical(argument_14)
strict_mode = argument_13
threshold_general = argument_4
threshold_expression_count = argument_5
threshold_majority_vote = as.numeric(argument_15)

pdf("deg_analysis_graphs.pdf", paper = "a4")

tool_selection <- c("baySeq", "DESeq2", "edgeR", "limma", "NOISeq", "sleuth")

# Weights of the prediction of individual tools
# Are listed in alphabetical order (important)
weights <- as.numeric(c(argument_6, argument_7, argument_8, argument_9, argument_10, argument_11))

# This is where the script either performs a normal run or a benchmarking run
# Benchmarking is useful for using simulated data from tools like polyester
# It does not perform low-expression-filtering
if (benchmark_mode == FALSE){
  metadata_file = argument_1
  metadata = read.table(file = paste("raw/", metadata_file, sep = ""), sep = "\t", header = FALSE, comment.char = "@")
  metadata_experiments <- metadata$V2
  
  gene_names <- create_gene_list(metadata$V1[1])
  filtered_gene_counts <- create_count_matrix(metadata$V1, gene_names, threshold_expression_count)
  filtered_gene_names <- rownames(filtered_gene_counts)
  
  results_folder = "../../../../deg/"
} else {
  metadata_file = "libraries/miscellaneous/simulation_data/sim_rep_info.txt"
  metadata = read.table(file = metadata_file, sep = "\t", header = TRUE)
  metadata_experiments <- factor(metadata$group)
  
  simulation_table = read.table(file = "libraries/miscellaneous/simulation_data/sim_tx_info.txt", quote = "", sep = "\t", header = TRUE)
  
  load("libraries/miscellaneous/simulation_data/sim_counts_matrix.rda")
  filtered_gene_counts <- counts_matrix
  
  c <- 1
  filtered_gene_names = ""
  
  for (id in rownames(counts_matrix)){
    split_names <- strsplit(id, " +")
    filtered_gene_names[c] <- split_names[[1]][1]
    c <- c + 1
  }

  # Reading the list of true positives (TPs) from the polyester simulation
  simulation_table$transcriptid <- filtered_gene_names
  
  d <- 1
  DE_list = ""
  
  for (gene in simulation_table$transcriptid){
    DE_list[d] <- FALSE
    if (simulation_table[d, 4] == TRUE){
      DE_list[d] <- TRUE
    }
    if (simulation_table[d, 5] == TRUE){
      DE_list[d] <- TRUE
    }
    d <- d + 1
  }
  
  simulation_table_updated <- simulation_table
  simulation_table_updated$foldchange.1 <- DE_list
  simulation_table_updated <- simulation_table_updated[1:2]
  colnames(simulation_table_updated) <- c("IDs", "DE")
  
  gene_names <- filtered_gene_names
  rownames(filtered_gene_counts) <- filtered_gene_names
  results_folder = "deg/"
}

setwd(results_folder)
  
# Performs DEG analyis using the individual tools
# TO-DO: Also possible to try baySeq in edgeR mode?
# TO-DO: Verify normalization for sequencing depth?
time_start <- Sys.time()
results_bayseq = run_bayseq(filtered_gene_names, filtered_gene_counts, metadata_experiments, genes_background)
time_bayseq <- Sys.time()
results_deseq2 = run_deseq2(filtered_gene_names, filtered_gene_counts, metadata_experiments)
time_deseq2 <- Sys.time()
results_edger = run_edger(filtered_gene_counts, metadata_experiments)
time_edger <- Sys.time()
results_limma = run_limma(filtered_gene_counts, metadata_experiments)
time_limma <- Sys.time()
results_noiseq = run_noiseq(filtered_gene_names, filtered_gene_counts, metadata_experiments, threshold_general, noiseq_biological_mode)
time_noiseq <- Sys.time()
results_sleuth = run_sleuth(metadata[1:2], benchmark_mode)
time_sleuth <- Sys.time()

results = export_results(results_bayseq, results_deseq2, results_edger, results_limma, results_noiseq, results_sleuth, filtered_gene_names)
results_binary = probabilities_to_binaries(threshold_general)
results_consensus = smart_consensus(results_binary, weights, strict_mode, threshold_majority_vote)
visualization(results_binary, merge_images)

deg_frame <- c(sum(results_binary$baySeq), sum(results_binary$DESeq2), sum(results_binary$edgeR), sum(results_binary$limma), sum(results_binary$NOISeq), sum(results_binary$sleuth), nrow(results_consensus))
time_frame <- c(difftime(time_bayseq, time_start, units = "secs"), difftime(time_deseq2, time_bayseq, units = "secs"), difftime(time_edger, time_deseq2, units = "secs"), difftime(time_limma, time_edger, units = "secs"), difftime(time_noiseq, time_limma, units = "secs"), difftime(time_sleuth, time_noiseq, units = "secs"), "NA")

# Combine runtime and DEG counts into one summary file for the analysis
summary_frame <- data.frame(Runtimes = time_frame, DEGS = deg_frame)
rownames(summary_frame) <- c(tool_selection, "SCORE")

if (benchmark_mode == TRUE){
  results_binary$DE <- simulation_table_updated$DE
  statistics_frame <- calculate_statistics(results_binary, results_consensus)
  summary_frame$FN <- statistics_frame$FN
  summary_frame$FP <- statistics_frame$FP
  summary_frame$TN <- statistics_frame$TN
  summary_frame$TP <- statistics_frame$TP
  summary_frame$TPR <- statistics_frame$TPR
  summary_frame$TNR <- statistics_frame$TNR
  summary_frame$ACC <- statistics_frame$ACC
  summary_frame$FDR <- statistics_frame$FDR
  summary_frame$FPR <- statistics_frame$FPR
  summary_frame$FNR <- statistics_frame$FNR
  summary_frame$PRE <- statistics_frame$PRE
}

tpm_gene_counts = calculate_tpm(filtered_gene_counts)

write.csv(filtered_gene_counts, file = "filtered_gene_counts_raw.csv")
write.csv(results_consensus, file = "consensus_diffexpr_results.csv")
write.csv(tpm_gene_counts, file = "filtered_gene_counts_tpm.csv")
write.csv(summary_frame, file = "deg_summary.csv")

dev.off()

# Move result graphs file to the other output files
setwd("../")
system("mv deg_analysis_graphs.pdf deg/")