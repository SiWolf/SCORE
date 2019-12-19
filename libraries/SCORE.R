# --------------------------------------------
# Title: SCORE.R
# Author: Silver A. Wolf
# Last Modified: Fr, 13.12.2019
# Version: 0.7.4
# --------------------------------------------

# Installers
#install.packages("BiocManager")
#install.packages("devtools")
#install.packages("plyr")
#install.packages("stringr")
#install.packages("UpSetR")
#BiocManager::install("baySeq")
#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("limma")
#BiocManager::install("NOISeq")
#BiocManager::install("rhdf5")
#devtools::install_github("pachterlab/sleuth")

# Imports
library("baySeq")
library("DESeq2")
library("edgeR")
library("ggplot2")
library("limma")
library("NOISeq")
library("plyr")
library("rhdf5")
library("sleuth")
library("stringr")
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
    l = transcript_lengths[transcript_lengths$Transcript_ID == transcript, 2] / 1000
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
  
  # Principal component analysis (PCA) (old)
  # cpm_log_filtered <- cpm(count_matrix, log = TRUE)
  # pca <- prcomp(t(cpm_log_filtered), scale. = TRUE)
  # plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
  # text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log_filtered))
  # summary(pca)
  
  return(count_matrix)
}

# Generate a histogram of the frequency of CPM values
# For all filtered genes
# The count cutoff is visualized as a vertical line
create_cpm_histogram <- function(counts, threshold){
  cpm_log <- cpm(counts, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  hist(median_log2_cpm, col = "grey", main = "Histogram of CPM", xlab = "Median log2CPM", ylab = "Frequency")
  abline(v = threshold, col = "red", lwd = 3)
  sum(median_log2_cpm > threshold)
}

# Removes genes with low counts
# According to a cutoff for low-expressed genes
# Filter out all genes which do not have rowsum >= cutoff
# TO-DO: Remove ubiquitous genes?
filter_matrix <- function(input_matrix, cutoff){
  output_matrix <- input_matrix[rowSums(input_matrix) >= cutoff, ]
  return(output_matrix)
}

# Merges the results of all individual tools into a CSV file
export_results <- function(bayseq_result, deseq2_result, edger_result, limma_result, noiseq_result, sleuth_result, gene_names_list){
  # Write individual results to files
  write.csv(bayseq_result, file = "diffexpr_results_bayseq.csv", row.names = FALSE)
  write.csv(deseq2_result, file = "diffexpr_results_deseq2.csv", row.names = FALSE)
  write.csv(edger_result, file = "diffexpr_results_edger.csv")
  write.csv(limma_result, file = "diffexpr_results_limma.csv")
  write.csv(noiseq_result, file = "diffexpr_results_noiseq.csv")
  
  # Filtering sleuth -> gene subset
  # This will remove genes which did not pass the low expression filter
  sleuth_result <- subset(sleuth_result, target_id %in% gene_names_list)
  write.csv(sleuth_result, file = "diffexpr_results_sleuth.csv", row.names = FALSE)
  
  # Old reordering steps
  #bayseq_result <- bayseq_result[order(match(bayseq_result$gene_list, gene_names_list)), ]
  #limma_result <- limma_result[order(match(rownames(limma_result), gene_names_list)), ]
  
  # Fetch probabilities/p-values of differential expression
  bayseq_FDR <- bayseq_result$FDR.DE
  #deseq2_pvalues <- deseq2_result[, "padj"]
  deseq2_pvalues <- deseq2_result$padj
  #edger_FDR <- unlist(edger_result[, "FDR"])[1:length(gene_names_list)]
  edger_FDR <- edger_result$FDR
  limma_pvalues <- limma_result$adj.P.Val
  
  rownames(sleuth_result) <- sleuth_result$target_id
  
  # Adding missing genes as NA values to NOISeq and sleuth
  for (gene in gene_names_list){
    if ((gene %in% rownames(noiseq_result)) == FALSE){
      new_row_noiseq <- matrix(c(NA, NA, NA, NA, NA, NA), nrow = 1, dimnames = list(gene))
      noiseq_result[nrow(noiseq_result) + 1, ] = as.data.frame(new_row_noiseq)
    }
    if ((gene %in% rownames(sleuth_result)) == FALSE){
      new_row_sleuth <- matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), nrow = 1, dimnames = list(gene))
      sleuth_result[nrow(sleuth_result) + 1, ] = as.data.frame(new_row_sleuth)
    }
  }
  
  sleuth_result <- sleuth_result[2:12]
  
  noiseq_reordered <- noiseq_result[order(match(rownames(noiseq_result), gene_names_list)), ]
  noiseq_probabilities <- noiseq_reordered$prob
  
  sleuth_reordered <- sleuth_result[order(match(rownames(sleuth_result), gene_names_list)), ]
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
  # This means all entrys which are not equal to NA should be true DEGs
  # Filter by 1 - p_value threshold to ensure this
  noiseq_column <- binary_results$NOISeq
  noiseq_column[is.na(noiseq_column)] <- 0
  noiseq_column[noiseq_column < (1 - cutoff_general)] <- 0
  noiseq_column[noiseq_column >= (1 - cutoff_general)] <- 1
  
  # baySeq, DESeq, edgeR, limma and sleuth results are simply transformed
  # By comparing p_values to cutoff
  # Includes the NOISeq column which will be replaced later
  binary_results[is.na(binary_results)] <- 10
  binary_results[binary_results > cutoff_general] <- 10
  binary_results[binary_results <= cutoff_general] <- 20
  binary_results[binary_results == 20] <- 1
  binary_results[binary_results == 10] <- 0
  
  # Overwrite the noiseq results with the precomputed ones
  binary_results$NOISeq <- noiseq_column
  
  write.csv(binary_results, file = "diffexpr_results_all_binary.csv") 
  return(binary_results)
}

# Function to call baySeq
run_bayseq <- function(gene_list, bayseq_gene_counts, raw_replicates_list, total_genes_background){
  DE <- as.numeric(raw_replicates_list == unique(raw_replicates_list)[2])
  groups <- list(NDE = rep(1, length(raw_replicates_list)), DE = DE + 1)
  CD <- new("countData", data = bayseq_gene_counts, replicates = raw_replicates_list, groups = groups)
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- as.data.frame(gene_list)
  cl <- NULL
  CDP.NBML <- getPriors.NB(CD, samplesize = total_genes_background, estimation = "QL", cl = cl)
  CDPost.NBML <- getLikelihoods(CDP.NBML, pET = "BIC", cl = cl)
  #CDPost.NBML@estProps
  #topCounts(CDPost.NBML, group = 2)
  #NBML.TPs <- getTPs(CDPost.NBML, group = 2, TPs= 1:100)
  # Fetching (all) top hits
  bayseq_de = topCounts(CDPost.NBML, group = 2, number = length(gene_list))
  bayseq_de <- bayseq_de[order(bayseq_de$gene_list), ]
  rownames(bayseq_de) <- c()
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
  # Merge with normalized count data and rename first column
  resdata <- merge(as.matrix(res), as.matrix(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
  colnames(resdata)[1] <- "Genes"
  rownames(resdata) <- c()
  
  # Plots
  # Normalized counts across groups for most significant gene
  plotCounts(dds, gene = which.min(res$padj), intgroup = "sample_conditions", xlab = "groups")
  # Histogram of adjusted p-values
  hist(res$padj, breaks = 50, col = "grey", main = "Histogram of adjusted p-values", xlab = "p_adjust", ylab = "frequency")
  # Principal component analysis (PCA)
  # This will not work with low-count genes and should be silenced in these cases
  vsd <- vst(dds, blind = FALSE)
  pca <- plotPCA(vsd, intgroup = c("sample_conditions"), ntop = length(list_of_gene_names), returnData = TRUE)
  pca_gg_plot <- ggplot(pca, aes(x = PC1, y = PC2, colour = sample_conditions)) + geom_point()
  print(pca_gg_plot)
  
  return(resdata)
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
  plotBCV(dgList, main = "Biological Coefficient of Variation vs. Gene Abundance (log2CPM)")
  # Test for differential expression between the two classes
  et <- exactTest(dgList)
  edgeR_results <- topTags(et, n = nrow(read_counts), sort.by = "none")
  
  # Plots
  # How many genes are differentially expressed at an FDR of 5%?
  sum(edgeR_results$table$FDR < .05)
  plotSmear(et, de.tags = rownames(edgeR_results)[edgeR_results$table$FDR < .05], main = "Significant DEGs: logFC vs. logCPM")
  abline(h = c(-2, 2), col = "blue")
  # The following lines allow labeling of individual genes
  # This has been moved to an experimental setting as analyses with many DEGs might cause overlaps between labels
  #extreme_values = edgeR_results[abs(edgeR_results$table$logFC)>2,]
  #n = nrow(extreme_values)
  #for (i in 1:n){
    #text(extreme_values$table$logCPM[i], extreme_values$table$logFC[i], labels = rownames(extreme_values)[i], cex = 0.7, pos = 4)
  #}

  return(edgeR_results$table)
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
  limma_results <- limma_results[order(row.names(limma_results)), ]
  return(limma_results)
}

# Function to call NOIseq
run_noiseq <- function(names_noiseq, counts_noiseq, groups_noiseq, threshold_noiseq, use_biological_samples){
  # Import and format list of lengths
  # Lengths corresponding to fragments of the same gene (identifier) are summed
  list_of_lengths <- read.table(file = "transcript_lengths.csv", sep = ",", header = TRUE)
  list_of_lengths <- ddply(list_of_lengths, "Transcript_ID", numcolwise(sum))
  list_of_lengths <- list_of_lengths[order(list_of_lengths$Transcript_ID), ]
  rownames(list_of_lengths) <- c()
  
  # Convert list to data frame
  lengths_DF <- as.data.frame(list_of_lengths$Length, levels(list_of_lengths$Transcript_ID))
  colnames(lengths_DF) <- c("Lengths")
  
  # Filter DF according to pre-filtered gene list
  lengths_DF_new <- lengths_DF[rownames(lengths_DF) %in% names_noiseq, ]
  lengths_DF_new <- as.data.frame(lengths_DF_new, names_noiseq)
  #lengths_DF_turned <- as.data.frame(lengths_DF_new$lengths_DF_new,c("ID", "Length"))
  
  # DEG Prediction
  DE_noiseq <- as.data.frame(as.numeric(groups_noiseq == unique(groups_noiseq)[2]) + 1)
  colnames(DE_noiseq) <- c("Group")
  internal_threshold = 1 - threshold_noiseq
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
  noiseq_results <- noiseq_results[order(rownames(noiseq_results)), ]
  #noiseq_results <- mynoiseqbio@results[[1]]
  return(noiseq_results)
}

# Function to call sleuth
run_sleuth <- function(metadata_sleuth, benchmarking, identifier_sleuth){
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
  
  # Renaming long names to make table easier to read (old)
  # Increases compatibility with merging later
  #first_split <- str_split_fixed(sleuth_table$target_id, paste(identifier_sleuth, "=", sep = ""), 2)
  #first_split <- as.data.frame(first_split)
  #second_split <- str_split_fixed(first_split$V2, "]", 2)
  #second_split <- as.data.frame(second_split)
  #sleuth_table$target_id <- second_split$V1
  
  # Filters all NA-value-genes from table
  # Checks for duplicate IDs and will only keep the first hit
  sleuth_table <- sleuth_table[!is.na(sleuth_table$pval), ]
  sleuth_table <- sleuth_table[!duplicated(sleuth_table$target_id), ]
  sleuth_table <- sleuth_table[order(sleuth_table$target_id), ]
  rownames(sleuth_table) <- c()
  
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
argument_16 = args[16]

# Special case if this script is executed manually without any given parameters
# Example: RStudio
if (is.na(argument_1)){
  argument_1 = "raw/Metadata.tsv"
  argument_2 = 5000
  argument_3 = TRUE
  argument_4 = 0.05
  argument_5 = 10
  argument_6 = 1.0
  argument_7 = 1.0
  argument_8 = 1.0
  argument_9 = 1.0
  argument_10 = 1.0
  argument_11 = 1.0
  argument_12 = FALSE
  argument_13 = TRUE
  argument_14 = TRUE
  argument_15 = 0.5
  argument_16 = "locus_tag"
  setwd("../")
}

# Formatting options
# Deactivates e-values
# This is better for the in-script comparisons but makes output less readable
options(scipen = 999)

# Percentage of expected DEGs
# Might need to be adjusted per experiment
benchmark_mode = as.logical(argument_12)
genes_background = as.numeric(argument_2)
identifier = argument_16
merge_images = as.logical(argument_3)
noiseq_biological_mode = as.logical(argument_14)
strict_mode = as.logical(argument_13)
threshold_general = as.numeric(argument_4)
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
  metadata = read.table(file = metadata_file, sep = "\t", header = FALSE, comment.char = "@")
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

# Ensure lexicographical ordering of genes
filtered_gene_counts <- filtered_gene_counts[order(row.names(filtered_gene_counts)), ]
filtered_gene_names <- sort(filtered_gene_names)
gene_names <- sort(gene_names)

# Initial visualization
# Histogram
create_cpm_histogram(filtered_gene_counts, threshold_expression_count)

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
results_sleuth = run_sleuth(metadata[1:2], benchmark_mode, identifier)
time_sleuth <- Sys.time()

results = export_results(results_bayseq, results_deseq2, results_edger, results_limma, results_noiseq, results_sleuth, filtered_gene_names)
results_binary = probabilities_to_binaries(threshold_general)
results_consensus = smart_consensus(results_binary, weights, strict_mode, threshold_majority_vote)
visualization(results_binary, merge_images)

deg_frame <- c(sum(results_binary$baySeq), sum(results_binary$DESeq2), sum(results_binary$edgeR), sum(results_binary$limma), sum(results_binary$NOISeq), sum(results_binary$sleuth), nrow(results_consensus))
time_frame <- c(difftime(time_bayseq, time_start, units = "secs"), difftime(time_deseq2, time_bayseq, units = "secs"), difftime(time_edger, time_deseq2, units = "secs"), difftime(time_limma, time_edger, units = "secs"), difftime(time_noiseq, time_limma, units = "secs"), difftime(time_sleuth, time_noiseq, units = "secs"), " ")

# Combine runtime and DEG counts into one summary file for the analysis
summary_frame <- data.frame(Runtimes = time_frame, DEGs = deg_frame)
rownames(summary_frame) <- c(tool_selection, "SCORE")

# Output additional stats such as the total amount of genes and filtered genes
percentage_degs = round((nrow(results_consensus)/length(gene_names))*100, digits = 2)
percentage_genes = round((length(filtered_gene_names)/length(gene_names))*100, digits = 2)
additional_stats_c1 <- c("Total", "Filtered", "DEGs-SCORE")
additional_stats_c2 <- c(length(gene_names), paste(length(filtered_gene_names), " (", percentage_genes, "%)", sep = ""), paste(nrow(results_consensus), " (", percentage_degs, "%)", sep = ""))
additional_stats_df <- data.frame(Genes = additional_stats_c1, Amount = additional_stats_c2)

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

write.csv(additional_stats_df, file = "deg_stats.csv", row.names = FALSE)
write.csv(filtered_gene_counts, file = "filtered_gene_counts_raw.csv")
write.csv(results_consensus, file = "consensus_diffexpr_results.csv")
write.csv(summary_frame, file = "deg_summary.csv")
write.csv(tpm_gene_counts, file = "filtered_gene_counts_tpm.csv")

dev.off()

# Move result graphs file to the other output files
setwd("../")
system("mv deg_analysis_graphs.pdf deg/")