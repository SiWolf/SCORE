# --------------------------------------------
# Title: SCORE.R
# Author: Silver A. Wolf
# Last Modified: Wed, 12.09.2018
# Version: 0.1.5
# --------------------------------------------

#source("https://bioconductor.org/biocLite.R")
#biocLite("baySeq")
#biocLite("DESeq2")
#biocLite("edgeR")

# Imports
library("baySeq")
library("DESeq2")
library("edgeR")

# Functions

# Reads the first counts file to fetch a list of all gene symbols
create_gene_list <- function(sample){
  sample_path <- paste("mapped/bowtie2/featureCounts/", sample, sep = "")
  setwd(sample_path)
  gene_list = read.csv("counts", sep = "", head = T, skip = 1)[,c("Geneid")]
  return(gene_list)
}

# Reads the individual count-files into a matrix
# Applies cutoff to reduce computational time
create_count_matrix <- function(sample_list, gene_list){
  sample_nr = 0
  # Lets receive the rest of the counts...
  for (sample in sample_list){
    sample_nr = sample_nr + 1
    path <- paste("../", sample, sep = "")
    setwd(path)
    counts = read.csv("counts", sep = "", head = T, skip = 1, row.names = 1)
    # Remove first six columns (genesymbol, chr, start, end, strand, length)
    counts <- counts[ ,6:ncol(counts)]
    if (sample_nr == 1){
      count_matrix <- counts
    } else {
      count_matrix <- cbind(count_matrix, counts)  
    }
  }
  # Convert to matrix
  count_matrix <- as.matrix(count_matrix)
  # Rename columns
  colnames(count_matrix) <- paste(sample_list)
  rownames(count_matrix) <- paste(gene_list)
  # unfiltered_count_matrix <- count_matrix
  
  # Cutoff for low-expressed genes
  # Using cpm (edgeR) to calculate counts per gene per million
  # Might add additional filters later on
  # TO-DO: Remove ubiquitous genes?
  cpm_log <- cpm(count_matrix, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  expr_cutoff <- -1
  # Filter out all genes which do not have median(gene_across_samples) > cutoff
  count_matrix <- count_matrix[median_log2_cpm > expr_cutoff, ]
  
  # Initial visualization
  # Histogram
  hist(median_log2_cpm)
  abline(v = expr_cutoff, col = "red", lwd = 3)
  sum(median_log2_cpm > expr_cutoff)
  # Heatmap
  # TO-DO: Heatmap is distorted
  cpm_log_filtered <- cpm(count_matrix, log = TRUE)
  heatmap(cor(cpm_log_filtered))
  # PCA
  pca <- prcomp(t(cpm_log_filtered), scale. = TRUE)
  plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
  text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log_filtered))
  summary(pca)
  
  return(count_matrix)
}

# Merges the results of all individual tools into a CSV file
export_results <- function(bayseq_result, deseq2_result, edgeR_result, gene_names_list){
  # Write results
  results_folder = "../../../../deg/"
  setwd(results_folder)
  write.csv(bayseq_result, file = "bayseq_diffexpr_results_extended.csv")
  write.csv(deseq2_result, file = "deseq2_diffexpr_results_extended.csv")
  write.csv(edgeR_result, file = "edger_diffexpr_results_extended.csv")
  
  bayseq_reordered <- bayseq_result[order(match(bayseq_result$annotation, gene_names_list)), ]
  bayseq_likelihood <- bayseq_reordered$Likelihood
  edger_pvalues <- unlist(edgeR_result[, "PValue"])[1:length(gene_names_list)]
  deseq2_pvalues <- deseq2_result[, "padj"]
  final_results <- structure(list(baySeq = bayseq_likelihood, DESeq2 = deseq2_pvalues, edgeR = edger_pvalues), row.names = gene_names_list, class = "data.frame")
  #final_results <- merge(as.data.frame(edger_pvalues), as.data.frame(deseq2_pvalues))
  
  write.csv(final_results, file = "all_diffexpr_results.csv")
  return(final_results)
}

# Transforms the individual predictions to a binary table
# Uses several cutoff values
probabilities_to_binaries <- function(cutoff_bayseq, cutoff_general, total_number_of_genes){
  # Reads the results file and sets the first column as the rownames
  raw_binary_results <- read.csv(file = "all_diffexpr_results.csv", header = TRUE, sep = ",")
  binary_results <- raw_binary_results[,-1]
  rownames(binary_results) <- raw_binary_results[,1]
  
  # BaySeq uses a dynamic treshold
  # Labels the first 3% (default) of all nonfiltered genes as DE
  # Also possible to use length(bayseq_column) instead of total_number_of_genes
  # Results will be more specific then
  # Possibly will have to change back to this if filters become more strict
  bayseq_column <- binary_results$baySeq
  degs_expected <- round(total_number_of_genes*cutoff_bayseq)
  bayseq_column <- as.data.frame(bayseq_column)
  rownames(bayseq_column) <- raw_binary_results[,1]
  bayseq_column <- bayseq_column[order(-bayseq_column$bayseq_column), , drop = FALSE]
  for (i in 1:degs_expected){
    bayseq_column[1][i, ] <- 1
  }
  bayseq_column[bayseq_column != 1] <- 0
  bayseq_column <- bayseq_column[order(row.names(bayseq_column)), , drop = FALSE]
  
  # DESeq and edgeR results are simply transformed
  # Compare p_values to cutoff
  # This includes the bayseq column but this will be replaced later
  binary_results[is.na(binary_results)] <- 100
  binary_results[binary_results > cutoff_general] <- 100
  binary_results[binary_results <= cutoff_general] <- 0
  binary_results[binary_results == 0] <- 1
  binary_results[binary_results == 100] <- 0
  # Overwrite the current baySeq results with the precomputed ones
  binary_results$baySeq <- bayseq_column$bayseq_column
  return(binary_results)
}

# Function to call baySeq
run_bayseq <- function(gene_list, gene_counts, raw_replicates_list){
  DE <- as.numeric(raw_replicates_list == unique(raw_replicates_list)[2])
  groups <- list(NDE = rep(1, length(metadata$V2)), DE = DE + 1)
  CD <- new("countData", data = gene_counts, replicates = raw_replicates_list, groups = groups)
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- as.data.frame(gene_list)
  cl <- NULL
  CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
  CDPost.NBML <- getLikelihoods(CDP.NBML, pET = "BIC", cl = cl)
  #CDPost.NBML@estProps
  #topCounts(CDPost.NBML, group = 2)
  #NBML.TPs <- getTPs(CDPost.NBML, group = 2, TPs= 1:100)
  # Fetching the top hits
  bayseq_de = topCounts(CDPost.NBML, group = 2, number = length(gene_list))
  return(bayseq_de)
}

# Function to call DESeq2
run_deseq2 <- function(list_of_gene_names, sample_counts, sample_conditions){
  # How to manually assign conditions
  # condition <- factor(c(rep("normal", 2), rep("treated", 2)))
  
  # First create a coldata frame and instantiate the DESeqDataSet
  coldata <- data.frame(row.names = colnames(sample_counts), sample_conditions)
  dds <- DESeqDataSetFromMatrix(countData = sample_counts, colData = coldata, design =~ sample_conditions)
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  # Receive the differential expression results
  res <- results(dds)
  table(res$padj<0.05)
  # Order them by adjusted p-value
  # res <- res[order(res$padj), ]
  # Merge with normalized count data and gene symbols
  resdata <- merge(as.matrix(res), as.matrix(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
  new_resdata <- merge(as.matrix(resdata), as.matrix(list_of_gene_names), by = "row.names", sort = FALSE)
  new_resdata <- new_resdata[3:13]
  names(new_resdata)[11] <- "Gene"
  head(resdata)
  
  # Plots
  # Normalized counts across groups for most significant gene
  plotCounts(dds, gene = which.min(res$padj), intgroup = "sample_conditions")
  # As well as marR
  # plotCounts(dds, gene = "FHHAHKNL_04047", intgroup = "sample_conditions")
  # Histogram of adjusted p-values
  hist(res$padj, breaks = 50, col = "grey", main = "Histogram of adjusted p-values", xlab = "p_adjust")
  # Principal component analysis
  # This will not work with low-count genes and should be silenced in these cases
  vsd <- vst(dds, blind = FALSE)
  plotPCA(vsd, intgroup = c("sample_conditions"))
  
  return(new_resdata)
}

# Function to call edgeR
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
  results_edgeR <- topTags(et, n = nrow(read_counts), sort.by = "none")
  head(results_edgeR$table)
  
  # How many genes are differentially expressed at an FDR of 10%?
  sum(results_edgeR$table$FDR < .1)
  plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
  abline(h = c(-2, 2), col = "blue")
  
  return(results_edgeR)
}

# Export consensus list (uses a majority vote of methods)
smart_consensus <- function(binary_file){
  consensus_degs <- subset(binary_file, rowSums(binary_file)/ncol(binary_file) >= 0.5)
  return(consensus_degs)
}

# Visualization of the DEGs as a venn diagram
visualization_vennDiagram <- function(binary_table){
  v <- vennCounts(binary_table)
  vennDiagram(v, circle.col = c("blue", "red", "green"))
}

# Main
args <- commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]
argument_3 = args[3]

# Special case if this script is run manually using RStudio
if (is.na(argument_1)){
  argument_1 = "Metadata.tsv"
  argument_2 = 0.03
  argument_3 = 0.05
  setwd("../")
}

# Percentage of expected DEGs
# Might need to be adjusted per experiment
threshold_bayseq = argument_2
threshold_general = argument_3

pdf("deg_analysis_graphs.pdf")

metadata = read.table(file = paste("raw/", argument_1, sep = ""), sep = "\t", header = FALSE, comment.char = "@")
gene_names <- create_gene_list(metadata$V1[1])
filtered_gene_counts <- create_count_matrix(metadata$V1, gene_names)
filtered_gene_names <- rownames(filtered_gene_counts)

# TO-DO: Also possible to try baySeq in edgeR mode?
results_bayseq = run_bayseq(filtered_gene_names, filtered_gene_counts, metadata$V2)
results_deseq2 = run_deseq2(filtered_gene_names, filtered_gene_counts, metadata$V2)
results_edger = run_edger(filtered_gene_counts, metadata$V2)

results = export_results(results_bayseq, results_deseq2, results_edger, filtered_gene_names)
results_binary = probabilities_to_binaries(threshold_bayseq, threshold_general, length(gene_names))
results_consensus = smart_consensus(results_binary)
visualization_vennDiagram(results_binary)
write.csv(results_consensus, file = "consensus_diffexpr_results.csv")

dev.off()