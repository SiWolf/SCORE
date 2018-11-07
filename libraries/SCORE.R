# --------------------------------------------
# Title: SCORE.R
# Author: Silver A. Wolf
# Last Modified: Wed, 07.11.2018
# Version: 0.3.8
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

# Reads the first counts-file in order to fetch a list of all gene symbols
create_gene_list <- function(sample){
  sample_path <- paste("mapped/bowtie2/featureCounts/", sample, sep = "")
  setwd(sample_path)
  gene_list = read.csv(paste("counts", sample, sep = "_"), sep = "", head = T, skip = 1)[,c("Geneid")]
  return(gene_list)
}

# Reads the individual count-files into a matrix
# Applies a low-expression cutoff to reduce computational time
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
  
  # Cutoff for low-expressed genes
  # Removed genes with low counts
  # Filter out all genes which do not have rowsum >= cutoff
  # TO-DO: Remove ubiquitous genes?
  count_matrix <- count_matrix[rowSums(count_matrix) >= as.numeric(low_expression_cutoff), ]
  
  # Initial visualization
  # Histogram
  cpm_log <- cpm(count_matrix, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  hist(median_log2_cpm)
  abline(v = low_expression_cutoff, col = "red", lwd = 3)
  sum(median_log2_cpm > low_expression_cutoff)
  # Heatmap
  # TO-DO: Heatmap is distorted!
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
export_results <- function(bayseq_result, deseq2_result, edgeR_result, limma_result, noiseq_result, sleuth_result, gene_names_list){
  # Write results
  write.csv(bayseq_result, file = "diffexpr_results_bayseq.csv")
  write.csv(deseq2_result, file = "diffexpr_results_deseq2.csv")
  write.csv(edgeR_result, file = "diffexpr_results_edger.csv")
  write.csv(limma_result, file = "diffexpr_results_limma.csv")
  write.csv(noiseq_result, file = "diffexpr_results_noiseq.csv")
  sleuth_result <- subset(sleuth_result, target_id %in% gene_names_list)
  write.csv(sleuth_result, file = "diffexpr_results_sleuth.csv")
  
  bayseq_reordered <- bayseq_result[order(match(bayseq_result$annotation, gene_names_list)), ]
  bayseq_likelihood <- bayseq_reordered$Likelihood
  limma_reordered <- limma_result[order(match(rownames(limma_result), gene_names_list)), ]
  limma_pvalues <- limma_reordered$adj.P.Val
  edger_pvalues <- unlist(edgeR_result[, "PValue"])[1:length(gene_names_list)]
  deseq2_pvalues <- deseq2_result[, "padj"]
  noiseq_probabilities <- noiseq_result$prob
  # TO-DO: Replace pval with qval
  sleuth_reordered <- sleuth_result[order(match(sleuth_result$target_id, gene_names_list)), ]
  sleuth_pvalues <- sleuth_reordered$pval
  final_results <- structure(list(baySeq = bayseq_likelihood, DESeq2 = deseq2_pvalues, edgeR = edger_pvalues, limma = limma_pvalues, NOISeq = noiseq_probabilities, sleuth = sleuth_pvalues), row.names = gene_names_list, class = "data.frame")
  #final_results <- merge(as.data.frame(edger_pvalues), as.data.frame(deseq2_pvalues))
  
  write.csv(final_results, file = "diffexpr_results_all.csv")
  return(final_results)
}

# Transforms the individual predictions to a binary table
# Uses several set cutoff values
probabilities_to_binaries <- function(cutoff_bayseq, cutoff_general, total_number_of_genes){
  # Reads the results file and sets the first column as the rownames
  raw_binary_results <- read.csv(file = "diffexpr_results_all.csv", header = TRUE, sep = ",")
  binary_results <- raw_binary_results[,-1]
  rownames(binary_results) <- raw_binary_results[,1]
  
  # BaySeq uses a dynamic treshold
  # Labels the first 3% (default) of all nonfiltered genes as DE
  # Also possible to use length(bayseq_column) instead of total_number_of_genes
  # Results will be more specific then
  # Possibly will have to change back to this if filters become more strict
  bayseq_column <- binary_results$baySeq
  degs_expected <- round(total_number_of_genes*as.numeric(cutoff_bayseq))
  bayseq_column <- as.data.frame(bayseq_column)
  rownames(bayseq_column) <- raw_binary_results[,1]
  bayseq_column <- bayseq_column[order(-bayseq_column$bayseq_column), , drop = FALSE]
  for (i in 1:degs_expected){
    bayseq_column[1][i, ] <- 1
  }
  bayseq_column[bayseq_column != 1] <- 0
  bayseq_column <- bayseq_column[order(row.names(bayseq_column)), , drop = FALSE]
  
  # Special case for NOISeq since it uses probabilities not p_values
  noiseq_column <- binary_results$NOISeq
  noiseq_column[noiseq_column >= 1 - as.numeric(cutoff_general)] <- 1
  noiseq_column[noiseq_column < 1 - as.numeric(cutoff_general)] <- 0
  noiseq_column[is.na(noiseq_column)] <- 0
  
  # DESeq, edgeR, limma and sleuth results are simply transformed
  # Compare p_values to cutoff
  # This includes the bayseq and noiseq columns which will be replaced later
  binary_results[is.na(binary_results)] <- 100
  binary_results[binary_results > cutoff_general] <- 100
  binary_results[binary_results <= cutoff_general] <- 0
  binary_results[binary_results == 0] <- 1
  binary_results[binary_results == 100] <- 0
  # Overwrite the current bayseq and noiseq results with the precomputed ones
  binary_results$baySeq <- bayseq_column$bayseq_column
  binary_results$NOISeq <- noiseq_column
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
  # How to manually assign conditions:
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
  head(edgeR_results$table)
  
  # How many genes are differentially expressed at an FDR of 10%?
  sum(edgeR_results$table$FDR < .1)
  plotSmear(et, de.tags = rownames(edgeR_results)[edgeR_results$table$FDR < .1])
  abline(h = c(-2, 2), col = "blue")
  
  return(edgeR_results)
}

# Function to call voom and limma
run_limma <- function(counts, groups){
  #groups <- metadata$V2
  #counts <- filtered_gene_counts
  DE <- as.matrix(as.numeric(groups == unique(groups)[2]) + 1)
  DE <- as.data.frame(DE)
  vec <- as.data.frame(c(1))
  DE <- cbind(vec, DE)
  colnames(DE) <- c("Intercept", "Groups")
  #design <- model.matrix(~ Groups, DE)
  DE <- as.matrix(DE)
  
  #DE <- model.matrix(~0+factor(c(1,1,1,2,2,2)))
  
  dge <- DGEList(counts = counts)
  
  #dge <- calcNormFactors(dge)
  
  #logCPM <- cpm(dge, log=TRUE, prior.count=3)
  #fit <- lmFit(logCPM, DE)
  #fit <- eBayes(fit, trend=TRUE)
  #s <- topTable(fit, coef=ncol(DE) -1, number=length(counts))
  
  #s <- topTable(fit, adjust.method="fdr", number=length(counts))
  
  v <- voom(counts, design = DE, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, DE)
  fit <- eBayes(fit)
  #limma_results <- topTable(fit, coef = ncol(DE), number = length(counts))
  limma_results <- topTable(fit, coef = ncol(DE) - 1, number = length(counts))
  return(limma_results)
}

# Function to call NOIseq
run_noiseq <- function(counts_noiseq, groups_noiseq){
  list_of_lengths <- read.table(file = "transcript_lengths.csv", sep = ",", header = TRUE)
  
  lengths_DF <- as.data.frame(list_of_lengths$Length, levels(list_of_lengths$Transcript.ID))
  colnames(lengths_DF) <- c("Lengths")
  
  lengths_DF_new <- lengths_DF[rownames(lengths_DF) %in% filtered_gene_names, ]
  lengths_DF_new <- as.data.frame(lengths_DF_new, filtered_gene_names)
  
  #lengths_DF_turned <- as.data.frame(lengths_DF_new$lengths_DF_new,c("ID", "Length"))
  
  DE_noiseq <- as.data.frame(as.numeric(groups_noiseq == unique(groups_noiseq)[2]) + 1)
  colnames(DE_noiseq) <- c("Group")
  mydata <- readData(data = counts_noiseq, length = lengths_DF_new, factors = DE_noiseq)

  mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "rpkm", factor = "Group", lc = 0, r = 20, adj = 1.5, plot = TRUE, a0per = 0.9, random.seed = 12345, filter = 1)
  
  noiseq_results <- mynoiseqbio@results[[1]]
  
  return(noiseq_results)
}

# Function to call sleuth
run_sleuth <- function(metadata_sleuth){
  colnames(metadata_sleuth) <- c("sample", "condition")
  path_list = ""
  for (sample in metadata_sleuth$sample){
    path_raw <- paste("../mapped/kallisto/", sample, sep = "")
    path_list <- cbind(path_list, path_raw)
  }
  path_list = path_list[2:(length(metadata$V1) + 1)]
  metadata_sleuth_updated <- dplyr::mutate(metadata_sleuth, path = path_list)
  so <- sleuth_prep(metadata_sleuth_updated, num_cores = 1)
  so <- sleuth_fit(so, ~condition, "full")
  so <- sleuth_fit(so, ~1, "reduced")
  so <- sleuth_lrt(so, "reduced", "full")
  sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE)
  return(sleuth_table)
}

# Creates a consensus list of DEGs
# DEGs predicted using a weighted average (majority vote) of individual predicitons
smart_consensus <- function(binary_file, w){
  consensus_degs <- subset(binary_file, rowWeightedMeans(as.matrix(binary_file), w = w) >= 0.5)
  return(consensus_degs)
}

# Visualization of the DEGs as a venn diagram and using upsetr
# TO-DO: Prioritize overlapping 49 genes with all tools
# TO-DO: Then rank rest of genes according to overlaps
# TO-DO: What is the difference between the 88 and the 18 groups of genes detected by single tools?
visualization <- function(binary_table){
  # Venn diagrams
  v1 <- vennCounts(binary_table[1:3])
  vennDiagram(v1, circle.col = c("blue", "red", "green"))
  v2 <- vennCounts(binary_table[4:6])
  vennDiagram(v2, circle.col = c("grey", "orange", "cyan"))
  v3 <- vennCounts(binary_table[1:4])
  vennDiagram(v3, circle.col = c("blue", "red", "green", "grey"))
  # UpsetR images
  upset(binary_table, nsets = 6, mainbar.y.label = "DEG Intersections", sets.x.label = "DEGs Per Tool", order.by = "freq")
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

# Special case if this script is run manually using RStudio
if (is.na(argument_1)){
  argument_1 = "Metadata_C1.tsv"
  argument_2 = 0.03
  argument_3 = 0.05
  argument_4 = 5
  argument_5 = 1.0
  argument_6 = 1.0
  argument_7 = 1.0
  argument_8 = 1.0
  argument_9 = 1.0
  argument_10 = 1.0
  setwd("../")
}

# Percentage of expected DEGs
# Might need to be adjusted per experiment
threshold_bayseq = argument_2
threshold_general = argument_3
threshold_expression_count = argument_4
# Weights of the prediction of individual tools
# Are listed in alphabetical order (highly important)
weights <- as.numeric(c(argument_5, argument_6, argument_7, argument_8, argument_9, argument_10))

pdf("deg_analysis_graphs.pdf", paper = "a4")

metadata = read.table(file = paste("raw/", argument_1, sep = ""), sep = "\t", header = FALSE, comment.char = "@")
gene_names <- create_gene_list(metadata$V1[1])
filtered_gene_counts <- create_count_matrix(metadata$V1, gene_names, threshold_expression_count)
filtered_gene_names <- rownames(filtered_gene_counts)

results_folder = "../../../../deg/"
setwd(results_folder)

# TO-DO: Also possible to try baySeq in edgeR mode?
# TO-DO: Normalization for sequencing depth?
time_start <- Sys.time()
results_bayseq = run_bayseq(filtered_gene_names, filtered_gene_counts, metadata$V2)
time_bayseq <- Sys.time()
results_deseq2 = run_deseq2(filtered_gene_names, filtered_gene_counts, metadata$V2)
time_deseq2 <- Sys.time()
results_edger = run_edger(filtered_gene_counts, metadata$V2)
time_edger <- Sys.time()
results_limma = run_limma(filtered_gene_counts, metadata$V2)
time_limma <- Sys.time()
results_noiseq = run_noiseq(filtered_gene_counts, metadata$V2)
time_noiseq <- Sys.time()
results_sleuth = run_sleuth(metadata)
time_sleuth <- Sys.time()

time_frame <- c(difftime(time_bayseq, time_start, units = "secs"), difftime(time_deseq2, time_bayseq, units = "secs"), difftime(time_edger, time_deseq2, units = "secs"), difftime(time_limma, time_edger, units = "secs"), difftime(time_noiseq, time_limma, units = "secs"), difftime(time_sleuth, time_noiseq, units = "secs"))

results = export_results(results_bayseq, results_deseq2, results_edger, results_limma, results_noiseq, results_sleuth, filtered_gene_names)
results_binary = probabilities_to_binaries(threshold_bayseq, threshold_general, length(gene_names))
results_consensus = smart_consensus(results_binary, weights)
visualization(results_binary)
write.csv(results_consensus, file = "consensus_diffexpr_results.csv")
write.csv(filtered_gene_counts, file = "filtered_gene_counts.csv")

# Combine runtime and DEG counts into one summary file for the analysis
deg_frame <- c(sum(results_binary$baySeq), sum(results_binary$DESeq2), sum(results_binary$edgeR), sum(results_binary$limma), sum(results_binary$NOISeq), sum(results_binary$sleuth))
summary_frame <- data.frame(DEGS = deg_frame, Runtimes = time_frame)
rownames(summary_frame) <- c("baySeq", "DESeq2", "edgeR", "limma", "NOISeq", "sleuth")
write.csv(summary_frame, file = "deg_summary.csv")

# TO-DO: Add raw counts to final output file
# Possibly in Python file?

dev.off()

# TO-DO: Which (open) licence will be used? GPL?