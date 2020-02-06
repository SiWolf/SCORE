# -----------------------------
# Title: visualize_kegg.R
# Author: Silver A. Wolf
# Last Modified: Thur, 06.02.2020
# Version: 0.0.4
# -----------------------------

# Installers
#install.packages("BiocManager")
#install.packages("ggplot2")
#BiocManager::install("pathview")

# Imports
library("ggplot2")
library("pathview")

# Functions

# Better GeneSCF images
visualize_genescf <- function(KEGG_data, GO_data){
  for (i in 1:2) {
    if (i == 1) {
      name = "KEGG"
      tmp_data <- KEGG_data[order(KEGG_data$P.value), ]
    }
    else {
      name = "GO"
      tmp_data <- GO_data[order(GO_data$P.value), ]
    }
    #size <- length(tmp_data[, 1])
    tmp_data[, "Rank_in_order"] <- c(seq(1:length(tmp_data[, 1])))
    tmp_data <- tmp_data[, cbind("Rank_in_order", "Process.name", "percentage.", "P.value")]
    colnames(tmp_data) <- c("Rank", "Process", "Genes", "Pval")
    tmp_data <- tmp_data[1:20, ]
    newdata <- cbind(gsub( "~.*$", "", tmp_data[, "Process"]), tmp_data[, "Rank"])
    colnames(newdata) <- c("IDs", "Rank")
    png(paste("../../pathway_analysis_", name, "/enrichment_plot.png", sep = ""), width = 1200, height = 800)
    # In case we find a zero p-value, we replace it by its larger neighbor and shrink it by an arbitrary value for better visualization
    tmp_data$Pval[tmp_data$Pval == 0] <- min(tmp_data$Pval[tmp_data$Pval > 0]) * 0.05
    fill_data = paste0("[", sprintf("%02d", as.numeric(tmp_data$Rank)), "]\t", tmp_data[, "Process"])
    plot <- ggplot(tmp_data, aes(x = Rank, y = -log10(Pval), size = Genes, label = newdata[, "IDs"], fill = fill_data, guide = FALSE)) +
      geom_point(colour = "#2E2E2E", shape = 21) +
      scale_size_area(max_size = 5) +
      labs(fill = "[Rank] Process (Top 20)", title = paste("Top 20 ", name, " Hits", sep = "")) +
      scale_x_continuous(name = "Rank", limits = c(0, 22)) +
      scale_y_continuous(name = "-log10(Pval)", limits = c(1.1, max(-log10(tmp_data$Pval)) + 1)) +
      #scale_size_continuous(range=c(1, 15)) +
      geom_text(size = 5, color = "#2E2E2E", hjust = -0.1, vjust = 0, angle = 45) +
      geom_hline(yintercept = 1.3) +
      theme(text = element_text(size = 14), plot.title = element_text(hjust = 0.5))
    print(plot)
    dev.off()
  }
}

# Visualize KEGG pathways using pathview
visualize_pathview <- function(d1, d2, d3, d4){
  DEGs <- d1[d1[5] != 0, ]
  
  c <- 1
  fc <- c()
  
  for (gene in d2$V2) {
    #gene_id <- d2[c, 1]
    hits <- DEGs[DEGs$gene.name == gene, ]
    pos = 1
    if (nrow(hits) > 1) {
      pos <- which(hits$log2FC == max(abs(hits$log2FC)) | hits$log2FC == -max(abs(hits$log2FC)))
    }
    fc[c] <- hits$log2FC[pos]
    c = c + 1
  }
  
  #df <- data.frame(IDs = d2[, 1], logFC = fc)
  names(fc) <- d2[, 1]
  
  pathways <- sub("~.*", "", d3$Process.name)
  
  # Go through pathways in a loop (try) since this is more stable than parsing the full pathway list
  for (pathway in pathways) {
    pv <- try(pathview(gene.data = fc, species = d4, pathway.id = pathway, min.nnodes = 0, gene.idtype = "KEGG"))
  }
}

# Main
args <- commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]

if (is.na(argument_1)) {
  argument_1 = "ecocyc"
  argument_2 = "eco"
}

go_species = argument_1
kegg_species = argument_2

dir.create("deg/pathway_analysis_KEGG/pathview/")
setwd("deg/pathway_analysis_KEGG/pathview/")

data_1 <- read.delim(file = "../../summary.tsv", header = TRUE)
data_2 <- read.delim(file = "../deg_gene_symbols_user_mapped.list", header = FALSE)
data_3 <- read.delim(file = paste("../deg_gene_symbols_KEGG_", kegg_species, "_functional_classification.tsv", sep = ""), header = TRUE)
data_4 <- read.delim(file = paste("../../pathway_analysis_GO/deg_gene_symbols_GO_all_", go_species, "_functional_classification.tsv", sep = ""), header = TRUE)

visualize_genescf(data_3, data_4)
visualize_pathview(data_1, data_2, data_3, kegg_species)