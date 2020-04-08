# --------------------------------------------
# Title: generate_additional_images.R
# Author: Silver A. Wolf
# Last Modified: Wed, 08.04.2020
# Version: 0.3.2
# --------------------------------------------

# This script is used to generate additional images for publications, etc.
# It will need to be heavily adapted for any new data and is used to test new plots
# Finished plots are moved to the main SCORE scripts during development

# Libraries

library("dplyr")
library("ggplot2")
library("limma")
library("pheatmap")
library("reshape2")
library("stringdist")
library("stringr")
library("viridis")

# Variables

export = FALSE

# Preprocessing

setwd("../../deg/V1/")

final_summary <- read.csv(file = "summary.tsv", header = TRUE, sep = "\t", quote = "")
final_summary_V1_extended <- read.csv(file = "summary_V1.tsv", header = TRUE, sep = "\t", quote = "", dec = ",")

for (i in 1:nrow(final_summary)){
  if (as.character(final_summary$TIGRFAM.main.role[i]) == ""){
    final_summary$TIGRFAM.main.role[i] <- "Unknown function"
  }
}

degs = final_summary[final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. != 0, ]
tpm_values <- read.csv(file = "filtered_gene_counts_tpm.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)
tpm_values$TIGRFAM <- final_summary$TIGRFAM.main.role

if (export == TRUE){
  pdf("deg_analysis_graphs_extended.pdf", paper = "a4")
}

# V1

# Heatmaps

# Apply log2FC cutoff to narrow down list of targets
high_FC <- degs[abs(degs$log2FC) > 2, ]
tpm_degs <- tpm_values[tpm_values$ID %in% high_FC$ID, ]
tpm_degs_ordered <- tpm_degs[order(tpm_degs$TIGRFAM, tpm_degs$gene.name), ]
tpm_degs_ordered$gene.name <- replace(tpm_degs_ordered$gene.name, tpm_degs_ordered$gene.name == "", "hypothetical_gene")
rownames(tpm_degs_ordered) <- tpm_degs_ordered$gene.name

groups = data.frame(Groups = c(rep("Control_Group", 3), rep("Zn_Group", 3)))
rownames(groups) <- c("X6122_C1", "X6122_C5", "X6122_C6", "X6122_Z1", "X6122_Z5", "X6122_Z6")

anno_colours = list(Groups = c(Control_Group = "#96ca00", Zn_Group = "#ff82cb"))

g1 <- pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
        annotation_col = groups, annotation_colors = anno_colours, fontsize_row = 9, fontsize_col = 9,
        cellheight = 10, cellwidth = 20,
        border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
        main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
        scale = "none", breaks = seq(0, 10000, 100))

g2 <- pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
         annotation_col = groups, annotation_colors = anno_colours, fontsize_row = 9, fontsize_col = 9,
         cellheight = 10, cellwidth = 20,
         border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
         main = "Heatmap of the gene expression (TPM)\nbetween sample groups (scaled)\n",
         scale = "row")

# Boxplots

g3 <- ggplot(data = final_summary, aes(x = TIGRFAM.main.role, y = log2FC)) + 
  geom_boxplot(fill = "lightgrey") +
  stat_boxplot(geom = "errorbar") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Boxplot of the log2FC of all genes within each TIGRFAM Main Role (Baseline)") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC")

g4 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = log2FC)) + 
  geom_boxplot(fill = "lightgrey") +
  stat_boxplot(geom = "errorbar") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Boxplot of the log2FC of DEGs within each TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC")

# Histograms

g5 <- ggplot(data = final_summary, aes(final_summary$log2FC)) +
  geom_histogram(breaks = seq(-2.5, 2.5, 0.05), col = "black", fill = "darkgrey", alpha = 0.2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Histogram of the log2FC (all genes)") +
  xlab("log2FC") +
  ylab("counts (genes)")

g6 <- ggplot(data = degs, aes(degs$log2FC)) +
  geom_histogram(breaks = seq(-2.5, 2.5, 0.05), col = "black", fill = "darkgrey", alpha = 0.2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Histogram of the log2FC (DEGs)") +
  xlab("log2FC") +
  ylab("counts (DEGs)")

# Barcharts

BF1 <- data.frame(degs$ID, degs$Mean.TPM..Wildtype., "G1", degs$TIGRFAM.main.role)
BF2 <- data.frame(degs$ID, degs$Mean.TPM..ZnGroup., "G2", degs$TIGRFAM.main.role)
BF3 <- data.frame(TIGRFAM = degs$TIGRFAM.main.role, WT = degs$Mean.TPM..Wildtype., Zn = degs$Mean.TPM..ZnGroup.)
BF3_hot <- melt(BF3)
BFtranspose <- cbind(t(BF1), t(BF2))
BFtranspose <- as.data.frame(t(BFtranspose))
colnames(BFtranspose) <- c("ID", "TPM", "Group", "TIGRFAM_Main_Role")

g7 <- ggplot(data = final_summary, aes(x = factor(TIGRFAM.main.role))) +
  geom_bar(stat = "count") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the counts (all genes) per TIGRFAM Main Role (Baseline)") +
  xlab("TIGRFAM Main Roles") +
  ylab("counts (genes)")

g8 <- ggplot(data = degs, aes(x = factor(TIGRFAM.main.role))) +
  geom_bar(stat = "count") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the counts (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("counts (DEGs)")

g9 <- ggplot(data = degs, aes(fill = degs$TIGRFAM.main.role, y = degs$Mean.TPM..Wildtype., x = degs$TIGRFAM.main.role)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the max. TPM (G1) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("max. TPM (G1)") +
  labs(fill = "TIGRFAM Main Roles")

g10 <- ggplot(data = degs, aes(fill = degs$TIGRFAM.main.role, y = degs$Mean.TPM..ZnGroup., x = degs$TIGRFAM.main.role)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the max. TPM (G2) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("max. TPM (G2)") +
  labs(fill = "TIGRFAM Main Roles")

g11 <- ggplot(data = BF3_hot, aes(x = TIGRFAM, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the max. TPM (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("max. TPM (DEGs)") +
  labs(fill = "Groups")

terms_list <- sort(unique(degs$TIGRFAM.main.role))
g1_list <- c()
g2_list <- c()

for (term in terms_list){
  current_list <- degs[degs$TIGRFAM.main.role == term, ]
  sum_g1 <- sum(current_list$Mean.TPM..Wildtype.)
  sum_g2 <- sum(current_list$Mean.TPM..ZnGroup.)
  g1_list <- c(g1_list, sum_g1/length(current_list$Mean.TPM..Wildtype.))
  g2_list <- c(g2_list, sum_g2/length(current_list$Mean.TPM..ZnGroup.))
}

deg_data_frame <- data.frame(TIGRFAM = rep(terms_list, 2), tpm_sum = c(g1_list, g2_list), group = c(rep("G1", length(unique(degs$TIGRFAM.main.role))), rep("G2", length(unique(degs$TIGRFAM.main.role)))))

g12 <- ggplot(deg_data_frame, aes(fill = group, y = tpm_sum, x = TIGRFAM)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the total TPM (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("Normalized TPM") +
  labs(fill = "Groups")

g13 <- ggplot(final_summary, aes(y = log2FC, x = TIGRFAM.main.role)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the log2FC (all genes) per TIGRFAM Main Role (Baseline)") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC") +
  labs(fill = "Groups")

g14 <- ggplot(degs, aes(y = log2FC, x = TIGRFAM.main.role)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the log2FC (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC") +
  labs(fill = "Groups")

# Old Plots

g15 <- ggplot(data = BFtranspose, aes(x = ID, y = TPM, colour = Group)) +
  geom_bar(stat = "identity") +
  theme_void()

g16 <- ggplot(data = BFtranspose, aes(x = TIGRFAM_Main_Role, y = TPM, colour = Group)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the total TPM per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("Total TPM")

g17 <- ggplot(data = BFtranspose, aes(x = factor(ID), y = TPM, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM per Gene") +
  xlab("Genes") +
  ylab("Total TPM")

g18 <- ggplot(data = BFtranspose, aes(x = factor(TIGRFAM_Main_Role), y = TPM, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the total TPM per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("Total TPM")

g19 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = Mean.TPM..Wildtype.)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM (WT) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("TPM")

g20 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = Mean.TPM..ZnGroup.)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM (ZnGroup) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("TPM")

g21 <- ggplot(data = final_summary, aes(x = TIGRFAM.main.role, y = log2FC)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the cumulative log2FC (all genes) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC (all genes)")

g22 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = log2FC)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the cumulative log2FC (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC (DEGs)")

# V2

setwd("../V2/")

final_summary_V2 <- read.csv(file = "summary.tsv", header = TRUE, sep = "\t", quote = "")
final_summary_V2_extended <- read.csv(file = "summary_V2.tsv", header = TRUE, sep = "\t", quote = "", dec = ",")

for (i in 1:nrow(final_summary_V2)){
  if (as.character(final_summary_V2$TIGRFAM.main.role[i]) == ""){
    final_summary_V2$TIGRFAM.main.role[i] <- "Unknown function"
  }
}

degs_V2 = final_summary_V2[final_summary_V2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. != 0, ]
tpm_values_V2 <- read.csv(file = "filtered_gene_counts_tpm.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)
tpm_values_V2$TIGRFAM <- final_summary_V2$TIGRFAM.main.role

# Heatmaps

# Apply log2FC cutoff to narrow down list of targets
high_FC_V2 <- degs_V2[abs(degs_V2$log2FC) > 2, ]
tpm_degs_V2 <- tpm_values_V2[tpm_values_V2$ID %in% high_FC_V2$ID, ]
tpm_degs_V2$gene.name <- replace(tpm_degs_V2$gene.name, tpm_degs_V2$gene.name == "", "hypothetical_gene")
tpm_degs_V2_ordered <- tpm_degs_V2[order(tpm_degs_V2$TIGRFAM, tpm_degs_V2$gene.name), ]
rownames(tpm_degs_V2_ordered) <- make.names(tpm_degs_V2_ordered$gene.name, unique = TRUE)

groups_V2 = data.frame(Groups = c(rep("Control_Group", 3), rep("Zn_Group", 3)))
rownames(groups_V2) <- c("X5974_C2", "X5974_C3", "X5974_C6", "X5974_Z2", "X5974_Z3", "X5974_Z6")

anno_colours_V2 = list(Groups = c(Control_Group = "#96ca00", Zn_Group = "#ff82cb"))

g23 <- pheatmap(tpm_degs_V2_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_V2_ordered[10],
               annotation_col = groups_V2, annotation_colors = anno_colours_V2, fontsize_row = 9, fontsize_col = 9,
               cellheight = 10, cellwidth = 20,
               border_color = "black", gaps_row = c(2, 5, 10, 15, 17, 19, 21, 22, 29),
               main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
               scale = "none", breaks = seq(0, 10000, 100))

g24 <- pheatmap(tpm_degs_V2_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_V2_ordered[10],
               annotation_col = groups_V2, annotation_colors = anno_colours_V2, fontsize_row = 9, fontsize_col = 9,
               cellheight = 10, cellwidth = 20,
               border_color = "black", gaps_row = c(2, 5, 10, 15, 17, 19, 21, 22, 29),
               main = "Heatmap of the gene expression (TPM)\nbetween sample groups (scaled)\n",
               scale = "row")

# Cross-checking with summary table
# Removing any hypothetical proteins
final_summary_V1_short <- data.frame(gene.name = final_summary$gene.name, DE = final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup., log2FC = final_summary$log2FC, Mean_TPM_Base = final_summary$Mean.TPM..Wildtype., Mean_TPM_V1 = final_summary$Mean.TPM..ZnGroup., TIGRFAM = final_summary$TIGRFAM.main.role)
final_summary_V1_short <- final_summary_V1_short[final_summary_V1_short$gene.name != "", ]
final_summary_V1_short <- final_summary_V1_short[final_summary_V1_short$DE == "1" | final_summary_V1_short$DE == "-1", ]
final_summary_V1_short <- final_summary_V1_short[abs(final_summary_V1_short$log2FC) > 2, ]

final_summary_V2_short <- data.frame(gene.name = final_summary_V2$gene.name, DE = final_summary_V2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup., log2FC = final_summary_V2$log2FC, Mean_TPM_Base = final_summary_V2$Mean.TPM..Wildtype., Mean_TPM_V2 = final_summary_V2$Mean.TPM..ZnGroup., TIGRFAM = final_summary_V2$TIGRFAM.main.role)
final_summary_V2_short <- final_summary_V2_short[final_summary_V2_short$gene.name != "", ]
final_summary_V2_short <- final_summary_V2_short[final_summary_V2_short$DE == "1" | final_summary_V2_short$DE == "-1", ]
final_summary_V2_short <- final_summary_V2_short[abs(final_summary_V2_short$log2FC) > 2, ]

merged_df <- merge(final_summary_V1_short, final_summary_V2_short, by = "gene.name")
merged_df <- merged_df[merged_df$DE.x == "1" | merged_df$DE.x == "-1" | merged_df$DE.y == "1" | merged_df$DE.y == "-1", ]
merged_df <- merged_df[abs(merged_df$log2FC.x) > 2 | abs(merged_df$log2FC.y) > 2, ]
#merged_df <- data.frame(gene.name = merged_df$gene.name, TPM1 = merged_df$log2FC.x, TPM2 = merged_df$log2FC.y, TIGRFAM = merged_df$TIGRFAM.x)
merged_df <- data.frame(gene.name = merged_df$gene.name, TPM1 = merged_df$Mean_TPM_V1, TPM2 = merged_df$Mean_TPM_V2, TIGRFAM = merged_df$TIGRFAM.x)
rownames(merged_df) <- merged_df$gene.name
merged_df <- merged_df[order(merged_df$TIGRFAM, merged_df$gene.name), ]

groups_test = data.frame(Groups = c(rep("Control_Group", 3), rep("Zn_Group", 3)))
rownames(groups_test) <- c("X5974_C2", "X5974_C3", "X5974_C6", "X5974_Z2", "X5974_Z3", "X5974_Z6")

anno_colours_test = list(Groups = c(Control_Group = "#96ca00", Zn_Group = "#ff82cb"))

g25 <- pheatmap(merged_df[, c("TPM1", "TPM2")], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = merged_df[4],
                annotation_col = groups_test, annotation_colors = anno_colours_test, fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(1, 2, 4, 5, 9),
                main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
                scale = "none", breaks = seq(-5, 5, 0.25))

# Updated attempt
final_summary_V1_cut <- final_summary[final_summary$gene.name != "", ]
final_summary_V2_cut <- final_summary_V2[final_summary_V2$gene.name != "", ]
merged_df_2 <- merge(final_summary_V1_cut, final_summary_V2_cut, by = "gene.name")
merged_df_2 <- data.frame(gene.name = merged_df_2$gene.name, log2FC_V1 = merged_df_2$log2FC.x, DE_V1 = merged_df_2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..x, TPM_WT_V1 = merged_df_2$Mean.TPM..Wildtype..x, TPM_V1 = merged_df_2$Mean.TPM..ZnGroup..x, log2FC_V2 = merged_df_2$log2FC.y, DE_V2 = merged_df_2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..y, TPM_WT_V2 = merged_df_2$Mean.TPM..Wildtype..y, TPM_V2 = merged_df_2$Mean.TPM..ZnGroup..y, TIGRFAM = merged_df_2$TIGRFAM.main.role.x)
merged_df_2 <- merged_df_2[merged_df_2$DE_V1 == "1" | merged_df_2$DE_V1 == "-1" | merged_df_2$DE_V2 == "1" | merged_df_2$DE_V2 == "-1", ]
merged_df_2 <- merged_df_2[abs(merged_df_2$log2FC_V1) > 2 | abs(merged_df_2$log2FC_V2) > 2, ]
merged_df_3 <- data.frame(gene.name = merged_df_2$gene.name, log2FC_V1 = merged_df_2$log2FC_V1, log2FC_V2 = merged_df_2$log2FC_V2, "C_5974" = merged_df_2$TPM_WT_V2, "C_6122" = merged_df_2$TPM_WT_V1, "Z_5974" = merged_df_2$TPM_V2, "Z_6122" = merged_df_2$TPM_V1, TIGRFAM = merged_df_2$TIGRFAM)
rownames(merged_df_3) <- merged_df_3$gene.name
merged_df_3 <- merged_df_3[order(merged_df_3$TIGRFAM, merged_df_3$gene.name), ]

merged_df_2 <- data.frame(gene.name = merged_df_2$gene.name, log2FC_V1 = merged_df_2$log2FC_V1, log2FC_V2 = merged_df_2$log2FC_V2, "C_5974_6122" = ((merged_df_2$TPM_WT_V1 + merged_df_2$TPM_WT_V2) / 2), "Z_5974" = merged_df_2$TPM_V2, "Z_6122" = merged_df_2$TPM_V1, TIGRFAM = merged_df_2$TIGRFAM)
rownames(merged_df_2) <- merged_df_2$gene.name
merged_df_2 <- merged_df_2[order(merged_df_2$TIGRFAM, merged_df_2$gene.name), ]

groups_test_new = data.frame(Groups = c(rep("Control_Group", 1), rep("Zn_Group", 2)))
rownames(groups_test_new) <- c("C_5974_6122", "Z_5974", "Z_6122")

anno_colours_test_new = list(Groups = c(Control_Group = "#96ca00", Zn_Group = "#ff82cb"))

g26 <- pheatmap(merged_df_2[, c("C_5974_6122", "Z_5974", "Z_6122")], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = merged_df_2[7],
                annotation_col = groups_test_new, annotation_colors = anno_colours_test_new, fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
                main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
                scale = "none", breaks = seq(0, 5000, 50))

groups_test_new_new = data.frame(Groups = c(rep("Control_Group", 2), rep("Zn_Group", 2)))
rownames(groups_test_new_new) <- c("C_5974" , "C_6122", "Z_5974", "Z_6122")

g27 <- pheatmap(merged_df_3[, c("C_5974", "C_6122", "Z_5974", "Z_6122")], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = merged_df_3[8],
                annotation_col = groups_test_new_new, annotation_colors = anno_colours_test_new, fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
                main = "Heatmap of the\ngene expression (TPM)\nbetween sample groups\n",
                scale = "none", breaks = seq(0, 5000, 50))

g28 <- pheatmap(merged_df_3[, c("Z_5974", "Z_6122")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_3[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
                main = "Heatmap of the gene\nexpression (TPM)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(0, 5000, 50))

merged_df_4 <- merged_df_3
colnames(merged_df_4) <- c("gene.name", "Z_5974", "Z_6122", "TPM1", "TPM2", "TPM3", "TPM4", "TIGRFAM")
tmp_1 <- merged_df_4$Z_5974
tmp_2 <- merged_df_4$Z_6122
merged_df_4$Z_5974 <- tmp_2
merged_df_4$Z_6122 <- tmp_1

g29 <- pheatmap(merged_df_4[, c("Z_5974", "Z_6122")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_4[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
                main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57))

# Updated heatmap from 02.10.2019
# Switched position of columns
# TO-DO: Include genes only found in one of the strains
# TO-DO: Include hypothetical genes (e.g. name them by ID)
# TO-DO: Move exporting part of code to the exporting function below

#final_summary_V1_known <- final_summary[final_summary$gene.name != "", ]
#final_summary_V2_known <- final_summary_V2[final_summary_V2$gene.name != "", ]

# Problem with filtering before merging: might loose potential genes
# Merge first?
final_summary_V1_unknown <- final_summary[which(final_summary$gene.name == "" & (final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. == "1" | final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. == "-1") & abs(final_summary$log2FC) > 2), ]
final_summary_V2_unknown <- final_summary_V2[which(final_summary_V2$gene.name == "" & (final_summary_V2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. == "1" | final_summary_V2$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. == "-1") & abs(final_summary_V2$log2FC) > 2), ]
merged_df_5 <- merge(final_summary_V1_unknown, final_summary_V2_unknown, by = "nucleotide.sequence", all = TRUE)
#temp_names <- sprintf("hypothetical_%s", seq(1:nrow(merged_df_5)))
temp_names <- merged_df_5$ID.y
merged_df_6 <- data.frame(gene.name = temp_names,
                          "Z_6122" = merged_df_5$log2FC.x,
                          "Z_5974" = merged_df_5$log2FC.y,
                          "TPM1" = merged_df_5$Mean.TPM..Wildtype..y,
                          "TPM2" = merged_df_5$Mean.TPM..Wildtype..x,
                          "TPM3" = merged_df_5$Mean.TPM..ZnGroup..y,
                          "TPM4" = merged_df_5$Mean.TPM..ZnGroup..x,
                          TIGRFAM = merged_df_5$TIGRFAM.main.role.x)
rownames(merged_df_6) <- temp_names
# This will need to be manually updated for new data
merged_df_6$TIGRFAM = c("Unknown function", "Unknown function", "Unknown function", "Unknown function", "Purines, pyrimidines, nucleosides, and nucleotides")
merged_df_7 <- rbind(merged_df_6, merged_df_4)
merged_df_7 <- merged_df_7[order(merged_df_7$TIGRFAM, tolower(merged_df_7$gene.name)), ]

g30 <- pheatmap(merged_df_7[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_7[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 31, 38),
                main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57))

t1 <- merge(final_summary, final_summary_V2, by = "nucleotide.sequence", all = TRUE)
t2 <- t1[which((abs(t1$log2FC.x) > 2 & (t1$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..x == "1" | t1$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..x == "-1")) | (abs(t1$log2FC.y) > 2) & (t1$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..y == "1" | t1$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup..y == "-1")), ]
t3 <- data.frame(t2$gene.name.x, t2$gene.name.y)

v <- as.vector(t3$t2.gene.name.x)
v[3] <- "borD_1"
v[51] <- "puuA"
v[64] <- "borD_2"
v[65] <- "cysA"
v[66] <- "ompF"
v[67] <- "hpf"
v[68] <- "zinT"
v[69] <- "pgaB"
# Check if these IDs are correct
v[v == ""] <- as.character(temp_names)

merged_df_8 <- data.frame(gene.name = v,
                          "Z_6122" = t2$log2FC.x,
                          "Z_5974" = t2$log2FC.y,
                          "TPM1" = t2$Mean.TPM..Wildtype..y,
                          "TPM2" = t2$Mean.TPM..Wildtype..x,
                          "TPM3" = t2$Mean.TPM..ZnGroup..y,
                          "TPM4" = t2$Mean.TPM..ZnGroup..x,
                          TIGRFAM = t2$TIGRFAM.main.role.x)
merged_df_8$TIGRFAM[is.na(merged_df_8$TIGRFAM)] <- "Unknown function"
merged_df_8 <- merged_df_8[order(merged_df_8$TIGRFAM, tolower(merged_df_8$gene.name)), ]

# zinT is duplicated for some reason...
merged_df_8[67, 3] <- merged_df_8[68, 3]
merged_df_8[67, 4] <- merged_df_8[68, 4]
merged_df_8[67, 6] <- merged_df_8[68, 6]
merged_df_8[68, ] <- merged_df_8[69, ]
merged_df_8 <- merged_df_8[1:68, ]

rownames(merged_df_8) <- merged_df_8$gene.name

g31 <- pheatmap(merged_df_8[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_8[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 11, 17, 23, 25, 29, 32, 33, 39),
                main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57))

merged_df_9 <- merged_df_8
#merged_df_9[merged_df_9$gene.name == "yeaG",] <- c("yeaG", 2.0487, 1.1456, NA, NA, NA, NA, "Protein fate")
#merged_df_9[merged_df_9$gene.name == "borD_2",] <- c("borD_2", -3.6996, -3.169, NA, NA, NA, NA, "Unknown function")
#merged_df_9[merged_df_9$gene.name == "cysA",] <- c("cysA", -1.2269, -2.2662, NA, NA, NA, NA, "Unknown function")
#merged_df_9[merged_df_9$gene.name == "hpf",] <- c("hpf", -0.017, 2.3889, NA, NA, NA, NA, "Unknown function")
#merged_df_9[merged_df_9$gene.name == "ompF",] <- c("ompF", -1.1785, -4.6426, NA, NA, NA, NA, "Unknown function")
#merged_df_9[merged_df_9$gene.name == "pgaB",] <- c("pgaB", -0.0532, -2.5346, NA, NA, NA, NA, "Unknown function")

merged_df_9[29, 3] <- 1.1456
merged_df_9[42, 2] <- -3.6996
merged_df_9[45, 2] <- -1.2269
merged_df_9[52, 2] <- -0.017
merged_df_9[57, 2] <- -1.1785
merged_df_9[58, 2] <- -0.0532

colfunc <- colorRampPalette(c("blue", "aliceblue"))
col1 <- colfunc(11)
col2 <- "white"
colfunc <- colorRampPalette(c("pink", "red3"))
col3 <- colfunc(11)
colours_new <- c(col1, col2, col3)

g32 <- pheatmap(merged_df_9[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_9[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 11, 17, 23, 25, 29, 32, 33, 39),
                main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(-5.5, 5.5, 0.5), color = colours_new)

# Scatter Plots

temp_df <- data.frame(genes = c(as.character(merged_df_9$gene.name), as.character(merged_df_9$gene.name)), log_FC = c(merged_df_9$Z_6122, merged_df_9$Z_5974), Groups = c(rep("Z_6122", 68), rep("Z_5974", 68)), Regression = c(rep("Z_6122", 68), rep("Z_5974", 68)), TIGRFAM = c(as.character(merged_df_9$TIGRFAM), as.character(merged_df_9$TIGRFAM)), gene_number = c(seq(1, 68, 1), seq(1, 68, 1)))

ggplot(temp_df, aes(x = gene_number, y = (log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2))

ggsave("cloud01.png", width = 15, height = 10)

ggplot(temp_df, aes(x = gene_number, y = (log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = ""), method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2))

ggsave("cloud02.png", width = 15, height = 10)

ggplot(temp_df, aes(x = gene_number, y = abs(log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("absolute log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02))

ggsave("cloud03.png", width = 15, height = 10)

ggplot(temp_df, aes(x = gene_number, y = abs(log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("absolute log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = ""), method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02))

ggsave("cloud04.png", width = 15, height = 10)

ggplot(temp_df, aes(x = gene_number, y = abs(log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("absolute log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(6.5, 10.5, 11.5, 17.5, 23.5, 25.5, 29.5, 32.5, 33.5, 39.5), linetype = "dashed", color = "lightgrey", size = 1)

ggsave("cloud05.png", width = 15, height = 10)

ggplot(temp_df, aes(x = gene_number, y = log_FC, shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("relative log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(6.5, 10.5, 11.5, 17.5, 23.5, 25.5, 29.5, 32.5, 33.5, 39.5), linetype = "dashed", color = "lightgrey", size = 1)

ggsave("cloud06.png", width = 15, height = 10)

new_data <- read.csv(file = "summary_V1_V2_filtered.csv", header = TRUE, sep = "\t", quote = "")

# New Revision of data from Vanessa

temp_df_2 <- temp_df

temp_genes <- as.character(temp_df_2$genes)
temp_genes[33] <- "hiuH"
temp_genes[101] <- "hiuH"
temp_genes[53] <- "yciG"
temp_genes[121] <- "yciG"
temp_genes[54] <- "virK"
temp_genes[122] <- "virK"
temp_genes[56] <- "yebE "
temp_genes[124] <- "yebE"

temp_df_2$genes <- temp_genes

temp_tigrfams <- as.character(temp_df_2$TIGRFAM)

i = 0

for (gene in temp_df_2$genes){
  i = i + 1
  if (gene %in% new_data$gene.name){
    tigrfam <- as.character(new_data[new_data$gene.name == gene, 4])
    if (is.na(tigrfam)){
      temp_tigrfams[i] <- "Unknown function"
    } else{
      if (tigrfam == "Unknown"){
        temp_tigrfams[i] <- "Unknown function"
      } else{
        temp_tigrfams[i] <- tigrfam
      }
    }
  }
}

temp_df_2$TIGRFAM <- as.character(temp_tigrfams)

temp_df_2 <- temp_df_2[order(temp_df_2$TIGRFAM, temp_df_2$genes), ]

l = c(seq(1, 68, 1), seq(1, 68, 1))
l = l[order(l)]

temp_df_2$gene_number <- l

ggplot(temp_df_2, aes(x = gene_number, y = abs(log_FC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("absolute log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02))

# WARNING - ISSUE WITH GENE YEBE - WRONG TIGRFAM

ggsave("cloud07.png", width = 15, height = 10)

merged_df_10 <- merged_df_9

merged_genes <- as.character(merged_df_10$gene.name)
merged_genes[33] <- "hiuH"
merged_genes[53] <- "yciG"
merged_genes[54] <- "virK"
merged_genes[56] <- "yebE "

merged_df_10$gene.name <- merged_genes
rownames(merged_df_10) <- merged_genes

merged_tigrfams <- as.character(merged_df_10$TIGRFAM)

i = 0

for (gene in merged_df_10$gene.name){
  i = i + 1
  if (gene %in% new_data$gene.name){
    tigrfam <- as.character(new_data[new_data$gene.name == gene, 4])
    if (is.na(tigrfam)){
      merged_tigrfams[i] <- "Unknown function"
    } else{
      if (tigrfam == "Unknown"){
        merged_tigrfams[i] <- "Unknown function"
      } else{
        merged_tigrfams[i] <- tigrfam
      }
    }
  }
}

merged_df_10$TIGRFAM <- as.character(merged_tigrfams)

merged_df_10 <- merged_df_10[order(merged_df_10$TIGRFAM, merged_df_10$gene.name), ]

g33 <- pheatmap(merged_df_10[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_10[8],
                fontsize_row = 9, fontsize_col = 9,
                cellheight = 10, cellwidth = 20,
                border_color = "black", gaps_row = c(6, 10, 15, 16, 21, 28, 30, 36, 41, 42, 47, 64, 67),
                main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
                scale = "none", breaks = seq(-5.5, 5.5, 0.5), color = colours_new, file = "new_pheatmap_12.png")

genes = c()
logFC = c()
groups = c()
TIGRFAM = c()
gene_number = c()
n = 0

#for (i in final_summary$ID){
#  n = n + 1
#  V1_gene_name = as.character(final_summary$gene.name[final_summary$ID == i])
#  V1_nucleotide_seq = as.character(final_summary$nucleotide.sequence[final_summary$ID == i])
#  for (j in final_summary_V2$ID){
#    V2_gene_name = as.character(final_summary_V2$gene.name[final_summary_V2$ID == j])
#    V2_nucleotide_seq = as.character(final_summary_V2$nucleotide.sequence[final_summary_V2$ID == j])
#    if (((V1_gene_name == V2_gene_name) && (V1_gene_name != "")) || V1_nucleotide_seq == V2_nucleotide_seq){
#      genes = c(genes, V1_gene_name, V1_gene_name)
#      logFC = c(logFC, final_summary$log2FC[final_summary$ID == i], final_summary_V2$log2FC[final_summary_V2$ID == j])
#      groups = c(groups, "Z_6122", "Z_5974")
#      TIGRFAM = c(TIGRFAM, final_summary$TIGRFAM.main.role, final_summary$TIGRFAM.main.role)
#      gene_number = c(gene_number, n, n)
#    }
#  }
#}

#regression = groups

#final_df = data.frame(Genes = genes, logFC = logFC, Groups = groups, Regression = regression, TIGRFAM = TIGRAM, gene_number = gene_number)

#ggplot(final_df, aes(x = gene_number, y = abs(logFC), shape = Groups, colour = TIGRFAM)) +
#  geom_point() +
#  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) +
#  xlab("genes") +
#  ylab("absolute log2 FC") +
#  ggtitle("Z_5974 vs. Z_6122") +
#  geom_smooth(aes(linetype = Regression, group = ""), method = "lm", formula = y ~ x, se = FALSE) +
#  scale_x_continuous(breaks = seq(1, 68), labels = temp_df$genes[1:68], expand = c(0, 0.2)) +
#  scale_y_continuous(expand = c(0, 0.02))

#ggsave("cloud10.png", width = 15, height = 10)

final_summary_V1_V2_merged <- merge(final_summary_V1_extended, final_summary_V2_extended, by = "nucleotide.sequence")
final_genes = c(as.character(final_summary_V1_V2_merged$gene.name.x), as.character(final_summary_V1_V2_merged$gene.name.x))
final_number_genes = nrow(final_summary_V1_V2_merged)
final_logFC = c(final_summary_V1_V2_merged$log2FC.x, final_summary_V1_V2_merged$log2FC.y)
final_groups = c(rep("Z_6122", final_number_genes), rep("Z_5974", final_number_genes))
final_regression = final_groups
final_TIGRFAM = c(as.character(final_summary_V1_V2_merged$TIGRFAM.main.role.x), as.character(final_summary_V1_V2_merged$TIGRFAM.main.role.y))
final_TIGRFAM[final_TIGRFAM == ""] <- "Unknown function"
final_gene_number = c(seq(1, final_number_genes, 1), seq(1, final_number_genes, 1))
final_summary_V1_V2_output <- data.frame(Genes = final_genes, logFC = final_logFC, Groups = final_groups, Regression = final_regression, TIGRFAM = final_TIGRFAM, gene_number = final_gene_number)
final_summary_V1_V2_output <- final_summary_V1_V2_output[order(final_summary_V1_V2_output$TIGRFAM, final_summary_V1_V2_output$Genes), ]
s = c(seq(1, final_number_genes, 1), seq(1, final_number_genes, 1))
s = s[order(s)]
final_summary_V1_V2_output$gene_number <- s

ggplot(final_summary_V1_V2_output, aes(x = gene_number, y = abs(logFC), shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("absolute log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, final_number_genes), labels = final_summary_V1_V2_output$Genes[1:final_number_genes], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(93.5, 252.5, 364.5, 518.5, 587.5, 711.5, 1031.5, 1069.5, 1097.5, 1125.5, 1287.5, 1471.5, 1531.5, 1669.5, 1692.5, 1726.5, 2103.5), linetype = "dashed", color = "lightgrey", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 0.6)

ggsave("cloud08.png", width = 15, height = 10)

ggplot(final_summary_V1_V2_output, aes(x = gene_number, y = logFC, shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("relative log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, final_number_genes), labels = final_summary_V1_V2_output$Genes[1:final_number_genes], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(93.5, 252.5, 364.5, 518.5, 587.5, 711.5, 1031.5, 1069.5, 1097.5, 1125.5, 1287.5, 1471.5, 1531.5, 1669.5, 1692.5, 1726.5, 2103.5), linetype = "dashed", color = "lightgrey", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1)

ggsave("cloud09.png", width = 15, height = 10)

ggplot(final_summary_V1_V2_output[final_summary_V1_V2_output$logFC < 0, ], aes(x = gene_number, y = logFC, shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("relative log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, final_number_genes), labels = final_summary_V1_V2_output$Genes[1:final_number_genes], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(93.5, 252.5, 364.5, 518.5, 587.5, 711.5, 1031.5, 1069.5, 1097.5, 1125.5, 1287.5, 1471.5, 1531.5, 1669.5, 1692.5, 1726.5, 2103.5), linetype = "dashed", color = "lightgrey", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 0.6)

ggsave("cloud10.png", width = 15, height = 10)

ggplot(final_summary_V1_V2_output[final_summary_V1_V2_output$logFC > 0, ], aes(x = gene_number, y = logFC, shape = Groups, colour = TIGRFAM)) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("genes") +
  ylab("relative log2 FC") +
  ggtitle("Z_5974 vs. Z_6122") +
  geom_smooth(aes(linetype = Regression, group = Groups), method = "loess", formula = y ~ x, se = FALSE) +
  scale_x_continuous(breaks = seq(1, final_number_genes), labels = final_summary_V1_V2_output$Genes[1:final_number_genes], expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  geom_vline(xintercept = c(93.5, 252.5, 364.5, 518.5, 587.5, 711.5, 1031.5, 1069.5, 1097.5, 1125.5, 1287.5, 1471.5, 1531.5, 1669.5, 1692.5, 1726.5, 2103.5), linetype = "dashed", color = "lightgrey", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 0.6)

ggsave("cloud11.png", width = 15, height = 10)

# NEW IMAGES DECEMBER

colfunc <- colorRampPalette(c("#1BB6E4", "#02131A"))
col1 <- colfunc(9)
col2 <- "black"
colfunc <- colorRampPalette(c("#1B1B00", "#FDFE00"))
col3 <- colfunc(9)
colours_new_new <- c(col1, col2, col3)

temp_df_3 <- data.frame(Z_6122 = final_summary_V1_V2_output$logFC[final_summary_V1_V2_output$Groups == "Z_6122"],
                       Z_5974 = final_summary_V1_V2_output$logFC[final_summary_V1_V2_output$Groups == "Z_5974"],
                       gene_number = final_summary_V1_V2_output$gene_number[final_summary_V1_V2_output$Groups == "Z_5974"])

colnames(temp_df_3) <- c("6122", "5974")

#temp_df_3 <- temp_df_3[order(rowsum(abs(tmp_df_3$Z_5974), abs(tmp_df_3$Z_6122))), ]
#tmp_df_new <- tmp_df_new[as.integer(rownames(tmp_df_new)) < 10, ]

test_1 <- pheatmap(temp_df_3[, c("6122", "5974")],
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 annotation_row = gene_number,
                 show_rownames = FALSE,
                 #fontsize_row = 9,
                 #fontsize_col = 9,
                 #cellheight = 10,
                 #cellwidth = 20,
                 #border_color = "black",
                 #gaps_row = c(6, 10, 15, 16, 21, 28, 30, 36, 41, 42, 47, 64, 67),
                 main = "Heatmap of all DEGs",
                 scale = "none",
                 breaks = seq(-4.5, 4.5, 0.5),
                 color = colours_new_new,
                 file = "pheatmap_dark_all_genes.png")

colnames(temp_df_3) <- c("Z_6122", "Z_5974", "gene_number")

down <- temp_df_3[1:2]
up <- temp_df_3[1:2]

down[down$Z_6122 > 0, 1] <- 0
down[down$Z_5974 > 0, 2] <- 0
down[down$Z_6122 < 0, 1] <- 1
down[down$Z_5974 < 0, 2] <- 1

up[up$Z_6122 > 0, 1] <- 1
up[up$Z_5974 > 0, 2] <- 1
up[up$Z_6122 < 0, 1] <- 0
up[up$Z_5974 < 0, 2] <- 0

v1 <- vennCounts(down)
v2 <- vennCounts(up)

png(filename = "venn_diagram_downregulated.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
vennDiagram(v1, circle.col = c("red", "blue"))
dev.off()

png(filename = "venn_diagram_upregulated.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
vennDiagram(v2, circle.col = c("red", "blue"))
dev.off()

colnames(temp_df_3) <- c("6122", "5974", "gene_number")

temp_df_3_filtered <- temp_df_3[abs(temp_df_3$"6122") > 1 | abs(temp_df_3$"5974") > 1, ]

test_2 <- pheatmap(temp_df_3_filtered[, c("6122", "5974")],
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 annotation_row = gene_number,
                 show_rownames = FALSE,
                 #fontsize_row = 9,
                 #fontsize_col = 9,
                 #cellheight = 10,
                 #cellwidth = 20,
                 #border_color = "black",
                 #gaps_row = c(6, 10, 15, 16, 21, 28, 30, 36, 41, 42, 47, 64, 67),
                 main = "Heatmap of filtered DEGs",
                 scale = "none",
                 breaks = seq(-4.5, 4.5, 0.5),
                 color = colours_new_new,
                 file = "pheatmap_dark_filtered_genes.png")

# NEW IMAGES FEBRUARY

setwd("../../deg/")
R1 <- read.csv(file = "I1/summary.tsv", header = TRUE, sep = "\t", quote = "")
R2 <- read.csv(file = "I2/summary.tsv", header = TRUE, sep = "\t", quote = "")
R3 <- read.csv(file = "I3/summary.tsv", header = TRUE, sep = "\t", quote = "")
R4 <- read.csv(file = "I4/summary.tsv", header = TRUE, sep = "\t", quote = "")

R1_DE <- R1[R1$DE..SCORE....1..Input...Mock..1..Input...Mock. != 0 & abs(R1$log2FC) > 2, ]
R2_DE <- R2[R2$DE..SCORE....1..Input...WT..1..Input...WT. != 0 & abs(R2$log2FC) > 2, ]
R3_DE <- R3[R3$DE..SCORE....1..Input...dXCL1..1..Input...dXCL1. != 0 & abs(R3$log2FC) > 2, ]
R4_DE <- R4[R4$DE..SCORE....1..Input...UV..1..Input...UV. != 0 & abs(R4$log2FC) > 2, ]

# Ignore R1 since we do not focus on that
DEGs <- unique(c(as.character(R1_DE$ID), c(as.character(R2_DE$ID), as.character(R3_DE$ID), as.character(R4_DE$ID))))

f1 <- R1[R1$ID %in% DEGs, ]
f2 <- R2[R2$ID %in% DEGs, ]
f3 <- R3[R3$ID %in% DEGs, ]
f4 <- R4[R4$ID %in% DEGs, ]

d1 <- data.frame(ID = f1$ID, Gene = f1$gene.name, log2FC = f1$log2FC, DE = f1$DE..SCORE....1..Input...Mock..1..Input...Mock., TPM = f1$Presence.Absence..Input.)
d2 <- data.frame(ID = f2$ID, Gene = f2$gene.name, log2FC = f2$log2FC, DE = f2$DE..SCORE....1..Input...WT..1..Input...WT., TPM = f2$Presence.Absence..WT.)
d3 <- data.frame(ID = f3$ID, Gene = f3$gene.name, log2FC = f3$log2FC, DE = f3$DE..SCORE....1..Input...dXCL1..1..Input...dXCL1., TPM = f3$Presence.Absence..dXCL1.)
d4 <- data.frame(ID = f4$ID, Gene = f4$gene.name, log2FC = f4$log2FC, DE = f4$DE..SCORE....1..Input...UV..1..Input...UV., TPM = f4$Presence.Absence..UV.)

m1 <- merge(d1, d2, by = "ID", all = TRUE)
m2 <- merge(d3, d4, by = "ID", all = TRUE)
m3 <- merge(m1, m2, by = "ID", all = TRUE)

m4 <- data.frame(ID = m3$ID, Gene = m3$Gene.x.y, log2FC_R1 = m3$log2FC.x.x, log2FC_R2 = m3$log2FC.y.x, log2FC_R3 = m3$log2FC.x.y, log2FC_R4 = m3$log2FC.y.y, DE_R1 = m3$DE.x.x, DE_R2 = m3$DE.y.x, DE_R3 = m3$DE.x.y, DE_R4 = m3$DE.y.y, TPM_R1 = m3$TPM.x.x, TPM_R2 = m3$TPM.y.x, TPM_R3 = m3$TPM.x.y, TPM_R4 = m3$TPM.y.y)
m4[is.na(m4)] <- 0

names <- c("Mock", "RCMV-E wt", expression(paste("RCMV-E ", Delta, italic("vxcl1"), sep = "")), "UV")

pheatmap(m4[3:6],
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         file = "pheatmap_degs.png",
         main = "Heatmap of filtered DEGs",
         scale = "none",
         show_rownames = FALSE,
         labels_col = names,
         labels_row = m4$Gene,
         cellheight = 0.25,
         cellwidth = 50
)

m4_down <- m4
m4_up <- m4

m4_down <- m4_down[7:10]
m4_down[m4_down > 0 ] <- 0
m4_down[m4_down < 0 ] <- 1

m4_up <- m4_up[7:10]
m4_up[m4_up < 0 ] <- 0

v3 <- vennCounts(m4_down)
v4 <- vennCounts(m4_up)

png(filename = "venn_diagram_downregulated.png", width = 25, height = 25, units = "cm", res = 500, pointsize = 20)
vennDiagram(v3, circle.col = c("#0808AF", "#a00000", "#f2b77d", "#addfae"), main = "Downregulated Genes", names = names, cex = 1.1)
dev.off()

png(filename = "venn_diagram_upregulated.png", width = 25, height = 25, units = "cm", res = 500, pointsize = 20)
vennDiagram(v4, circle.col = c("#0808AF", "#a00000", "#f2b77d", "#addfae"), main = "Upregulated Genes", names = names, cex = 1.1)
dev.off()

# Print overlap of experiments
#m4$Gene[as.integer(rownames(m4_down[m4_down$DE_R1 > 0 & m4_down$DE_R2 > 0 & m4_down$DE_R3 > 0 & m4_down$DE_R4 > 0, ]))]
#m4$Gene[as.integer(rownames(m4_up[m4_up$DE_R1 > 0 & m4_up$DE_R2 > 0 & m4_up$DE_R3 > 0 & m4_up$DE_R4 > 0, ]))]

r1_specific_down <- m4[rownames(m4_down[m4_down$DE_R1 > 0 & m4_down$DE_R2 == 0 & m4_down$DE_R3 == 0 & m4_down$DE_R4 == 0, ]), ]
r1_specific_up <- m4[rownames(m4_up[m4_up$DE_R1 > 0 & m4_up$DE_R2 == 0 & m4_up$DE_R3 == 0 & m4_up$DE_R4 == 0, ]), ]
r1_degs <- rbind(r1_specific_down, r1_specific_up)

r2_specific_down <- m4[rownames(m4_down[m4_down$DE_R1 == 0 & m4_down$DE_R2 > 0 & m4_down$DE_R3 == 0 & m4_down$DE_R4 == 0, ]), ]
r2_specific_up <- m4[rownames(m4_up[m4_up$DE_R1 == 0 & m4_up$DE_R2 > 0 & m4_up$DE_R3 == 0 & m4_up$DE_R4 == 0, ]), ]
r2_degs <- rbind(r2_specific_down, r2_specific_up)

r3_specific_down <- m4[rownames(m4_down[m4_down$DE_R1 == 0 & m4_down$DE_R2 == 0 & m4_down$DE_R3 > 0 & m4_down$DE_R4 == 0, ]), ]
r3_specific_up <- m4[rownames(m4_up[m4_up$DE_R1 == 0 & m4_up$DE_R2 == 0 & m4_up$DE_R3 > 0 & m4_up$DE_R4 == 0, ]), ]
r3_degs <- rbind(r3_specific_down, r3_specific_up)

r4_specific_down <- m4[rownames(m4_down[m4_down$DE_R1 == 0 & m4_down$DE_R2 == 0 & m4_down$DE_R3 == 0 & m4_down$DE_R4 > 0, ]), ]
r4_specific_up <- m4[rownames(m4_up[m4_up$DE_R1 == 0 & m4_up$DE_R2 == 0 & m4_up$DE_R3 == 0 & m4_up$DE_R4 > 0, ]), ]
r4_degs <- rbind(r4_specific_down, r4_specific_up)

# R05

#R5 <- read.csv(file = "R05-2020-02-07-11-18-270e8/summary.tsv", header = TRUE, sep = "\t", quote = "")
#R5_DE <- R5[R5$DE..SCORE....1..WT...dXCL..1..WT...dXCL. != 0 & abs(R5$log2FC) > 0, ]

#pheatmap(R5_DE$log2FC,
#         border_color = "black",
#         cellheight = 10,
#         cellwidth = 20,
#         cluster_cols = FALSE,
#         cluster_rows = TRUE,
#         file = "pheatmap_r5.png",
#         fontsize_col = 5,
#         fontsize_row = 5,
#         main = "Heatmap (R5)",
#         show_rownames = TRUE,
#         labels_col = names[3],
#         labels_row = R5$gene.name
#)

# Visualization of genes of interest
genes_of_interest <- data.frame(Gene = c("Cd80", "Cd86", "Cd40", "Ccr7", "Akt3", "Pik3ca", "Mapk3", "Ccl2", "Ccl17", "Ccl6", "Rab27a", "Sec22b", "Tap1", "Tap2", "Tapbp", "Calr", "Pdia3", "Rac2", "Cybb", "Snap23"),
                                Group = c(rep("Maturation Markers", 4), rep("Migration Markers", 3), rep("Chemokines", 3), rep("Antigen Presentation", 10)))
genes_of_interest <- genes_of_interest[order(genes_of_interest$Gene), ]

R1_GOI <- R1[R1$gene.name %in% genes_of_interest$Gene, ]
R2_GOI <- R2[R2$gene.name %in% genes_of_interest$Gene, ]
R3_GOI <- R3[R3$gene.name %in% genes_of_interest$Gene, ]
R4_GOI <- R4[R4$gene.name %in% genes_of_interest$Gene, ]

R1_GOI <- R1_GOI[order(R1_GOI$gene.name), ]
R2_GOI <- R2_GOI[order(R2_GOI$gene.name), ]
R3_GOI <- R3_GOI[order(R3_GOI$gene.name), ]
R4_GOI <- R4_GOI[order(R4_GOI$gene.name), ]

df_GOI <- data.frame(Gene = R1_GOI$gene.name, Group = genes_of_interest$Group, log2FC_R1 = R1_GOI$log2FC, log2FC_R2 = R2_GOI$log2FC, log2FC_R3 = R3_GOI$log2FC, log2FC_R4 = R4_GOI$log2FC, TPM_R1 = R1_GOI$Mean.TPM..Mock., TPM_R2 = R2_GOI$Mean.TPM..WT., TPM_R3 = R3_GOI$Mean.TPM..dXCL1., TPM_R4 = R4_GOI$Mean.TPM..UV.)
df_GOI <- df_GOI[order(df_GOI$Group), ]
rownames(df_GOI) <- seq(1:nrow(genes_of_interest))

mycolor <- colorRampPalette(c("blue", "white", "red"))(50)
mybreaks <- c(seq(min(df_GOI[4:6]), 0, length.out = 25),
              seq(max(df_GOI[4:6])/50, max(df_GOI[4:6]), length.out = 25)
)

#pheatmap(df_GOI[3:6],
#         annotation_row = df_GOI[2],
#         border_color = "black",
#         #breaks = mybreaks,
#         cellheight = 50,
#         cellwidth = 100,
#         cluster_cols = TRUE,
#         cluster_rows = FALSE,
#         #color = mycolor,
#         file = "pheatmap_genes_of_interest_c1.png",
#         fontsize_col = 10,
#         fontsize_row = 10,
#         #fontsize = 10,
#         #gaps_row = c(4, 7),
#         main = "Heatmap of genes of interest",
#         #scale = "row",
#         show_rownames = TRUE,
#         labels_col = names,
#         labels_row = df_GOI$Gene
#)

pheatmap(df_GOI[3:6],
         annotation_row = df_GOI[2],
         border_color = "black",
         breaks = mybreaks,
         cellheight = 50,
         cellwidth = 50,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         color = mycolor,
         file = "pheatmap_genes_of_interest_c1.png",
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize = 12,
         #gaps_row = c(4, 7),
         main = "Heatmap of genes of interest",
         #scale = "row",
         show_rownames = TRUE,
         labels_col = names,
         labels_row = df_GOI$Gene
)

df_GOI_TPM <- data.frame(Genes = rep(df_GOI$Gene, 5), TPM = c(R1_GOI$Mean.TPM..Input., df_GOI$TPM_R1, df_GOI$TPM_R2, df_GOI$TPM_R3, df_GOI$TPM_R4), Sample = c(rep("Input", nrow(genes_of_interest)), rep("Mock", nrow(genes_of_interest)), rep("RCMV-E wt", nrow(genes_of_interest)), rep("RCMV-E dvxcl1", nrow(genes_of_interest)), rep("UV", nrow(genes_of_interest))))

ggplot(data = df_GOI_TPM, aes(x = Genes, y = log(TPM+0.8, 10), fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_manual(labels = c("Input", "Mock", bquote(paste("RCMV-E ", Delta, italic("vxcl1"))), "RCMV-E wt", "UV"), values = c("#ff7e00", "#010080", "#f2b77d", "#a00000", "#addfae")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Expression of selected genes (RNA-Seq)", x = "Genes of interest")

ggsave("bar_chart_01.png", width = 15, height = 5, dpi = 500)

#ggplot(data = df_GOI_TPM, aes(x = Genes, y = TPM, fill = Sample)) +
#  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
#  scale_fill_manual(labels = c("Input", "Mock", bquote(paste("RCMV-E ", Delta, italic("vxcl1"))), "RCMV-E wt", "UV"), values = c("#ff7e00", "#010080", "#f2b77d", "#a00000", "#addfae")) +
#  theme_minimal() +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  labs(title = "Expression of selected genes (RNA-Seq)", x = "Genes of interest")

#ggsave("bar_chart_02.png", width = 15, height = 5, dpi = 500)

# Visualization of pathway of interest
pathway_of_interest_g1 = c("Xcl1", "Xcr1",
                           "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7",
                           "Ccl17", "Ccl21", "Ccl22", "Ccl24", "Ccl28", "Ccr1",
                           "Ccr1l1", "Ccr3", "Ccr4", "Ccr7", "Ccr8", "Ccr9", "Ccr10",
                           "Cxcl1", "Cxcl2", "Cxcl9", "Cxcl16", "Cxcl17", "Cxcr2",
                           "Cxcr3", "Cxcr4", "Cxcr5", "Pf4", "Ppbp",
                           "Cx3cr1",
                           "Csf2", "Csf2ra", "Csf2rb", "Il3ra", "Il13", "Il13ra1",
                           "Clcf1", "Cntf", "Il6", "Il6r", "Il6st", "Il11", "Il12a", "Il12b", "Il12rb1",
                           "Il23a", "Il23r", "Il27", "Il31", "Il31ra", "Lif", "Osm", "Osmr",
                           "Ghr", "Lepr", "Prlr",
                           "Il2", "Il2ra", "Il2rg", "Il4", "Il7r", "Il9",
                           "Il9r", "Il15", "Il21", "Tslp",
                           "Il10", "Il10ra", "Il10rb", "Il19", "Il20ra", "Il20rb", "Il22", "Il22ra1", "Il24", "Ifnl3",
                           "Ifnar2", "Ifnb1", "Ifne", "Ifng", "Ifngr1", "Ifngr2",
                           "Il1a", "Il1b", "Il18",
                           "Il17a", "Il17c", "Il17d", "Il17f", "Il17ra", "Il17rb",
                           "Csf1r"
)

pathway_of_interest_g2 = c(rep("Chemokines - C Subfamily", 2),
                           rep("Chemokines - CC Subfamily", 20),
                           rep("Chemokines - CXC Subfamily", 11),
                           rep("Chemokines - CX3C Subfamily", 1),
                           rep("Cytokines - Class I Helical - IL4-like", 6),
                           rep("Cytokines - Class I Helical - IL6/12-like", 17),
                           rep("Cytokines - Class I Helical - Prolactin Family", 3),
                           rep("Cytokines - Class I Helical - Y-Chain Utilising", 10),
                           rep("Cytokines - Class II Helical - IL10/28-like", 10),
                           rep("Cytokines - Class II Helical - Interferon Family", 6),
                           rep("Cytokines - IL1-like", 3),
                           rep("Cytokines - IL17-like", 6),
                           rep("Non-classified", 1)
)

R1[R1$gene.name %in% "Csf1r", 2]
R2[R2$gene.name %in% "Csf1r", 2]
R3[R3$gene.name %in% "Csf1r", 2]
R4[R4$gene.name %in% "Csf1r", 2]

genes_of_interest <- data.frame(Gene = pathway_of_interest_g1, Group = pathway_of_interest_g2)
genes_of_interest <- genes_of_interest[order(genes_of_interest$Gene), ]

R1_GOI <- R1[R1$gene.name %in% genes_of_interest$Gene, ]
R2_GOI <- R2[R2$gene.name %in% genes_of_interest$Gene, ]
R3_GOI <- R3[R3$gene.name %in% genes_of_interest$Gene, ]
R4_GOI <- R4[R4$gene.name %in% genes_of_interest$Gene, ]

R1_GOI <- R1_GOI[order(R1_GOI$gene.name), ]
R2_GOI <- R2_GOI[order(R2_GOI$gene.name), ]
R3_GOI <- R3_GOI[order(R3_GOI$gene.name), ]
R4_GOI <- R4_GOI[order(R4_GOI$gene.name), ]

# Check if genes are available in all results
genes_of_interest[!(genes_of_interest$Gene %in% R1_GOI$gene.name), 1]
genes_of_interest[!(genes_of_interest$Gene %in% R2_GOI$gene.name), 1]
genes_of_interest[!(genes_of_interest$Gene %in% R3_GOI$gene.name), 1]
genes_of_interest[!(genes_of_interest$Gene %in% R4_GOI$gene.name), 1]

df_GOI <- data.frame(Gene = R1_GOI$gene.name, Group = genes_of_interest$Group, log2FC_R1 = R1_GOI$log2FC, log2FC_R2 = R2_GOI$log2FC, log2FC_R3 = R3_GOI$log2FC, log2FC_R4 = R4_GOI$log2FC, TPM_R1 = R1_GOI$Mean.TPM..Mock., TPM_R2 = R2_GOI$Mean.TPM..WT., TPM_R3 = R3_GOI$Mean.TPM..dXCL1., TPM_R4 = R4_GOI$Mean.TPM..UV.)
df_GOI <- df_GOI[order(df_GOI$Group), ]
rownames(df_GOI) <- seq(1:nrow(genes_of_interest))

mycolor <- colorRampPalette(c("blue", "white", "red"))(50)
mybreaks <- c(seq(min(df_GOI[4:6]), 0, length.out = 25),
              seq(max(df_GOI[4:6])/50, max(df_GOI[4:6]), length.out = 25)
)

pheatmap(df_GOI[3:6],
         annotation_row = df_GOI[2],
         border_color = "black",
         breaks = mybreaks,
         cellheight = 30,
         cellwidth = 50,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         color = mycolor,
         file = "pheatmap_genes_of_interest_c2.png",
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize = 12,
         #gaps_row = c(4, 7),
         main = "Heatmap of pathway of interest (rno04060)",
         #scale = "row",
         show_rownames = TRUE,
         labels_col = names,
         labels_row = df_GOI$Gene
)

# NEW IMAGES MARCH

setwd("../../deg/")
BA1 <- read.csv(file = "BA1/summary.tsv", header = TRUE, sep = "\t", quote = "")
BA2 <- read.csv(file = "BA2/summary.tsv", header = TRUE, sep = "\t", quote = "")
BC1 <- read.csv(file = "BC1/summary.tsv", header = TRUE, sep = "\t", quote = "")
BC2 <- read.csv(file = "BC2/summary.tsv", header = TRUE, sep = "\t", quote = "")

DE_Up_BA1 <- BA1[BA1$DE..SCORE....1..O2...CO2..1..O2...CO2. == 1 & abs(BA1$log2FC) >= 2, ]
DE_Up_BA2 <- BA2[BA2$DE..SCORE....1..O2...CO2..1..O2...CO2. == 1 & abs(BA2$log2FC) >= 2, ]
DE_Up_BC1 <- BC1[BC1$DE..SCORE....1..O2...CO2..1..O2...CO2. == 1 & abs(BC1$log2FC) >= 2, ]
DE_Up_BC2 <- BC2[BC2$DE..SCORE....1..O2...CO2..1..O2...CO2. == 1 & abs(BC2$log2FC) >= 2, ]

DE_Down_BA1 <- BA1[BA1$DE..SCORE....1..O2...CO2..1..O2...CO2. == -1 & abs(BA1$log2FC) >= 2, ]
DE_Down_BA2 <- BA2[BA2$DE..SCORE....1..O2...CO2..1..O2...CO2. == -1 & abs(BA2$log2FC) >= 2, ]
DE_Down_BC1 <- BC1[BC1$DE..SCORE....1..O2...CO2..1..O2...CO2. == -1 & abs(BC1$log2FC) >= 2, ]
DE_Down_BC2 <- BC2[BC2$DE..SCORE....1..O2...CO2..1..O2...CO2. == -1 & abs(BC2$log2FC) >= 2, ]

M1 <- merge(DE_Up_BA1, DE_Up_BA2, by = "AA.sequence", all = TRUE)
M2 <- merge(DE_Up_BC1, DE_Up_BC2, by = "AA.sequence", all = TRUE)
M3 <- merge(M1, M2, by = "AA.sequence", all = TRUE)

# Start with final merge here...

B1 <- DE_Up_BA1
B2 <- DE_Up_BA2
B3 <- DE_Up_BC1
B4 <- DE_Up_BC2

Gene_Symbol_BA1 = c()
Gene_Symbol_BA2 = c()
Gene_Symbol_BC1 = c()
Gene_Symbol_BC2 = c()
ID_BA1 = c()
ID_BA2 = c()
ID_BC1 = c()
ID_BC2 = c()
log2FC_BA1 = c()
log2FC_BA2 = c()
log2FC_BC1 = c()
log2FC_BC2 = c()
TPM_BA1 = c()
TPM_BA2 = c()
TPM_BC1 = c()
TPM_BC2 = c()
DE_BA1 = c()
DE_BA2 = c()
DE_BC1 = c()
DE_BC2 = c()

# BA1

i = 1
for (id in B1$ID) {
  t2 <- NA
  t3 <- NA
  t4 <- NA
  B2_ID <- ""
  B2_gene_name <- ""
  B2_fc <- 0
  B2_de <- 0
  B2_tpm <- 0
  B3_ID <- ""
  B3_gene_name <- ""
  B3_fc <- 0
  B3_de <- 0
  B3_tpm <- 0
  B4_ID <- ""
  B4_gene_name <- ""
  B4_fc <- 0
  B4_de <- 0
  B4_tpm <- 0

  B1_ID <- id
  B1_gene_name <- as.character(B1[B1$ID == id, 2])
  B1_fc <- B1[B1$ID == id, 4]
  B1_de <- B1[B1$ID == id, 5]
  B1_tpm <- B1[B1$ID == id, 15]
  B1_NS <- B1[B1$ID == id, 20]
  
  t2 <- amatch(B1_NS, B2$nucleotide.sequence, maxDist = length(B1_NS)*0.1)
  t3 <- amatch(B1_NS, B3$nucleotide.sequence, maxDist = length(B1_NS)*0.1)
  t4 <- amatch(B1_NS, B4$nucleotide.sequence, maxDist = length(B1_NS)*0.1)
  
  if (is.na(t2) == FALSE) {
    B2_ID <- as.character(B2[t2, 1])
    B2_gene_name <- as.character(B2[t2, 2])
    B2_fc <- B2[t2, 4]
    B2_de <- B2[t2, 5]
    B2_tpm <- B2[t2, 15]
    B2 <- B2[-c(t2), ]
  }

  if (is.na(t3) == FALSE) {
    B3_ID <- as.character(B3[t3, 1])
    B3_gene_name <- as.character(B3[t3, 2])
    B3_fc <- B3[t3, 4]
    B3_de <- B3[t3, 5]
    B3_tpm <- B3[t3, 15]
    B3 <- B3[-c(t3), ]
  }

  if (is.na(t4) == FALSE) {
    B4_ID <- as.character(B4[t4, 1])
    B4_gene_name <- as.character(B4[t4, 2])
    B4_fc <- B4[t4, 4]
    B4_de <- B4[t4, 5]
    B4_tpm <- B4[t4, 15]
    B4 <- B4[-c(t4), ]
  }
  
  Gene_Symbol_BA1[i] <- B1_gene_name
  Gene_Symbol_BA2[i] <- B2_gene_name
  Gene_Symbol_BC1[i] <- B3_gene_name
  Gene_Symbol_BC2[i] <- B4_gene_name
  ID_BA1[i] <- B1_ID
  ID_BA2[i] <- B2_ID
  ID_BC1[i] <- B3_ID
  ID_BC2[i] <- B4_ID
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1[i] <- B1_tpm
  TPM_BA2[i] <- B2_tpm
  TPM_BC1[i] <- B3_tpm
  TPM_BC2[i] <- B4_tpm
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de

  i = i + 1
}

# BA2

for (id in B2$ID) {
  t3 <- NA
  t4 <- NA
  B1_ID <- ""
  B1_gene_name <- ""
  B1_fc <- 0
  B1_de <- 0
  B1_tpm <- 0
  B3_ID <- ""
  B3_gene_name <- ""
  B3_fc <- 0
  B3_de <- 0
  B3_tpm <- 0
  B4_ID <- ""
  B4_gene_name <- ""
  B4_fc <- 0
  B4_de <- 0
  B4_tpm <- 0
  
  B2_ID <- id
  B2_gene_name <- as.character(B2[B2$ID == id, 2])
  B2_fc <- B2[B2$ID == id, 4]
  B2_de <- B2[B2$ID == id, 5]
  B2_tpm <- B2[B2$ID == id, 15]
  B2_NS <- B2[B2$ID == id, 20]
  
  t3 <- amatch(B2_NS, B3$nucleotide.sequence, maxDist = length(B2_NS)*0.1)
  t4 <- amatch(B2_NS, B4$nucleotide.sequence, maxDist = length(B2_NS)*0.1)
  
  if (is.na(t3) == FALSE) {
    B3_ID <- as.character(B3[t3, 1])
    B3_gene_name <- as.character(B3[t3, 2])
    B3_fc <- B3[t3, 4]
    B3_de <- B3[t3, 5]
    B3_tpm <- B3[t3, 15]
    B3 <- B3[-c(t3), ]
  }
  
  if (is.na(t4) == FALSE) {
    B4_ID <- as.character(B4[t4, 1])
    B4_gene_name <- as.character(B4[t4, 2])
    B4_fc <- B4[t4, 4]
    B4_de <- B4[t4, 5]
    B4_tpm <- B4[t4, 15]
    B4 <- B4[-c(t4), ]
  }
  
  Gene_Symbol_BA1[i] <- B1_gene_name
  Gene_Symbol_BA2[i] <- B2_gene_name
  Gene_Symbol_BC1[i] <- B3_gene_name
  Gene_Symbol_BC2[i] <- B4_gene_name
  ID_BA1[i] <- B1_ID
  ID_BA2[i] <- B2_ID
  ID_BC1[i] <- B3_ID
  ID_BC2[i] <- B4_ID
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1[i] <- B1_tpm
  TPM_BA2[i] <- B2_tpm
  TPM_BC1[i] <- B3_tpm
  TPM_BC2[i] <- B4_tpm
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de
  
  i = i + 1
}

# BC1

for (id in B3$ID) {
  t4 <- NA
  B1_ID <- ""
  B1_gene_name <- ""
  B1_fc <- 0
  B1_de <- 0
  B1_tpm <- 0
  B2_ID <- ""
  B2_gene_name <- ""
  B2_fc <- 0
  B2_de <- 0
  B2_tpm <- 0
  B4_ID <- ""
  B4_gene_name <- ""
  B4_fc <- 0
  B4_de <- 0
  B4_tpm <- 0
  
  B3_ID <- id
  B3_gene_name <- as.character(B3[B3$ID == id, 2])
  B3_fc <- B3[B3$ID == id, 4]
  B3_de <- B3[B3$ID == id, 5]
  B3_tpm <- B3[B3$ID == id, 15]
  B3_NS <- B3[B3$ID == id, 20]
  
  t4 <- amatch(B3_NS, B4$nucleotide.sequence, maxDist = length(B3_NS)*0.1)
  
  if (is.na(t4) == FALSE) {
    B4_ID <- as.character(B4[t4, 1])
    B4_gene_name <- as.character(B4[t4, 2])
    B4_fc <- B4[t4, 4]
    B4_de <- B4[t4, 5]
    B4_tpm <- B4[t4, 15]
    B4 <- B4[-c(t4), ]
  }
  
  Gene_Symbol_BA1[i] <- B1_gene_name
  Gene_Symbol_BA2[i] <- B2_gene_name
  Gene_Symbol_BC1[i] <- B3_gene_name
  Gene_Symbol_BC2[i] <- B4_gene_name
  ID_BA1[i] <- B1_ID
  ID_BA2[i] <- B2_ID
  ID_BC1[i] <- B3_ID
  ID_BC2[i] <- B4_ID
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1[i] <- B1_tpm
  TPM_BA2[i] <- B2_tpm
  TPM_BC1[i] <- B3_tpm
  TPM_BC2[i] <- B4_tpm
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de
  
  i = i + 1
}

# BC2

for (id in B4$ID) {
  B1_ID <- ""
  B1_gene_name <- ""
  B1_fc <- 0
  B1_de <- 0
  B1_tpm <- 0
  B2_ID <- ""
  B2_gene_name <- ""
  B2_fc <- 0
  B2_de <- 0
  B2_tpm <- 0
  B3_ID <- ""
  B3_gene_name <- ""
  B3_fc <- 0
  B3_de <- 0
  B3_tpm <- 0
  
  B4_ID <- id
  B4_gene_name <- as.character(B4[B4$ID == id, 2])
  B4_fc <- B4[B4$ID == id, 4]
  B4_de <- B4[B4$ID == id, 5]
  B4_tpm <- B4[B4$ID == id, 15]
  B4_NS <- B4[B4$ID == id, 20]
  
  Gene_Symbol_BA1[i] <- B1_gene_name
  Gene_Symbol_BA2[i] <- B2_gene_name
  Gene_Symbol_BC1[i] <- B3_gene_name
  Gene_Symbol_BC2[i] <- B4_gene_name
  ID_BA1[i] <- B1_ID
  ID_BA2[i] <- B2_ID
  ID_BC1[i] <- B3_ID
  ID_BC2[i] <- B4_ID
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1[i] <- B1_tpm
  TPM_BA2[i] <- B2_tpm
  TPM_BC1[i] <- B3_tpm
  TPM_BC2[i] <- B4_tpm
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de
  
  i = i + 1
}

M5 <- data.frame(Gene_BA1 = Gene_Symbol_BA1, Gene_BA2 = Gene_Symbol_BA2, Gene_BC1 = Gene_Symbol_BC1, Gene_BC2 = Gene_Symbol_BC2)
M5$Gene_BA1 <- ifelse(M5$Gene_BA1 == "", as.character(M5$Gene_BA2), as.character(M5$Gene_BA1))
M5$Gene_BA1 <- ifelse(M5$Gene_BA1 == "", as.character(M5$Gene_BC2), as.character(M5$Gene_BA1))
M5$Gene_BA1 <- ifelse(M5$Gene_BA1 == "", as.character(M5$Gene_BC1), as.character(M5$Gene_BA1))

M5$Gene_BA1[211] <- "dnaN"
M5$Gene_BA1[238] <- "nhaC"
M5$Gene_BA1[245] <- "proS"
M5$Gene_BA1[252] <- "opuD"
M5$Gene_BA1[348] <- "ansA"
M5$Gene_BA1[367] <- "ilvE"
M5$Gene_BA1[369] <- "ilvB_1"
M5$Gene_BA1[414] <- "abc-f"
M5$Gene_BA1[423] <- "prsA_2"
M5$Gene_BA1[431] <- "hisC"
M5$Gene_BA1[456] <- "fsa"
M5$Gene_BA1[476] <- "dhbF"
M5$Gene_BA1[496] <- "fabG"
M5$Gene_BA1[542] <- "tyrS"
M5$Gene_BA1[545] <- "hpt"
M5$Gene_BA1[616] <- "yidC"
M5$Gene_BA1[623] <- "pagA"
M5$Gene_BA1[629] <- "galU"
M5$Gene_BA1[635] <- "hfq"

M6 <- data.frame(BA1 = ID_BA1, BA2 = ID_BA2, BC1 = ID_BC1, BC2 = ID_BC2)
M6$BA1 <- ifelse(M6$BA1 == "", as.character(M6$BA2), as.character(M6$BA1))
M6$BA1 <- ifelse(M6$BA1 == "", as.character(M6$BC2), as.character(M6$BA1))
M6$BA1 <- ifelse(M6$BA1 == "", as.character(M6$BC1), as.character(M6$BA1))

M4_old <- data.frame(BA1_ID = ID_BA1, BA2_ID = ID_BA2, BC1_ID = ID_BC1, BC2_ID = ID_BC2,
                     Gene = M5$Gene_BA1,
                     log2FC_BA1 = log2FC_BA1, log2FC_BA2 = log2FC_BA2, log2FC_BC1 = log2FC_BC1, log2FC_BC2 = log2FC_BC2,
                     TPM_BA1 = TPM_BA1, TPM_BA2 = TPM_BA2, TPM_BC1 = TPM_BC1, TPM_BC2 = TPM_BC2,
                     DE_BA1 = DE_BA1, DE_BA2 = DE_BA2, DE_BC1 = DE_BC1, DE_BC2 = DE_BC2)

M4 <- data.frame(ID = M6$BA1, Gene = M5$Gene_BA1,
                 log2FC_BA1 = log2FC_BA1, log2FC_BA2 = log2FC_BA2, log2FC_BC1 = log2FC_BC1, log2FC_BC2 = log2FC_BC2,
                 TPM_BA1 = TPM_BA1, TPM_BA2 = TPM_BA2, TPM_BC1 = TPM_BC1, TPM_BC2 = TPM_BC2,
                 DE_BA1 = DE_BA1, DE_BA2 = DE_BA2, DE_BC1 = DE_BC1, DE_BC2 = DE_BC2)

# Venn Diagram
names_V5 <- c("Vollum", "Dobichau", "CA", "CI")
V5 <- vennCounts(M4[,11:14])
png(filename = "venn_diagram_upregulated.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
vennDiagram(V5, circle.col = c("red", "blue", "green", "grey"), main = "Upregulated Genes", names = names_V5, cex = 1.1)
dev.off()

# Upregulated_BA_Ids
Up_DEGs_BA <- M4_old[(M4_old$DE_BA1 == 1 & M4_old$DE_BA2 == 1) & (M4_old$DE_BC1 == 0 & M4_old$DE_BC2 == 0), ]
Up_DEGs_BA_short <- data.frame(Up_DEGs_BA$BA1_ID, Up_DEGs_BA$BA2_ID)
Up_DEGs_BA_short$Up_DEGs_BA.BA1_ID <- ifelse(Up_DEGs_BA_short$Up_DEGs_BA.BA1_ID == "", as.character(Up_DEGs_BA_short$Up_DEGs_BA.BA2_ID), as.character(Up_DEGs_BA_short$Up_DEGs_BA.BA1_ID))
Up_DEGs_BA_ID <- Up_DEGs_BA_short$Up_DEGs_BA.BA1_ID

# Upregulated_BC_Ids
Up_DEGs_BC <- M4_old[(M4_old$DE_BA1 == 0 & M4_old$DE_BA2 == 0) & (M4_old$DE_BC1 == 1 & M4_old$DE_BC2 == 1), ]
Up_DEGs_BC_short <- data.frame(Up_DEGs_BC$BC1_ID, Up_DEGs_BC$BC2_ID)
Up_DEGs_BC_short$Up_DEGs_BC.BC1_ID <- ifelse(Up_DEGs_BC_short$Up_DEGs_BC.BC1_ID == "", as.character(Up_DEGs_BC_short$Up_DEGs_BC.BC2_ID), as.character(Up_DEGs_BC_short$Up_DEGs_BC.BC1_ID))
Up_DEGs_BC_ID <- Up_DEGs_BC_short$Up_DEGs_BC.BC1_ID

a1 <- BA1[BA1$ID %in% Up_DEGs_BA_ID | BA1$ID %in% Up_DEGs_BC_ID, ]
a2 <- BA2[BA2$ID %in% Up_DEGs_BA_ID | BA2$ID %in% Up_DEGs_BC_ID, ]
a3 <- BC1[BC1$ID %in% Up_DEGs_BC_ID | BC1$ID %in% Up_DEGs_BC_ID, ]
a4 <- BC2[BC2$ID %in% Up_DEGs_BC_ID | BC2$ID %in% Up_DEGs_BC_ID, ]

C1 <- BA1
C2 <- BA2
C3 <- BC1
C4 <- BC2

ID <- c()
log2FC_BA1 <- c()
log2FC_BA2 <- c()
log2FC_BC1 <- c()
log2FC_BC2 <- c()
TPM_BA1_O2 <- c()
TPM_BA1_CO2 <- c()
TPM_BA2_O2 <- c()
TPM_BA2_CO2 <- c()
TPM_BC1_O2 <- c()
TPM_BC1_CO2 <- c()
TPM_BC2_O2 <- c()
TPM_BC2_CO2 <- c()
DE_BA1 <- c()
DE_BA2 <- c()
DE_BC1 <- c()
DE_BC2 <- c()
TIGRFAM_Main <- c()
TIGRFAM_Sub <- c()

i = 1
for (seq in a1$nucleotide.sequence) {
  b2 <- NA
  b3 <- NA
  b4 <- NA
  B1_fc <- C1[C1$nucleotide.sequence == seq, 4]
  B1_de <- C1[C1$nucleotide.sequence == seq, 5]
  B1_tpm_o2 <- C1[C1$nucleotide.sequence == seq, 14]
  B1_tpm_co2 <- C1[C1$nucleotide.sequence == seq, 15]
  B2_fc <- NA
  B2_de <- NA
  B2_tpm_o2 <- NA
  B2_tpm_co2 <- NA
  B3_fc <- NA
  B3_de <- NA
  B3_tpm_o2 <- NA
  B3_tpm_co2 <- NA
  B4_fc <- NA
  B4_de <- NA
  B4_tpm_o2 <- NA
  B4_tpm_co2 <- NA
  ID[i] <- as.character(C1[C1$nucleotide.sequence == seq, 1])
  TIGRFAM_Main[i] <- as.character(C1[C1$nucleotide.sequence == seq, 16])
  TIGRFAM_Sub[i] <- as.character(C1[C1$nucleotide.sequence == seq, 17])

  b2 <- amatch(seq, C2$nucleotide.sequence, maxDist = length(seq)*0.1)
  b3 <- amatch(seq, C3$nucleotide.sequence, maxDist = length(seq)*0.1)
  b4 <- amatch(seq, C4$nucleotide.sequence, maxDist = length(seq)*0.1)
  
  if (is.na(b2) == FALSE) {
    B2_fc <- C2[b2, 4]
    B2_de <- C2[b2, 5]
    B2_tpm_o2 <- C2[b2, 14]
    B2_tpm_co2 <- C2[b2, 15]
    C2 <- C2[-c(b2), ]
  }
  
  if (is.na(b3) == FALSE) {
    B3_fc <- C3[b3, 4]
    B3_de <- C3[b3, 5]
    B3_tpm_o2 <- C3[b3, 14]
    B3_tpm_co2 <- C3[b3, 15]
    C3 <- C3[-c(b3), ]
  }
  
  if (is.na(b4) == FALSE) {
    B4_fc <- C4[b4, 4]
    B4_de <- C4[b4, 5]
    B4_tpm_o2 <- C4[b4, 14]
    B4_tpm_co2 <- C4[b4, 15]
    C4 <- C4[-c(b4), ]
  }
  
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1_O2[i] <- B1_tpm_o2
  TPM_BA1_CO2[i] <- B1_tpm_co2
  TPM_BA2_O2[i] <- B2_tpm_o2
  TPM_BA2_CO2[i] <- B2_tpm_co2
  TPM_BC1_O2[i] <- B3_tpm_o2
  TPM_BC1_CO2[i] <- B3_tpm_co2
  TPM_BC2_O2[i] <- B4_tpm_o2
  TPM_BC2_CO2[i] <- B4_tpm_co2
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de
  
  i = i + 1
}

for (seq in a3$nucleotide.sequence) {
  b1 <- NA
  b2 <- NA
  b4 <- NA
  B1_fc <- NA
  B1_de <- NA
  B1_tpm_o2 <- NA
  B1_tpm_co2 <- NA
  B2_fc <- NA
  B2_de <- NA
  B2_tpm_o2 <- NA
  B2_tpm_co2 <- NA
  B3_fc <- C3[C3$nucleotide.sequence == seq, 4]
  B3_de <- C3[C3$nucleotide.sequence == seq, 5]
  B3_tpm_o2 <- C3[C3$nucleotide.sequence == seq, 14]
  B3_tpm_co2 <- C3[C3$nucleotide.sequence == seq, 15]
  B4_fc <- NA
  B4_de <- NA
  B4_tpm_o2 <- NA
  B4_tpm_co2 <- NA
  ID[i] <- as.character(C3[C3$nucleotide.sequence == seq, 1])
  TIGRFAM_Main[i] <- as.character(C3[C3$nucleotide.sequence == seq, 16])
  TIGRFAM_Sub[i] <- as.character(C3[C3$nucleotide.sequence == seq, 17])

  b1 <- amatch(seq, C1$nucleotide.sequence, maxDist = length(seq)*0.1)
  b2 <- amatch(seq, C2$nucleotide.sequence, maxDist = length(seq)*0.1)
  b4 <- amatch(seq, C4$nucleotide.sequence, maxDist = length(seq)*0.1)
  
  if (is.na(b1) == FALSE) {
    B1_fc <- C1[b1, 4]
    B1_de <- C1[b1, 5]
    B1_tpm_o2 <- C1[b1, 14]
    B1_tpm_co2 <- C1[b1, 15]
    C1 <- C1[-c(b1), ]
  }
  
  if (is.na(b2) == FALSE) {
    B2_fc <- C2[b2, 4]
    B2_de <- C2[b2, 5]
    B2_tpm_o2 <- C2[b2, 14]
    B2_tpm_co2 <- C2[b2, 15]
    C2 <- C2[-c(b2), ]
  }
  
  if (is.na(b4) == FALSE) {
    B4_fc <- C4[b4, 4]
    B4_de <- C4[b4, 5]
    B4_tpm_o2 <- C4[b4, 14]
    B4_tpm_co2 <- C4[b4, 15]
    C4 <- C4[-c(b4), ]
  }
  
  log2FC_BA1[i] <- B1_fc
  log2FC_BA2[i] <- B2_fc
  log2FC_BC1[i] <- B3_fc
  log2FC_BC2[i] <- B4_fc
  TPM_BA1_O2[i] <- B1_tpm_o2
  TPM_BA1_CO2[i] <- B1_tpm_co2
  TPM_BA2_O2[i] <- B2_tpm_o2
  TPM_BA2_CO2[i] <- B2_tpm_co2
  TPM_BC1_O2[i] <- B3_tpm_o2
  TPM_BC1_CO2[i] <- B3_tpm_co2
  TPM_BC2_O2[i] <- B4_tpm_o2
  TPM_BC2_CO2[i] <- B4_tpm_co2
  DE_BA1[i] <- B1_de
  DE_BA2[i] <- B2_de
  DE_BC1[i] <- B3_de
  DE_BC2[i] <- B4_de
  
  i = i + 1
}

M7 <- data.frame(ID = ID,
                 log2FC_BA1 = log2FC_BA1, log2FC_BA2 = log2FC_BA2, log2FC_BC1 = log2FC_BC1, log2FC_BC2 = log2FC_BC2,
                 TPM_BA1_O2 = TPM_BA1_O2, TPM_BA2_O2 = TPM_BA2_O2, TPM_BC1_O2 = TPM_BC1_O2, TPM_BC2_O2 = TPM_BC2_O2,
                 TPM_BA1_CO2 = TPM_BA1_CO2, TPM_BA2_CO2 = TPM_BA2_CO2, TPM_BC1_CO2 = TPM_BC1_CO2, TPM_BC2_CO2 = TPM_BC2_CO2,
                 DE_BA1 = DE_BA1, DE_BA2 = DE_BA2, DE_BC1 = DE_BC1, DE_BC2 = DE_BC2,
                 TIGRFAM_Main = TIGRFAM_Main, TIGRFAM_Sub = TIGRFAM_Sub)

# Pie Chart of Genes
M8_S1_Full <- M7[1:28, ]
M8_S2_Full <- M7[29:89, ]

M8_S1 <- (M8_S1_Full[M8_S1_Full$TIGRFAM_Main != "" & M8_S1_Full$TIGRFAM_Main != "Unknown function", ])
M8_S2 <- (M8_S2_Full[M8_S2_Full$TIGRFAM_Main != "" & M8_S2_Full$TIGRFAM_Main != "Unknown function", ])

M8_S1_counts <- data.frame(table(M8_S1$TIGRFAM_Main))
M8_S2_counts <- data.frame(table(M8_S2$TIGRFAM_Main))

ggplot(M8_S1_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (G1)", subtitle = "TIGRFAM Groups")

ggplot(M8_S2_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (G2)", subtitle = "TIGRFAM Groups")

# Upregulated in all conditions
Up_DEGs_All <- M4_old[M4_old$DE_BA1 == 1 & M4_old$DE_BA2 == 1 & M4_old$DE_BC1 == 1 & M4_old$DE_BC2 == 1, ]

# Heatmap 1

new_rownames <- M8_S1_Full$ID
j = 1
for (gene in a1$gene.name) {
  if (gene != "") {
    levels(new_rownames) <- c(levels(new_rownames), gene)
    new_rownames[j] <- gene
  }
  j = j + 1
}

new_rownames <- droplevels(new_rownames)

rownames(M8_S1_Full) <- new_rownames
names <- c("Vollum", "Dobichau", "CA", "CI", "Vollum", "Dobichau", "CA", "CI")
annotations <- data.frame(condition = c(rep("O2", 4), rep("CO2", 4)))
rownames(annotations) <- colnames(M8_S1_Full[6:13])
levels(M8_S1_Full$TIGRFAM_Main) <- c(levels(M8_S1_Full$TIGRFAM_Main), "Unknown function") 
M8_S1_Full$TIGRFAM_Main[M8_S1_Full$TIGRFAM_Main == ""] <- "Unknown function"
colnames(M8_S1_Full)[18] <- "TIGRFAM"
my_colour = list(
    condition = c(O2 = "chartreuse2", CO2 = "cornflowerblue"),
    TIGRFAM = c("Biosynthesis of cofactors, prosthetic groups, and carriers" = "#2f4f4f",
                "Cell envelope" = "#a0522d",
                "Cellular processes" = "#228b22",
                "Central intermediary metabolism" = "#4b0082",
                "DNA metabolism" = "#ff0000",
                "Energy metabolism" = "#ffff00",
                "Fatty acid and phospholipid metabolism" = "#00ff00",
                "Protein fate" = "#00ffff",
                "Protein synthesis" = "#0000ff",
                "Regulatory functions" = "#ff00ff",
                "Transcription" = "#ff69b4",
                "Transport and binding proteins" = "#6495ed",
                "Unknown function" = "#eee8aa")
)

pheatmap(log(M8_S1_Full[6:13] + 1),
         annotation_row = M8_S1_Full[18],
         annotation_col = annotations,
         annotation_colors = my_colour,
         border_color = "black",
         #breaks = mybreaks,
         cellheight = 10,
         cellwidth = 10,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = rev(heat.colors(10)),
         file = "pheatmap_genes_of_interest_g1.png",
         fontsize_col = 8,
         fontsize_row = 8,
         #gaps_row = c(5, 9, 15, 30),
         gaps_col = c(4),
         main = "Heatmap of genes of interest (I)",
         #scale = "row",
         show_rownames = TRUE,
         labels_col = names
         #labels_row = M8_S1_Full[1]
)

# Heatmap 2

new_rownames <- M8_S2_Full$ID
j = 1
for (gene in a3$gene.name) {
  if (gene != "") {
    levels(new_rownames) <- c(levels(new_rownames), gene)
    new_rownames[j] <- gene
  }
  j = j + 1
}

new_rownames <- droplevels(new_rownames)

rownames(M8_S2_Full) <- new_rownames
names <- c("Vollum", "Dobichau", "CA", "CI", "Vollum", "Dobichau", "CA", "CI")
annotations <- data.frame(condition = c(rep("O2", 4), rep("CO2", 4)))
rownames(annotations) <- colnames(M8_S2_Full[6:13])
levels(M8_S2_Full$TIGRFAM_Main) <- c(levels(M8_S2_Full$TIGRFAM_Main), "Unknown function") 
M8_S2_Full$TIGRFAM_Main[M8_S2_Full$TIGRFAM_Main == ""] <- "Unknown function"
colnames(M8_S2_Full)[18] <- "TIGRFAM"
my_colour = list(
    condition = c(O2 = "chartreuse2", CO2 = "cornflowerblue"),
    TIGRFAM = c("Biosynthesis of cofactors, prosthetic groups, and carriers" = "#2f4f4f",
                "Cell envelope" = "#a0522d",
                "Cellular processes" = "#228b22",
                "Central intermediary metabolism" = "#4b0082",
                "DNA metabolism" = "#ff0000",
                "Energy metabolism" = "#ffff00",
                "Fatty acid and phospholipid metabolism" = "#00ff00",
                "Protein fate" = "#00ffff",
                "Protein synthesis" = "#0000ff",
                "Regulatory functions" = "#ff00ff",
                "Transcription" = "#ff69b4",
                "Transport and binding proteins" = "#6495ed",
                "Unknown function" = "#eee8aa")
)

pheatmap(log(M8_S2_Full[6:13] + 1),
         annotation_row = M8_S2_Full[18],
         annotation_col = annotations,
         annotation_colors = my_colour,
         border_color = "black",
         #breaks = mybreaks,
         cellheight = 10,
         cellwidth = 10,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = rev(heat.colors(10)),
         file = "pheatmap_genes_of_interest_g2.png",
         fontsize_col = 8,
         fontsize_row = 8,
         #gaps_row = c(5, 9, 15, 30),
         gaps_col = c(4),
         main = "Heatmap of genes of interest (II)",
         #scale = "row",
         show_rownames = TRUE,
         labels_col = names
         #labels_row = M8_S2_Full[1]
)

# Upregulated in single conditions
TIGRFAM_new <- c()
k = 1
for (id in M4$ID) {
  if (id %in% BA1$ID) {
    TIGRFAM_new[k] <- as.character(BA1[BA1$ID == id, 16])
  }
  if (id %in% BA2$ID) {
    TIGRFAM_new[k] <- as.character(BA2[BA2$ID == id, 16])
  }
  if (id %in% BC1$ID) {
    TIGRFAM_new[k] <- as.character(BC1[BC1$ID == id, 16])
  }
  if (id %in% BC2$ID) {
    TIGRFAM_new[k] <- as.character(BC2[BC2$ID == id, 16])
  }
  k = k + 1
}

TIGRFAM_new[TIGRFAM_new == ""] <- "Unknown function"

M9 <- data.frame(M4, TIGRFAM = TIGRFAM_new)

M9_BA1 <- M9[M9$DE_BA1 == 1 & M9$DE_BA2 == 0 & M9$DE_BC1 == 0 & M9$DE_BC2 == 0 & M9$TIGRFAM != "Unknown function" & M9$TIGRFAM != "Unclassified", ]
M9_BA2 <- M9[M9$DE_BA1 == 0 & M9$DE_BA2 == 1 & M9$DE_BC1 == 0 & M9$DE_BC2 == 0 & M9$TIGRFAM != "Unknown function" & M9$TIGRFAM != "Unclassified", ]
M9_BC1 <- M9[M9$DE_BA1 == 0 & M9$DE_BA2 == 0 & M9$DE_BC1 == 1 & M9$DE_BC2 == 0 & M9$TIGRFAM != "Unknown function" & M9$TIGRFAM != "Unclassified", ]
M9_BC2 <- M9[M9$DE_BA1 == 0 & M9$DE_BA2 == 0 & M9$DE_BC1 == 0 & M9$DE_BC2 == 1 & M9$TIGRFAM != "Unknown function" & M9$TIGRFAM != "Unclassified", ]

M9_BA1_counts <- data.frame(table(M9_BA1$TIGRFAM))
M9_BA2_counts <- data.frame(table(M9_BA2$TIGRFAM))
M9_BC1_counts <- data.frame(table(M9_BC1$TIGRFAM))
M9_BC2_counts <- data.frame(table(M9_BC2$TIGRFAM))

ggplot(M9_BA1_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (Vollum)", subtitle = "TIGRFAM Groups")

ggplot(M9_BA2_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (Dobichau)", subtitle = "TIGRFAM Groups")

ggplot(M9_BC1_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (CA)", subtitle = "TIGRFAM Groups")

ggplot(M9_BC2_counts, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "TIGRFAM", title = "DEGs (CI)", subtitle = "TIGRFAM Groups")

# GeneSCF

BA1_KEGG <- read.csv(file = "BA1/pathway_analysis_KEGG/deg_gene_symbols_KEGG_banv_functional_classification.tsv", header = TRUE, sep = "\t", quote = "")
BA2_KEGG <- read.csv(file = "BA2/pathway_analysis_KEGG/deg_gene_symbols_KEGG_banv_functional_classification.tsv", header = TRUE, sep = "\t", quote = "")
BC1_KEGG <- read.csv(file = "BC1/pathway_analysis_KEGG/deg_gene_symbols_KEGG_bal_functional_classification.tsv", header = TRUE, sep = "\t", quote = "")
BC2_KEGG <- read.csv(file = "BC2/pathway_analysis_KEGG/deg_gene_symbols_KEGG_bal_functional_classification.tsv", header = TRUE, sep = "\t", quote = "")

BA1_KEGG <- BA1_KEGG[BA1_KEGG$Benjamini.and.Hochberg..FDR. <= 0.05, ]
BA2_KEGG <- BA2_KEGG[BA2_KEGG$Benjamini.and.Hochberg..FDR. <= 0.05, ]
BC1_KEGG <- BC1_KEGG[BC1_KEGG$Benjamini.and.Hochberg..FDR. <= 0.05, ]
BC2_KEGG <- BC2_KEGG[BC2_KEGG$Benjamini.and.Hochberg..FDR. <= 0.05, ]

BA1_KEGG$Process.name <- paste("kegg", str_split_fixed(BA1_KEGG$Process.name, "banv", 2)[, 2], sep = "")
BA2_KEGG$Process.name <- paste("kegg", str_split_fixed(BA2_KEGG$Process.name, "banv", 2)[, 2], sep = "")
BC1_KEGG$Process.name <- paste("kegg", str_split_fixed(BC1_KEGG$Process.name, "bal", 2)[, 2], sep = "")
BC2_KEGG$Process.name <- paste("kegg", str_split_fixed(BC2_KEGG$Process.name, "bal", 2)[, 2], sep = "")

Group_All <- BA1_KEGG[BA1_KEGG$Process.name %in% BA2_KEGG$Process.name & BA1_KEGG$Process.name %in% BC1_KEGG$Process.name & BA1_KEGG$Process.name %in% BC2_KEGG$Process.name, ]
Group_BA_Exclusive <- BA1_KEGG[BA1_KEGG$Process.name %in% BA2_KEGG$Process.name & !(BA1_KEGG$Process.name %in% BC1_KEGG$Process.name) & !(BA1_KEGG$Process.name %in% BC2_KEGG$Process.name), ]
Group_BC_Exclusive <- BC2_KEGG[BC2_KEGG$Process.name %in% BC1_KEGG$Process.name & !(BC2_KEGG$Process.name %in% BA1_KEGG$Process.name) & !(BC2_KEGG$Process.name %in% BA2_KEGG$Process.name), ]
BA1_Exclusive <- BA1_KEGG[!(BA1_KEGG$Process.name %in% BA2_KEGG$Process.name) & !(BA1_KEGG$Process.name %in% BC1_KEGG$Process.name) & !(BA1_KEGG$Process.name %in% BC2_KEGG$Process.name), ]
BA2_Exclusive <- BA2_KEGG[!(BA2_KEGG$Process.name %in% BA1_KEGG$Process.name) & !(BA2_KEGG$Process.name %in% BC1_KEGG$Process.name) & !(BA2_KEGG$Process.name %in% BC2_KEGG$Process.name), ]
BC1_Exclusive <- BC1_KEGG[!(BC1_KEGG$Process.name %in% BA1_KEGG$Process.name) & !(BC1_KEGG$Process.name %in% BA2_KEGG$Process.name) & !(BC1_KEGG$Process.name %in% BC2_KEGG$Process.name), ]
BC2_Exclusive <- BC2_KEGG[!(BC2_KEGG$Process.name %in% BA1_KEGG$Process.name) & !(BC2_KEGG$Process.name %in% BA2_KEGG$Process.name) & !(BC2_KEGG$Process.name %in% BC1_KEGG$Process.name), ]

tmp_data <- Group_All[order(Group_All$P.value), ]
tmp_data[, "Rank_in_order"] <- c(seq(1:length(tmp_data[, 1])))
tmp_data <- tmp_data[, cbind("Rank_in_order", "Process.name", "percentage.", "P.value")]
colnames(tmp_data) <- c("Rank", "Process", "Genes", "Pval")
tmp_data <- tmp_data[1:14, ]
newdata <- cbind(gsub( "~.*$", "", tmp_data[, "Process"]), tmp_data[, "Rank"])
colnames(newdata) <- c("IDs", "Rank")
png("pathway_analysis_enrichment_plot_all.png", width = 1200, height = 800)
tmp_data$Pval[tmp_data$Pval == 0] <- min(tmp_data$Pval[tmp_data$Pval > 0]) * 0.05
fill_data = paste0("[", sprintf("%02d", as.numeric(tmp_data$Rank)), "]\t", tmp_data[, "Process"])
plot <- ggplot(tmp_data, aes(x = Rank, y = -log10(Pval), size = Genes, label = newdata[, "IDs"], fill = fill_data, guide = FALSE)) +
  geom_point(colour = "#2E2E2E", shape = 21) +
  scale_size_area(max_size = 5) +
  labs(fill = "[Rank] Process", title = "Hits (All)") +
  scale_x_continuous(name = "Rank", limits = c(0, 15)) +
  scale_y_continuous(name = "-log10(Pval)", limits = c(1.1, max(-log10(tmp_data$Pval)) + 1)) +
  geom_text(size = 5, color = "#2E2E2E", hjust = -0.1, vjust = 0, angle = 45) +
  geom_hline(yintercept = 1.3) +
  theme(text = element_text(size = 14), plot.title = element_text(hjust = 0.5))
print(plot)
dev.off()

tmp_data <- Group_BA_Exclusive[order(Group_BA_Exclusive$P.value), ]
tmp_data[, "Rank_in_order"] <- c(seq(1:length(tmp_data[, 1])))
tmp_data <- tmp_data[, cbind("Rank_in_order", "Process.name", "percentage.", "P.value")]
colnames(tmp_data) <- c("Rank", "Process", "Genes", "Pval")
tmp_data <- tmp_data[1:19, ]
newdata <- cbind(gsub( "~.*$", "", tmp_data[, "Process"]), tmp_data[, "Rank"])
colnames(newdata) <- c("IDs", "Rank")
png("pathway_analysis_enrichment_plot_BA.png", width = 1200, height = 800)
tmp_data$Pval[tmp_data$Pval == 0] <- min(tmp_data$Pval[tmp_data$Pval > 0]) * 0.05
fill_data = paste0("[", sprintf("%02d", as.numeric(tmp_data$Rank)), "]\t", tmp_data[, "Process"])
plot <- ggplot(tmp_data, aes(x = Rank, y = -log10(Pval), size = Genes, label = newdata[, "IDs"], fill = fill_data, guide = FALSE)) +
  geom_point(colour = "#2E2E2E", shape = 21) +
  scale_size_area(max_size = 5) +
  labs(fill = "[Rank] Process", title = "Hits (BA)") +
  scale_x_continuous(name = "Rank", limits = c(0, 22)) +
  scale_y_continuous(name = "-log10(Pval)", limits = c(1.1, max(-log10(tmp_data$Pval)) + 1)) +
  geom_text(size = 5, color = "#2E2E2E", hjust = -0.1, vjust = 0, angle = 45) +
  geom_hline(yintercept = 1.3) +
  theme(text = element_text(size = 14), plot.title = element_text(hjust = 0.5))
print(plot)
dev.off()

tmp_data <- Group_BC_Exclusive[order(Group_BC_Exclusive$P.value), ]
tmp_data[, "Rank_in_order"] <- c(seq(1:length(tmp_data[, 1])))
tmp_data <- tmp_data[, cbind("Rank_in_order", "Process.name", "percentage.", "P.value")]
colnames(tmp_data) <- c("Rank", "Process", "Genes", "Pval")
tmp_data <- tmp_data[1:5, ]
newdata <- cbind(gsub( "~.*$", "", tmp_data[, "Process"]), tmp_data[, "Rank"])
colnames(newdata) <- c("IDs", "Rank")
png("pathway_analysis_enrichment_plot_BC.png", width = 1200, height = 800)
tmp_data$Pval[tmp_data$Pval == 0] <- min(tmp_data$Pval[tmp_data$Pval > 0]) * 0.05
fill_data = paste0("[", sprintf("%02d", as.numeric(tmp_data$Rank)), "]\t", tmp_data[, "Process"])
plot <- ggplot(tmp_data, aes(x = Rank, y = -log10(Pval), size = Genes, label = newdata[, "IDs"], fill = fill_data, guide = FALSE)) +
  geom_point(colour = "#2E2E2E", shape = 21) +
  scale_size_area(max_size = 5) +
  labs(fill = "[Rank] Process", title = "Hits (BC)") +
  scale_x_continuous(name = "Rank", limits = c(0, 22)) +
  scale_y_continuous(name = "-log10(Pval)", limits = c(1.1, max(-log10(tmp_data$Pval)) + 1)) +
  geom_text(size = 5, color = "#2E2E2E", hjust = -0.1, vjust = 0, angle = 45) +
  geom_hline(yintercept = 1.3) +
  theme(text = element_text(size = 14), plot.title = element_text(hjust = 0.5))
print(plot)
dev.off()

# Exporting

if (export == TRUE) {
  pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
           annotation_col = groups, annotation_colors = anno_colours, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
           main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
           scale = "none", breaks = seq(0, 10000, 100), file = "new_pheatmap_01.png")
  
  pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
           annotation_col = groups, annotation_colors = anno_colours, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
           main = "Heatmap of the gene expression (TPM)\nbetween sample groups (scaled)\n",
           scale = "row", file = "new_pheatmap_02.png")
  
  pheatmap(tpm_degs_V2_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_V2_ordered[10],
           annotation_col = groups_V2, annotation_colors = anno_colours_V2, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(2, 5, 10, 15, 17, 19, 21, 22, 29),
           main = "Heatmap of the gene expression (TPM)\nbetween sample groups\n",
           scale = "none", breaks = seq(0, 10000, 100), file = "new_pheatmap_03.png")
  
  pheatmap(tpm_degs_V2_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_V2_ordered[10],
           annotation_col = groups_V2, annotation_colors = anno_colours_V2, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(2, 5, 10, 15, 17, 19, 21, 22, 29),
           main = "Heatmap of the gene expression (TPM)\nbetween sample groups (scaled)\n",
           scale = "row", file = "new_pheatmap_04.png")
  
  pheatmap(merged_df_2[, c("C_5974_6122", "Z_5974", "Z_6122")], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = merged_df_2[7],
           annotation_col = groups_test_new, annotation_colors = anno_colours_test_new, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
           main = "Heatmap of the\ngene expression (TPM)\nbetween sample groups\n",
           scale = "none", breaks = seq(0, 5000, 50), file = "new_pheatmap_05.png")
  
  pheatmap(merged_df_3[, c("C_5974", "C_6122", "Z_5974", "Z_6122")], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = merged_df_3[8],
           annotation_col = groups_test_new_new, annotation_colors = anno_colours_test_new, fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
           main = "Heatmap of the\ngene expression (TPM)\nbetween sample groups\n",
           scale = "none", breaks = seq(0, 5000, 50), file = "new_pheatmap_06.png")
  
  pheatmap(merged_df_3[, c("Z_5974", "Z_6122")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_3[8],
           fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
           main = "Heatmap of the gene\nexpression (TPM)\nbetween\nsample groups\n",
           scale = "none", breaks = seq(0, 5000, 50), file = "new_pheatmap_07.png")
  
  pheatmap(merged_df_4[, c("Z_5974", "Z_6122")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_4[8],
           fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 37),
           main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
           scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57), file = "new_pheatmap_08.png")
  
  pheatmap(merged_df_7[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_7[8],
          fontsize_row = 9, fontsize_col = 9,
          cellheight = 10, cellwidth = 20,
          border_color = "black", gaps_row = c(6, 10, 15, 21, 23, 27, 30, 31, 38),
          main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
          scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57), file = "new_pheatmap_09.png")
  
  pheatmap(merged_df_8[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_8[8],
          fontsize_row = 9, fontsize_col = 9,
          cellheight = 10, cellwidth = 20,
          border_color = "black", gaps_row = c(6, 10, 11, 17, 23, 25, 29, 32, 33, 39),
          main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
          scale = "none", breaks = seq(-5.5, 5.5, 0.2), color = plasma(57), file = "new_pheatmap_10.png")
  
  pheatmap(merged_df_9[, c("Z_6122", "Z_5974")], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = merged_df_9[8],
           fontsize_row = 9, fontsize_col = 9,
           cellheight = 10, cellwidth = 20,
           border_color = "black", gaps_row = c(6, 10, 11, 17, 23, 25, 29, 32, 33, 39),
           main = "Heatmap of the gene\nexpression (log2FC)\nbetween\nsample groups\n",
           scale = "none", breaks = seq(-5.5, 5.5, 0.5), color = colours_new, file = "new_pheatmap_11.png")
  
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  print(g5)
  print(g6)
  print(g7)
  print(g8)
  print(g9)
  print(g10)
  print(g11)
  print(g12)
  print(g13)
  print(g14)
  print(g15)
  print(g16)
  print(g17)
  print(g18)
  print(g19)
  print(g20)
  print(g21)
  print(g22)
  print(g23)
  print(g24)
  print(g25)
  print(g26)
  print(g27)
  print(g28)
  print(g29)
  print(g30)
  print(g31)
  print(g32)
  print(g33)
  
  dev.off()
}
