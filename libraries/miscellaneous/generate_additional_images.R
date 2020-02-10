# --------------------------------------------
# Title: generate_additional_images.R
# Author: Silver A. Wolf
# Last Modified: Mo, 10.02.2020
# Version: 0.2.7
# --------------------------------------------

# This script is used to generate additional images for publications, etc.
# It will need to be heavily adapted for any new data and is used to test new plots
# Finished plots are moved to the main SCORE scripts during development

# Libraries

library("ggplot2")
library("limma")
library("pheatmap")
library("reshape2")
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
R1 <- read.csv(file = "R01-2020-01-30-18-42-ecd9d/summary.tsv", header = TRUE, sep = "\t", quote = "")
R2 <- read.csv(file = "R02-2020-01-31-20-24-4a096/summary.tsv", header = TRUE, sep = "\t", quote = "")
R3 <- read.csv(file = "R03-2020-02-01-20-48-c2356/summary.tsv", header = TRUE, sep = "\t", quote = "")
R4 <- read.csv(file = "R04-2020-02-02-18-52-d31c7/summary.tsv", header = TRUE, sep = "\t", quote = "")

R1_DE <- R1[R1$DE..SCORE....1..Mock...Input..1..Mock...Input. != 0 & abs(R1$log2FC) > 3, ]
R2_DE <- R2[R2$DE..SCORE....1..Mock...WT..1..Mock...WT. != 0 & abs(R2$log2FC) > 3, ]
R3_DE <- R3[R3$DE..SCORE....1..Mock...dXCL1..1..Mock...dXCL1. != 0 & abs(R3$log2FC) > 3, ]
R4_DE <- R4[R4$DE..SCORE....1..Mock...UV..1..Mock...UV. != 0 & abs(R4$log2FC) > 3, ]

# Ignore R1 since we do not focus on that
DEGs <- unique(c(as.character(R2_DE$ID), as.character(R3_DE$ID), as.character(R4_DE$ID)))

f1 <- R1[R1$ID %in% DEGs, ]
f2 <- R2[R2$ID %in% DEGs, ]
f3 <- R3[R3$ID %in% DEGs, ]
f4 <- R4[R4$ID %in% DEGs, ]

d1 <- data.frame(ID = f1$ID, Gene = f1$gene.name, log2FC = f1$log2FC, DE = f1$DE..SCORE....1..Mock...Input..1..Mock...Input., TPM = f1$Presence.Absence..Input.)
d2 <- data.frame(ID = f2$ID, Gene = f2$gene.name, log2FC = f2$log2FC, DE = f2$DE..SCORE....1..Mock...WT..1..Mock...WT., TPM = f2$Presence.Absence..WT.)
d3 <- data.frame(ID = f3$ID, Gene = f3$gene.name, log2FC = f3$log2FC, DE = f3$DE..SCORE....1..Mock...dXCL1..1..Mock...dXCL1., TPM = f3$Presence.Absence..dXCL1.)
d4 <- data.frame(ID = f4$ID, Gene = f4$gene.name, log2FC = f4$log2FC, DE = f4$DE..SCORE....1..Mock...UV..1..Mock...UV., TPM = f4$Presence.Absence..UV.)

m1 <- merge(d1, d2, by = "ID", all = TRUE)
m2 <- merge(d3, d4, by = "ID", all = TRUE)
m3 <- merge(m1, m2, by = "ID", all = TRUE)

m4 <- data.frame(ID = m3$ID, Gene = m3$Gene.x.y, log2FC_R1 = m3$log2FC.x.x, log2FC_R2 = m3$log2FC.y.x, log2FC_R3 = m3$log2FC.x.y, log2FC_R4 = m3$log2FC.y.y, DE_R1 = m3$DE.x.x, DE_R2 = m3$DE.y.x, DE_R3 = m3$DE.x.y, DE_R4 = m3$DE.y.y, TPM_R1 = m3$TPM.x.x, TPM_R2 = m3$TPM.y.x, TPM_R3 = m3$TPM.x.y, TPM_R4 = m3$TPM.y.y)
m4[is.na(m4)] <- 0

names <- c("Not Cultivated", "RCMV-E wt", expression(paste("RCMV-E ", Delta, "vXCL1", sep = "")), "UV")

pheatmap(m4[3:6],
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         file = "pheatmap_degs.png",
         main = "Heatmap of filtered DEGs",
         scale = "none",
         show_rownames = FALSE,
         labels_col = names,
         labels_row = m4$Gene
         )

pheatmap(m4[4:6],
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         file = "pheatmap_degs_filtered_01.png",
         main = "Heatmap of filtered DEGs",
         scale = "none",
         show_rownames = FALSE,
         labels_col = names[-1],
         labels_row = m4$Gene
         )

pheatmap(m4[3],
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         file = "pheatmap_degs_filtered_02.png",
         main = "Heatmap of filtered DEGs",
         scale = "none",
         show_rownames = FALSE,
         labels_col = names[1],
         labels_row = m4$Gene
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

png(filename = "venn_diagram_downregulated.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
vennDiagram(v3, circle.col = c("red", "blue", "green", "grey"), main = "Downregulated Genes", names = names, cex = 1.1)
dev.off()

png(filename = "venn_diagram_upregulated.png", width = 30, height = 30, units = "cm", res = 600, pointsize = 20)
vennDiagram(v4, circle.col = c("red", "blue", "green", "grey"), main = "Upregulated Genes", names = names, cex = 1.1)
dev.off()

# Print overlap of all experiments
m4$Gene[as.integer(rownames(m4_down[m4_down$DE_R1 > 0 & m4_down$DE_R2 > 0 & m4_down$DE_R3 > 0 & m4_down$DE_R4 > 0, ]))]
m4$Gene[as.integer(rownames(m4_up[m4_up$DE_R1 > 0 & m4_up$DE_R2 > 0 & m4_up$DE_R3 > 0 & m4_up$DE_R4 > 0, ]))]

# Exporting

if (export == TRUE){
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
