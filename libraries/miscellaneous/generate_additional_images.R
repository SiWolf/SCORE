# --------------------------------------------
# Title: generate_additional_images.R
# Author: Silver A. Wolf
# Last Modified: Thur, 26.09.2019
# Version: 0.1.3
# --------------------------------------------

# This script is used to generate additional images for publications, etc.
# It will need to be heavily adapted for any new data

# Libraries

library("ggplot2")
library("pheatmap")
library("reshape2")
library("viridis")

# Variables

export = FALSE

# Preprocessing

setwd("../../deg/V1/")

final_summary <- read.csv(file = "summary.tsv", header = TRUE, sep = "\t", quote = "")

for (i in 1:nrow(final_summary)){
  if (as.character(final_summary$TIGRFAM.main.role[i]) == ""){
    final_summary$TIGRFAM.main.role[i] <- "Unknown function"
  }
}

degs = final_summary[final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. != 0, ]
tpm_values <- read.csv(file = "filtered_gene_counts_tpm.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)
tpm_values$TIGRFAM <- final_summary$TIGRFAM.main.role

if (export == True){
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
  
  dev.off()
}