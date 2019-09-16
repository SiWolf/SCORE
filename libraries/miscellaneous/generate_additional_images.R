# --------------------------------------------
# Title: generate_additional_images.R
# Author: Silver A. Wolf
# Last Modified: Mo, 16.09.2019
# Version: 0.0.9
# --------------------------------------------

# Libraries

library("ggplot2")
library("pheatmap")
library("reshape2")
#library("gplots")
#library("RColorBrewer")

# Preprocessing

setwd("../deg/")

final_summary <- read.csv(file = "summary.tsv", header = TRUE, sep = "\t", quote = "")

for (i in 1:nrow(final_summary)){
  if (as.character(final_summary$TIGRFAM.main.role[i]) == ""){
    final_summary$TIGRFAM.main.role[i] <- "Unknown function"
  }
}

degs = final_summary[final_summary$DE..SCORE....1..Wildtype...ZnGroup..1..Wildtype...ZnGroup. != 0, ]
tpm_values <- read.csv(file = "filtered_gene_counts_tpm.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)
tpm_values$TIGRFAM <- final_summary$TIGRFAM.main.role

pdf("deg_analysis_graphs_extended.pdf", paper = "a4")

# Boxplots

g1 <- ggplot(data = final_summary, aes(x = TIGRFAM.main.role, y = log2FC)) + 
  geom_boxplot(fill = "lightgrey") +
  stat_boxplot(geom = "errorbar") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Boxplot of the log2FC of all genes within each TIGRFAM Main Role (Baseline)") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC")

g2 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = log2FC)) + 
  geom_boxplot(fill = "lightgrey") +
  stat_boxplot(geom = "errorbar") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Boxplot of the log2FC of DEGs within each TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC")

# Heatmaps

high_FC <- degs[abs(degs$log2FC) > 2, ]
tpm_degs <- tpm_values[tpm_values$ID %in% high_FC$ID, ]
tpm_degs_ordered <- tpm_degs[order(tpm_degs$TIGRFAM, tpm_degs$gene.name), ]
tpm_degs_ordered$gene.name <- replace(tpm_degs_ordered$gene.name, tpm_degs_ordered$gene.name == "", "hypothetical_gene")
rownames(tpm_degs_ordered) <- tpm_degs_ordered$gene.name

groups = data.frame(Groups = c(rep("Control Group", 3), rep("Zn Group", 3)))
rownames(groups) <- c("X6122_C1", "X6122_C5", "X6122_C6", "X6122_Z1", "X6122_Z5", "X6122_Z6")

g3 <- pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
        annotation_col = groups, fontsize_row = 9, fontsize_col = 9,
        cellheight = 10, cellwidth = 20,
        border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
        main = "Heatmap of genes (TPM)\nwith high log2FC (> 2)\nbetween sample groups\n",
        scale = "none")

g4 <- pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
         annotation_col = groups, fontsize_row = 9, fontsize_col = 9,
         cellheight = 10, cellwidth = 20,
         border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
         main = "Heatmap of genes (TPM)\nwith high log2FC (> 2)\nbetween sample groups (scaled)\n",
         scale = "row")

# Exporting Pheatmaps

#pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
         #annotation_col = groups, fontsize_row = 9, fontsize_col = 9,
         #cellheight = 10, cellwidth = 20,
         #border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
         #main = "Heatmap of genes (TPM)\nwith high log2FC (> 2)\nbetween sample groups\n",
         #scale = "none", file = "new_pheatmap_01.png")

#pheatmap(tpm_degs_ordered[4:9], cluster_cols = TRUE, cluster_rows = FALSE, annotation_row = tpm_degs_ordered[10],
         #annotation_col = groups, fontsize_row = 9, fontsize_col = 9,
         #cellheight = 10, cellwidth = 20,
         #border_color = "black", gaps_row = c(5, 9, 10, 12, 14, 18, 20, 21, 26),
         #main = "Heatmap of genes (TPM)\nwith high log2FC (> 2)\nbetween sample groups (scaled)\n",
         #scale = "row", file = "new_pheatmap_02.png")

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

g21 <- ggplot(final_summary, aes(y = log2FC, x = TIGRFAM.main.role)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the log2FC (all genes) per TIGRFAM Main Role (Baseline)") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC") +
  labs(fill = "Groups")

g22 <- ggplot(degs, aes(y = log2FC, x = TIGRFAM.main.role)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the log2FC (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC") +
  labs(fill = "Groups")

# Old Plots

g13 <- ggplot(data = BFtranspose, aes(x = ID, y = TPM, colour = Group)) +
  geom_bar(stat = "identity") +
  theme_void()

g14 <- ggplot(data = BFtranspose, aes(x = TIGRFAM_Main_Role, y = TPM, colour = Group)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the total TPM per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("Total TPM")

g15 <- ggplot(data = BFtranspose, aes(x = factor(ID), y = TPM, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM per Gene") +
  xlab("Genes") +
  ylab("Total TPM")

g16 <- ggplot(data = BFtranspose, aes(x = factor(TIGRFAM_Main_Role), y = TPM, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the total TPM per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("Total TPM")

g17 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = Mean.TPM..Wildtype.)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM (WT) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("TPM")

g18 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = Mean.TPM..ZnGroup.)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the TPM (ZnGroup) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("TPM")

g19 <- ggplot(data = final_summary, aes(x = TIGRFAM.main.role, y = log2FC)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the cumulative log2FC (all genes) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC (all genes)")

g20 <- ggplot(data = degs, aes(x = TIGRFAM.main.role, y = log2FC)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Barchart of the cumulative log2FC (DEGs) per TIGRFAM Main Role") +
  xlab("TIGRFAM Main Roles") +
  ylab("log2FC (DEGs)")

print(g1)
print(g2)
#print(g3)
#print(g4)
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
      
dev.off()