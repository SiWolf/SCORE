# --------------------------------------------
# Title: compare_experiments_visualization.R
# Author: Silver A. Wolf
# Last Modified: Mo, 05.08.2019
# Version: 0.0.1
# --------------------------------------------

library("limma")

data_1 <- read.csv(file = "V1/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")
data_2 <- read.csv(file = "V2/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")
data_3 <- read.csv(file = "V3/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")
data_4 <- read.csv(file = "V4/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")
data_5 <- read.csv(file = "V5/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")
data_6 <- read.csv(file = "V6/consensus_diffexpr_results_extended.csv", header = TRUE, sep = ",")

# V1 vs. V2

data_1_genes <- data_1$gene.name[data_1$gene.name != ""]
data_2_genes <- data_2$gene.name[data_2$gene.name != ""]

merged_genes = unlist(list())
data_1_vector = unlist(list())
data_2_vector = unlist(list())

i = 0
for (gene in data_1_genes){
  merged_genes[i] = gene
  if (gene %in% data_2_genes){
    data_1_vector[i] = 1
    data_2_vector[i] = 1
  }
  else{
    data_1_vector[i] = 1
    data_2_vector[i] = 0
  }
  i = i + 1
}
for (gene in data_2_genes){
  if (gene %in% merged_genes){
    i = i
  }
  else{
    merged_genes[i] = gene
    data_1_vector[i] = 0
    data_2_vector[i] = 1
    i = i + 1
  }
}

data_frame_1 <- data.frame(genes = merged_genes, d1 = data_1_vector, d2 = data_2_vector)
write.csv(data_frame_1, file = "c1.csv")
data_frame_1_edit <- data.frame(d1 = data_1_vector, d2 = data_2_vector)

png(filename = "c1.png", width = 20, height = 20, units = "cm", res = 600)
vennDiagram(data_frame_1_edit, circle.col = c("red", "blue"), names = c("V1", "V2"))
dev.off()

# V3 vs. V4

data_3_genes <- data_3$gene.name[data_3$gene.name != ""]
data_4_genes <- data_4$gene.name[data_4$gene.name != ""]

merged_genes = unlist(list())
data_3_vector = unlist(list())
data_4_vector = unlist(list())

i = 0
for (gene in data_3_genes){
  merged_genes[i] = gene
  if (gene %in% data_4_genes){
    data_3_vector[i] = 1
    data_4_vector[i] = 1
  }
  else{
    data_3_vector[i] = 1
    data_4_vector[i] = 0
  }
  i = i + 1
}
for (gene in data_4_genes){
  if (gene %in% merged_genes){
    i = i
  }
  else{
    merged_genes[i] = gene
    data_3_vector[i] = 0
    data_4_vector[i] = 1
    i = i + 1
  }
}

data_frame_2 <- data.frame(genes = merged_genes, d3 = data_3_vector, d4 = data_4_vector)
write.csv(data_frame_2, file = "c2.csv")
data_frame_2_edit <- data.frame(d3 = data_3_vector, d4 = data_4_vector)

png(filename = "c2.png", width = 20, height = 20, units = "cm", res = 600)
vennDiagram(data_frame_2_edit, circle.col = c("red", "blue"), names = c("V3", "V4"))
dev.off()