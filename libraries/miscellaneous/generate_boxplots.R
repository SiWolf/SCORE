# --------------------------------------------
# Title: generate_boxplots.R
# Author: Silver A. Wolf
# Last Modified: Wed, 04.12.2019
# Version: 0.0.1
# --------------------------------------------

# Visualize DEGs as Boxplots and other diagrams
# Using ggplot2

# Installers
#install.packages("ggplot2")

# Imports
library("ggplot2")

# Data preprocessing
data_1 <- read.csv(file = "boxplot_data/deg_summary_c1.csv", header = TRUE, sep = ",")
data_2 <- read.csv(file = "boxplot_data/deg_summary_c2.csv", header = TRUE, sep = ",")
data_3 <- read.csv(file = "boxplot_data/deg_summary_c3.csv", header = TRUE, sep = ",")
data_4 <- read.csv(file = "boxplot_data/deg_summary_c4.csv", header = TRUE, sep = ",")

data_5 <- read.csv(file = "boxplot_data/consensus_diffexpr_results_c1.csv", header = TRUE, sep = ",")
data_6 <- read.csv(file = "boxplot_data/consensus_diffexpr_results_c2.csv", header = TRUE, sep = ",")
data_7 <- read.csv(file = "boxplot_data/consensus_diffexpr_results_c3.csv", header = TRUE, sep = ",")
data_8 <- read.csv(file = "boxplot_data/consensus_diffexpr_results_c4.csv", header = TRUE, sep = ",")

score_row <- c(length(data_5$X), length(data_6$X),length(data_7$X),length(data_8$X))

new_data_frame <- data.frame(data_1$DEGS, data_2$DEGS, data_3$DEGS, data_4$DEGS)

final_df <- as.data.frame(t(new_data_frame))

final_final_df <- data.frame(final_df, score_row)

colnames(final_final_df) <- c("baySeq", "DESeq2", "edgeR", "limma", "NOISeq", "sleuth", "SCORE")

# GGplot boxplot
a = data.frame(Tools = "baySeq", value = final_final_df$baySeq)
b = data.frame(Tools = "DESeq2", value = final_final_df$DESeq2)
c = data.frame(Tools = "edgeR", value = final_final_df$edgeR)
d = data.frame(Tools = "limma", value = final_final_df$limma)
e = data.frame(Tools = "NOISeq", value = final_final_df$NOISeq)
f = data.frame(Tools = "sleuth", value = final_final_df$sleuth)
g = data.frame(Tools = "SCORE", value = final_final_df$SCORE)

plot.data = rbind(a, b, c, d, e, f, g)

png(filename = "boxplot_data/deg_analysis_box_diagram_01.png", width = 20, height = 20, units = "cm", res = 600)
ggplot(plot.data, aes(x=Tools,y=value,fill=Tools)) + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + scale_fill_manual(values=c("blue", "red", "green", "grey", "orange", "cyan", "sienna")) + ylab("Amount of DEGs identified") + xlab("Tools") + guides(fill = FALSE)
dev.off()

# Other graphs
png(filename = "boxplot_data/deg_analysis_diagram_02.png", width = 20, height = 20, units = "cm", res = 600)
plot(c(data_1$DEGS, score_row[1]), type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 2500), axes = FALSE, xlab = "Tools", ylab = "Amount of identified DEGs", lwd = 2.5)
axis(1, 1:7, lab = c("baySeq", "DESeq2", "edgeR", "limma", "NOISeq", "sleuth", "SCORE"))
axis(2)
lines(c(data_2$DEGS, score_row[2]), type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(c(data_3$DEGS, score_row[3]), type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(c(data_4$DEGS, score_row[4]), type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
box()
legend(1, 2500, c("C1", "C2", "C3", "C4"), cex = 1.1, col = c("blue", "red", "green", "orange"), pch = 15:18, lty = 1)
dev.off()