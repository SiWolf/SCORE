# --------------------------------------------
# Title: generate_benchmarking_plots.R
# Author: Silver A. Wolf
# Last Modified: Thue, 07.01.2020
# Version: 0.0.5
# --------------------------------------------

# Visualize DEG prediction accuracy
# Uses mostly standard R libraries
# TO-DO: Visualize averages?

#library("PRROC")

# Data preprocessing
setwd("deg/")
folders <- list.files(path = ".", pattern = "B[0-9]*")
i = 1
tools <- c("baySeq", "DESeq2", "edgeR", "limma", "NOISeq", "sleuth", "SCORE")

for (f in folders) {
  assign(paste("data_", i, sep = ""), read.csv(file = paste(f, "/deg_summary.csv", sep = ""), header = TRUE, sep = ","))
  i = i + 1
}

data_frame_accuracy = data.frame(data_1$X, data_1$ACC, data_2$ACC, data_3$ACC, data_4$ACC, data_5$ACC, data_6$ACC, data_7$ACC, data_8$ACC, data_9$ACC, data_10$ACC)
data_frame_degs = data.frame(data_1$X, data_1$DEGs, data_2$DEGs, data_3$DEGs, data_4$DEGs, data_5$DEGs, data_6$DEGs, data_7$DEGs, data_8$DEGs, data_9$DEGs, data_10$DEGs)
data_frame_fns = data.frame(data_1$X, data_1$FN, data_2$FN, data_3$FN, data_4$FN, data_5$FN, data_6$FN, data_7$FN, data_8$FN, data_9$FN, data_10$FN)
data_frame_fps = data.frame(data_1$X, data_1$FP, data_2$FP, data_3$FP, data_4$FP, data_5$FP, data_6$FP, data_7$FP, data_8$FP, data_9$FP, data_10$FP)
data_frame_tns = data.frame(data_1$X, data_1$TN, data_2$TN, data_3$TN, data_4$TN, data_5$TN, data_6$TN, data_7$TN, data_8$TN, data_9$TN, data_10$TN)
data_frame_tps = data.frame(data_1$X, data_1$TP, data_2$TP, data_3$TP, data_4$TP, data_5$TP, data_6$TP, data_7$TP, data_8$TP, data_9$TP, data_10$TP)
data_frame_tpr = data.frame(data_1$X, data_1$TPR, data_2$TPR, data_3$TPR, data_4$TPR, data_5$TPR, data_6$TPR, data_7$TPR, data_8$TPR, data_9$TPR, data_10$TPR)
data_frame_tnr = data.frame(data_1$X, data_1$TNR, data_2$TNR, data_3$TNR, data_4$TNR, data_5$TNR, data_6$TNR, data_7$TNR, data_8$TNR, data_9$TNR, data_10$TNR)
data_frame_fdr = data.frame(data_1$X, data_1$FDR, data_2$FDR, data_3$FDR, data_4$FDR, data_5$FDR, data_6$FDR, data_7$FDR, data_8$FDR, data_9$FDR, data_10$FDR)
data_frame_fpr = data.frame(data_1$X, data_1$FPR, data_2$FPR, data_3$FPR, data_4$FPR, data_5$FPR, data_6$FPR, data_7$FPR, data_8$FPR, data_9$FPR, data_10$FPR)
data_frame_fnr = data.frame(data_1$X, data_1$FNR, data_2$FNR, data_3$FNR, data_4$FNR, data_5$FNR, data_6$FNR, data_7$FNR, data_8$FNR, data_9$FNR, data_10$FNR)
data_frame_pre = data.frame(data_1$X, data_1$PRE, data_2$PRE, data_3$PRE, data_4$PRE, data_5$PRE, data_6$PRE, data_7$PRE, data_8$PRE, data_9$PRE, data_10$PRE)

# First 8 Simulations

# Accuracy
png(filename = "benchmarking_diagram_acc.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$ACC, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0.95, 1), axes = FALSE, xlab = "Tools", ylab = "Accuracy", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$ACC, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$ACC, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$ACC, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$ACC, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$ACC, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$ACC, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$ACC, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.98, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.9, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# DEGs
png(filename = "benchmarking_diagram_degs.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$DEGs, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(50, 800), axes = FALSE, xlab = "Tools", ylab = "Amount of identified DEGs", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$DEGs, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$DEGs, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$DEGs, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$DEGs, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$DEGs, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$DEGs, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$DEGs, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1.5, 800, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.7, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 4)
dev.off()

# FN
png(filename = "benchmarking_diagram_fns.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$FN, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 4), axes = FALSE, xlab = "Tools", ylab = "Amount of FNs", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$FN, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$FN, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$FN, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$FN, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$FN, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$FN, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$FN, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 4, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# FP
png(filename = "benchmarking_diagram_fps.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$FP, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 300), axes = FALSE, xlab = "Tools", ylab = "Amount of FPs", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$FP, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$FP, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$FP, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$FP, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$FP, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$FP, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$FP, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 300, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# TN
png(filename = "benchmarking_diagram_tns.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$TN, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(4000, 4700), axes = FALSE, xlab = "Tools", ylab = "Amount of TNs", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$TN, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$TN, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$TN, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$TN, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$TN, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$TN, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$TN, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1.5, 4175, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 4)
dev.off()

# TP
png(filename = "benchmarking_diagram_tps.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$TP, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(50, 550), axes = FALSE, xlab = "Tools", ylab = "Amount of TPs", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$TP, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$TP, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$TP, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$TP, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$TP, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$TP, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$TP, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 485, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# TNR
png(filename = "benchmarking_diagram_tnr.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$TNR, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0.94, 1), axes = FALSE, xlab = "Tools", ylab = "TNR", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$TNR, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$TNR, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$TNR, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$TNR, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$TNR, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$TNR, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$TNR, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.98, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# TPR
png(filename = "benchmarking_diagram_tpr.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$TPR, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0.99, 1), axes = FALSE, xlab = "Tools", ylab = "TPR", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$TPR, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$TPR, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$TPR, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$TPR, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$TPR, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$TPR, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$TPR, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.999, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# FDR
png(filename = "benchmarking_diagram_fdr.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$FDR, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 0.4), axes = FALSE, xlab = "Tools", ylab = "FDR", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$FDR, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$FDR, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$FDR, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$FDR, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$FDR, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$FDR, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$FDR, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.4, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# FPR
png(filename = "benchmarking_diagram_fpr.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$FPR, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 0.1), axes = FALSE, xlab = "Tools", ylab = "FPR", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$FPR, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$FPR, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$FPR, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$FPR, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$FPR, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$FPR, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$FPR, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.1, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# FNR
png(filename = "benchmarking_diagram_fnr.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$FNR, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0, 0.1), axes = FALSE, xlab = "Tools", ylab = "FNR", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$FNR, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$FNR, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$FNR, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$FNR, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$FNR, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$FNR, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$FNR, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.1, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# PRE
png(filename = "benchmarking_diagram_pre.png", width = 20, height = 20, units = "cm", res = 600)
plot(data_1$PRE, type = "o", lty = 1, pch = 15, col = "blue", ylim = c(0.6, 1), axes = FALSE, xlab = "Tools", ylab = "PRE", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
lines(data_2$PRE, type = "o", lty = 1, pch = 16, col = "red", lwd = 2.5)
lines(data_3$PRE, type = "o", lty = 1, pch = 17, col = "green", lwd = 2.5)
lines(data_4$PRE, type = "o", lty = 1, pch = 18, col = "orange", lwd = 2.5)
lines(data_5$PRE, type = "o", lty = 1, pch = 18, col = "grey", lwd = 2.5)
lines(data_6$PRE, type = "o", lty = 1, pch = 18, col = "black", lwd = 2.5)
lines(data_7$PRE, type = "o", lty = 1, pch = 18, col = "cyan", lwd = 2.5)
lines(data_8$PRE, type = "o", lty = 1, pch = 18, col = "yellow", lwd = 2.5)
box()
legend(1, 0.85, c("1", "2", "3", "4", "5", "6", "7", "8"), cex = 0.8, col = c("blue", "red", "green", "orange", "grey", "black", "cyan", "yellow"), pch = 15:18, lty = 1, ncol = 2)
dev.off()

# Average Accuracy
average_data_frame = data.frame(Tools = tools, Av_ACC = rowMeans(data_frame_accuracy[2:i]), Av_DEG = rowMeans(data_frame_degs[2:i]), Av_FDR = rowMeans(data_frame_fdr[2:i]), Av_FNR = rowMeans(data_frame_fnr[2:i]), Av_FNS = rowMeans(data_frame_fns[2:i]), Av_FPR = rowMeans(data_frame_fpr[2:i]), Av_PRE = rowMeans(data_frame_pre[2:i]), Av_TNR = rowMeans(data_frame_tnr[2:i]), Av_TNS = rowMeans(data_frame_tns[2:i]), Av_TPR = rowMeans(data_frame_tpr[2:i]), Av_TPS = rowMeans(data_frame_tps[2:i]))
png(filename = "benchmarking_diagram_av_acc.png", width = 20, height = 20, units = "cm", res = 600)
plot(average_data_frame$Av_ACC, type = "o", lty = 1, pch = 15, col = "blue", axes = FALSE, xlab = "Tools", ylab = "Av_Acc", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
box()
dev.off()

# ROC Curve
#bayseq_tpr <- unlist(data_frame_tpr[1, 2:9])
#bayseq_fpr <- unlist(data_frame_fpr[1, 2:9])
#noiseq_tpr <- unlist(data_frame_tpr[5, 2:9])
#noiseq_fpr <- unlist(data_frame_fpr[5, 2:9])

#bayseq_pr <- pr.curve(bayseq_fpr, bayseq_tpr, curve = TRUE)
#bayseq_roc <- roc.curve(bayseq_fpr, bayseq_tpr, curve = TRUE)

#noiseq_pr <- pr.curve(noiseq_fpr, noiseq_tpr, curve = TRUE)
#noiseq_roc <- roc.curve(noiseq_fpr, noiseq_tpr, curve = TRUE)

#plot(bayseq_pr, color = "red", auc.main = FALSE)
#plot(noiseq_pr, color = "green", add = TRUE)

#plot(bayseq_roc, color = "red", auc.main = FALSE)
#plot(noiseq_roc, color = "green", add = TRUE)