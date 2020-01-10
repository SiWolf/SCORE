# --------------------------------------------
# Title: generate_benchmarking_plots.R
# Author: Silver A. Wolf
# Last Modified: Fr, 10.01.2020
# Version: 0.0.7
# --------------------------------------------

# Visualize DEG prediction accuracy
# Uses mostly standard R libraries

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

data_frame_accuracy = data.frame(data_1$X, data_1$ACC, data_2$ACC, data_3$ACC, data_4$ACC, data_5$ACC, data_6$ACC, data_7$ACC, data_8$ACC, data_9$ACC, data_10$ACC,
                                 data_11$ACC, data_12$ACC, data_13$ACC, data_14$ACC, data_15$ACC, data_16$ACC, data_17$ACC, data_18$ACC, data_19$ACC, data_20$ACC,
                                 data_21$ACC, data_22$ACC, data_23$ACC, data_24$ACC, data_25$ACC, data_26$ACC, data_27$ACC, data_28$ACC, data_29$ACC, data_30$ACC,
                                 data_31$ACC, data_32$ACC, data_33$ACC, data_34$ACC, data_35$ACC, data_36$ACC, data_37$ACC, data_38$ACC, data_39$ACC, data_40$ACC,
                                 data_41$ACC, data_42$ACC, data_43$ACC, data_44$ACC, data_45$ACC, data_46$ACC, data_47$ACC, data_48$ACC, data_49$ACC, data_50$ACC,
                                 data_51$ACC, data_52$ACC, data_53$ACC, data_54$ACC, data_55$ACC, data_56$ACC, data_57$ACC, data_58$ACC, data_59$ACC, data_60$ACC,
                                 data_61$ACC, data_62$ACC, data_63$ACC, data_64$ACC, data_65$ACC, data_66$ACC, data_67$ACC, data_68$ACC, data_69$ACC, data_70$ACC,
                                 data_71$ACC, data_72$ACC, data_73$ACC, data_74$ACC, data_75$ACC, data_76$ACC, data_77$ACC, data_78$ACC, data_79$ACC, data_80$ACC,
                                 data_81$ACC, data_82$ACC, data_83$ACC, data_84$ACC, data_85$ACC, data_86$ACC, data_87$ACC, data_88$ACC, data_89$ACC, data_90$ACC,
                                 data_91$ACC, data_92$ACC, data_93$ACC, data_94$ACC, data_95$ACC, data_96$ACC, data_97$ACC, data_98$ACC, data_99$ACC, data_100$ACC)
data_frame_degs = data.frame(data_1$X, data_1$DEGs, data_2$DEGs, data_3$DEGs, data_4$DEGs, data_5$DEGs, data_6$DEGs, data_7$DEGs, data_8$DEGs, data_9$DEGs, data_10$DEGs,
                             data_11$DEGs, data_12$DEGs, data_13$DEGs, data_14$DEGs, data_15$DEGs, data_16$DEGs, data_17$DEGs, data_18$DEGs, data_19$DEGs, data_20$DEGs,
                             data_21$DEGs, data_22$DEGs, data_23$DEGs, data_24$DEGs, data_25$DEGs, data_26$DEGs, data_27$DEGs, data_28$DEGs, data_29$DEGs, data_30$DEGs,
                             data_31$DEGs, data_32$DEGs, data_33$DEGs, data_34$DEGs, data_35$DEGs, data_36$DEGs, data_37$DEGs, data_38$DEGs, data_39$DEGs, data_40$DEGs,
                             data_41$DEGs, data_42$DEGs, data_43$DEGs, data_44$DEGs, data_45$DEGs, data_46$DEGs, data_47$DEGs, data_48$DEGs, data_49$DEGs, data_50$DEGs,
                             data_51$DEGs, data_52$DEGs, data_53$DEGs, data_54$DEGs, data_55$DEGs, data_56$DEGs, data_57$DEGs, data_58$DEGs, data_59$DEGs, data_60$DEGs,
                             data_61$DEGs, data_62$DEGs, data_63$DEGs, data_64$DEGs, data_65$DEGs, data_66$DEGs, data_67$DEGs, data_68$DEGs, data_69$DEGs, data_70$DEGs,
                             data_71$DEGs, data_72$DEGs, data_73$DEGs, data_74$DEGs, data_75$DEGs, data_76$DEGs, data_77$DEGs, data_78$DEGs, data_79$DEGs, data_80$DEGs,
                             data_81$DEGs, data_82$DEGs, data_83$DEGs, data_84$DEGs, data_85$DEGs, data_86$DEGs, data_87$DEGs, data_88$DEGs, data_89$DEGs, data_90$DEGs,
                             data_91$DEGs, data_92$DEGs, data_93$DEGs, data_94$DEGs, data_95$DEGs, data_96$DEGs, data_97$DEGs, data_98$DEGs, data_99$DEGs, data_100$DEGs)
data_frame_fdr = data.frame(data_1$X, data_1$FDR, data_2$FDR, data_3$FDR, data_4$FDR, data_5$FDR, data_6$FDR, data_7$FDR, data_8$FDR, data_9$FDR, data_10$FDR,
                            data_11$FDR, data_12$FDR, data_13$FDR, data_14$FDR, data_15$FDR, data_16$FDR, data_17$FDR, data_18$FDR, data_19$FDR, data_20$FDR,
                            data_21$FDR, data_22$FDR, data_23$FDR, data_24$FDR, data_25$FDR, data_26$FDR, data_27$FDR, data_28$FDR, data_29$FDR, data_30$FDR,
                            data_31$FDR, data_32$FDR, data_33$FDR, data_34$FDR, data_35$FDR, data_36$FDR, data_37$FDR, data_38$FDR, data_39$FDR, data_40$FDR,
                            data_41$FDR, data_42$FDR, data_43$FDR, data_44$FDR, data_45$FDR, data_46$FDR, data_47$FDR, data_48$FDR, data_49$FDR, data_50$FDR,
                            data_51$FDR, data_52$FDR, data_53$FDR, data_54$FDR, data_55$FDR, data_56$FDR, data_57$FDR, data_58$FDR, data_59$FDR, data_60$FDR,
                            data_61$FDR, data_62$FDR, data_63$FDR, data_64$FDR, data_65$FDR, data_66$FDR, data_67$FDR, data_68$FDR, data_69$FDR, data_70$FDR,
                            data_71$FDR, data_72$FDR, data_73$FDR, data_74$FDR, data_75$FDR, data_76$FDR, data_77$FDR, data_78$FDR, data_79$FDR, data_80$FDR,
                            data_81$FDR, data_82$FDR, data_83$FDR, data_84$FDR, data_85$FDR, data_86$FDR, data_87$FDR, data_88$FDR, data_89$FDR, data_90$FDR,
                            data_91$FDR, data_92$FDR, data_93$FDR, data_94$FDR, data_95$FDR, data_96$FDR, data_97$FDR, data_98$FDR, data_99$FDR, data_100$FDR)
data_frame_fnr = data.frame(data_1$X, data_1$FNR, data_2$FNR, data_3$FNR, data_4$FNR, data_5$FNR, data_6$FNR, data_7$FNR, data_8$FNR, data_9$FNR, data_10$FNR,
                            data_11$FNR, data_12$FNR, data_13$FNR, data_14$FNR, data_15$FNR, data_16$FNR, data_17$FNR, data_18$FNR, data_19$FNR, data_20$FNR,
                            data_21$FNR, data_22$FNR, data_23$FNR, data_24$FNR, data_25$FNR, data_26$FNR, data_27$FNR, data_28$FNR, data_29$FNR, data_30$FNR,
                            data_31$FNR, data_32$FNR, data_33$FNR, data_34$FNR, data_35$FNR, data_36$FNR, data_37$FNR, data_38$FNR, data_39$FNR, data_40$FNR,
                            data_41$FNR, data_42$FNR, data_43$FNR, data_44$FNR, data_45$FNR, data_46$FNR, data_47$FNR, data_48$FNR, data_49$FNR, data_50$FNR,
                            data_51$FNR, data_52$FNR, data_53$FNR, data_54$FNR, data_55$FNR, data_56$FNR, data_57$FNR, data_58$FNR, data_59$FNR, data_60$FNR,
                            data_61$FNR, data_62$FNR, data_63$FNR, data_64$FNR, data_65$FNR, data_66$FNR, data_67$FNR, data_68$FNR, data_69$FNR, data_70$FNR,
                            data_71$FNR, data_72$FNR, data_73$FNR, data_74$FNR, data_75$FNR, data_76$FNR, data_77$FNR, data_78$FNR, data_79$FNR, data_80$FNR,
                            data_81$FNR, data_82$FNR, data_83$FNR, data_84$FNR, data_85$FNR, data_86$FNR, data_87$FNR, data_88$FNR, data_89$FNR, data_90$FNR,
                            data_91$FNR, data_92$FNR, data_93$FNR, data_94$FNR, data_95$FNR, data_96$FNR, data_97$FNR, data_98$FNR, data_99$FNR, data_100$FNR)
data_frame_fns = data.frame(data_1$X, data_1$FN, data_2$FN, data_3$FN, data_4$FN, data_5$FN, data_6$FN, data_7$FN, data_8$FN, data_9$FN, data_10$FN,
                            data_11$FN, data_12$FN, data_13$FN, data_14$FN, data_15$FN, data_16$FN, data_17$FN, data_18$FN, data_19$FN, data_20$FN,
                            data_21$FN, data_22$FN, data_23$FN, data_24$FN, data_25$FN, data_26$FN, data_27$FN, data_28$FN, data_29$FN, data_30$FN,
                            data_31$FN, data_32$FN, data_33$FN, data_34$FN, data_35$FN, data_36$FN, data_37$FN, data_38$FN, data_39$FN, data_40$FN,
                            data_41$FN, data_42$FN, data_43$FN, data_44$FN, data_45$FN, data_46$FN, data_47$FN, data_48$FN, data_49$FN, data_50$FN,
                            data_51$FN, data_52$FN, data_53$FN, data_54$FN, data_55$FN, data_56$FN, data_57$FN, data_58$FN, data_59$FN, data_60$FN,
                            data_61$FN, data_62$FN, data_63$FN, data_64$FN, data_65$FN, data_66$FN, data_67$FN, data_68$FN, data_69$FN, data_70$FN,
                            data_71$FN, data_72$FN, data_73$FN, data_74$FN, data_75$FN, data_76$FN, data_77$FN, data_78$FN, data_79$FN, data_80$FN,
                            data_81$FN, data_82$FN, data_83$FN, data_84$FN, data_85$FN, data_86$FN, data_87$FN, data_88$FN, data_89$FN, data_90$FN,
                            data_91$FN, data_92$FN, data_93$FN, data_94$FN, data_95$FN, data_96$FN, data_97$FN, data_98$FN, data_99$FN, data_100$FN)
data_frame_fpr = data.frame(data_1$X, data_1$FPR, data_2$FPR, data_3$FPR, data_4$FPR, data_5$FPR, data_6$FPR, data_7$FPR, data_8$FPR, data_9$FPR, data_10$FPR,
                            data_11$FPR, data_12$FPR, data_13$FPR, data_14$FPR, data_15$FPR, data_16$FPR, data_17$FPR, data_18$FPR, data_19$FPR, data_20$FPR,
                            data_21$FPR, data_22$FPR, data_23$FPR, data_24$FPR, data_25$FPR, data_26$FPR, data_27$FPR, data_28$FPR, data_29$FPR, data_30$FPR,
                            data_31$FPR, data_32$FPR, data_33$FPR, data_34$FPR, data_35$FPR, data_36$FPR, data_37$FPR, data_38$FPR, data_39$FPR, data_40$FPR,
                            data_41$FPR, data_42$FPR, data_43$FPR, data_44$FPR, data_45$FPR, data_46$FPR, data_47$FPR, data_48$FPR, data_49$FPR, data_50$FPR,
                            data_51$FPR, data_52$FPR, data_53$FPR, data_54$FPR, data_55$FPR, data_56$FPR, data_57$FPR, data_58$FPR, data_59$FPR, data_60$FPR,
                            data_61$FPR, data_62$FPR, data_63$FPR, data_64$FPR, data_65$FPR, data_66$FPR, data_67$FPR, data_68$FPR, data_69$FPR, data_70$FPR,
                            data_71$FPR, data_72$FPR, data_73$FPR, data_74$FPR, data_75$FPR, data_76$FPR, data_77$FPR, data_78$FPR, data_79$FPR, data_80$FPR,
                            data_81$FPR, data_82$FPR, data_83$FPR, data_84$FPR, data_85$FPR, data_86$FPR, data_87$FPR, data_88$FPR, data_89$FPR, data_90$FPR,
                            data_91$FPR, data_92$FPR, data_93$FPR, data_94$FPR, data_95$FPR, data_96$FPR, data_97$FPR, data_98$FPR, data_99$FPR, data_100$FPR)
data_frame_fps = data.frame(data_1$X, data_1$FP, data_2$FP, data_3$FP, data_4$FP, data_5$FP, data_6$FP, data_7$FP, data_8$FP, data_9$FP, data_10$FP,
                            data_11$FP, data_12$FP, data_13$FP, data_14$FP, data_15$FP, data_16$FP, data_17$FP, data_18$FP, data_19$FP, data_20$FP,
                            data_21$FP, data_22$FP, data_23$FP, data_24$FP, data_25$FP, data_26$FP, data_27$FP, data_28$FP, data_29$FP, data_30$FP,
                            data_31$FP, data_32$FP, data_33$FP, data_34$FP, data_35$FP, data_36$FP, data_37$FP, data_38$FP, data_39$FP, data_40$FP,
                            data_41$FP, data_42$FP, data_43$FP, data_44$FP, data_45$FP, data_46$FP, data_47$FP, data_48$FP, data_49$FP, data_50$FP,
                            data_51$FP, data_52$FP, data_53$FP, data_54$FP, data_55$FP, data_56$FP, data_57$FP, data_58$FP, data_59$FP, data_60$FP,
                            data_61$FP, data_62$FP, data_63$FP, data_64$FP, data_65$FP, data_66$FP, data_67$FP, data_68$FP, data_69$FP, data_70$FP,
                            data_71$FP, data_72$FP, data_73$FP, data_74$FP, data_75$FP, data_76$FP, data_77$FP, data_78$FP, data_79$FP, data_80$FP,
                            data_81$FP, data_82$FP, data_83$FP, data_84$FP, data_85$FP, data_86$FP, data_87$FP, data_88$FP, data_89$FP, data_90$FP,
                            data_91$FP, data_92$FP, data_93$FP, data_94$FP, data_95$FP, data_96$FP, data_97$FP, data_98$FP, data_99$FP, data_100$FP)
data_frame_pre = data.frame(data_1$X, data_1$PRE, data_2$PRE, data_3$PRE, data_4$PRE, data_5$PRE, data_6$PRE, data_7$PRE, data_8$PRE, data_9$PRE, data_10$PRE,
                            data_11$PRE, data_12$PRE, data_13$PRE, data_14$PRE, data_15$PRE, data_16$PRE, data_17$PRE, data_18$PRE, data_19$PRE, data_20$PRE,
                            data_21$PRE, data_22$PRE, data_23$PRE, data_24$PRE, data_25$PRE, data_26$PRE, data_27$PRE, data_28$PRE, data_29$PRE, data_30$PRE,
                            data_31$PRE, data_32$PRE, data_33$PRE, data_34$PRE, data_35$PRE, data_36$PRE, data_37$PRE, data_38$PRE, data_39$PRE, data_40$PRE,
                            data_41$PRE, data_42$PRE, data_43$PRE, data_44$PRE, data_45$PRE, data_46$PRE, data_47$PRE, data_48$PRE, data_49$PRE, data_50$PRE,
                            data_51$PRE, data_52$PRE, data_53$PRE, data_54$PRE, data_55$PRE, data_56$PRE, data_57$PRE, data_58$PRE, data_59$PRE, data_60$PRE,
                            data_61$PRE, data_62$PRE, data_63$PRE, data_64$PRE, data_65$PRE, data_66$PRE, data_67$PRE, data_68$PRE, data_69$PRE, data_70$PRE,
                            data_71$PRE, data_72$PRE, data_73$PRE, data_74$PRE, data_75$PRE, data_76$PRE, data_77$PRE, data_78$PRE, data_79$PRE, data_80$PRE,
                            data_81$PRE, data_82$PRE, data_83$PRE, data_84$PRE, data_85$PRE, data_86$PRE, data_87$PRE, data_88$PRE, data_89$PRE, data_90$PRE,
                            data_91$PRE, data_92$PRE, data_93$PRE, data_94$PRE, data_95$PRE, data_96$PRE, data_97$PRE, data_98$PRE, data_99$PRE, data_100$PRE)
data_frame_tnr = data.frame(data_1$X, data_1$TNR, data_2$TNR, data_3$TNR, data_4$TNR, data_5$TNR, data_6$TNR, data_7$TNR, data_8$TNR, data_9$TNR, data_10$TNR,
                            data_11$TNR, data_12$TNR, data_13$TNR, data_14$TNR, data_15$TNR, data_16$TNR, data_17$TNR, data_18$TNR, data_19$TNR, data_20$TNR,
                            data_21$TNR, data_22$TNR, data_23$TNR, data_24$TNR, data_25$TNR, data_26$TNR, data_27$TNR, data_28$TNR, data_29$TNR, data_30$TNR,
                            data_31$TNR, data_32$TNR, data_33$TNR, data_34$TNR, data_35$TNR, data_36$TNR, data_37$TNR, data_38$TNR, data_39$TNR, data_40$TNR,
                            data_41$TNR, data_42$TNR, data_43$TNR, data_44$TNR, data_45$TNR, data_46$TNR, data_47$TNR, data_48$TNR, data_49$TNR, data_50$TNR,
                            data_51$TNR, data_52$TNR, data_53$TNR, data_54$TNR, data_55$TNR, data_56$TNR, data_57$TNR, data_58$TNR, data_59$TNR, data_60$TNR,
                            data_61$TNR, data_62$TNR, data_63$TNR, data_64$TNR, data_65$TNR, data_66$TNR, data_67$TNR, data_68$TNR, data_69$TNR, data_70$TNR,
                            data_71$TNR, data_72$TNR, data_73$TNR, data_74$TNR, data_75$TNR, data_76$TNR, data_77$TNR, data_78$TNR, data_79$TNR, data_80$TNR,
                            data_81$TNR, data_82$TNR, data_83$TNR, data_84$TNR, data_85$TNR, data_86$TNR, data_87$TNR, data_88$TNR, data_89$TNR, data_90$TNR,
                            data_91$TNR, data_92$TNR, data_93$TNR, data_94$TNR, data_95$TNR, data_96$TNR, data_97$TNR, data_98$TNR, data_99$TNR, data_100$TNR)
data_frame_tns = data.frame(data_1$X, data_1$TN, data_2$TN, data_3$TN, data_4$TN, data_5$TN, data_6$TN, data_7$TN, data_8$TN, data_9$TN, data_10$TN,
                            data_11$TN, data_12$TN, data_13$TN, data_14$TN, data_15$TN, data_16$TN, data_17$TN, data_18$TN, data_19$TN, data_20$TN,
                            data_21$TN, data_22$TN, data_23$TN, data_24$TN, data_25$TN, data_26$TN, data_27$TN, data_28$TN, data_29$TN, data_30$TN,
                            data_31$TN, data_32$TN, data_33$TN, data_34$TN, data_35$TN, data_36$TN, data_37$TN, data_38$TN, data_39$TN, data_40$TN,
                            data_41$TN, data_42$TN, data_43$TN, data_44$TN, data_45$TN, data_46$TN, data_47$TN, data_48$TN, data_49$TN, data_50$TN,
                            data_51$TN, data_52$TN, data_53$TN, data_54$TN, data_55$TN, data_56$TN, data_57$TN, data_58$TN, data_59$TN, data_60$TN,
                            data_61$TN, data_62$TN, data_63$TN, data_64$TN, data_65$TN, data_66$TN, data_67$TN, data_68$TN, data_69$TN, data_70$TN,
                            data_71$TN, data_72$TN, data_73$TN, data_74$TN, data_75$TN, data_76$TN, data_77$TN, data_78$TN, data_79$TN, data_80$TN,
                            data_81$TN, data_82$TN, data_83$TN, data_84$TN, data_85$TN, data_86$TN, data_87$TN, data_88$TN, data_89$TN, data_90$TN,
                            data_91$TN, data_92$TN, data_93$TN, data_94$TN, data_95$TN, data_96$TN, data_97$TN, data_98$TN, data_99$TN, data_100$TN)
data_frame_tpr = data.frame(data_1$X, data_1$TPR, data_2$TPR, data_3$TPR, data_4$TPR, data_5$TPR, data_6$TPR, data_7$TPR, data_8$TPR, data_9$TPR, data_10$TPR,
                            data_11$TPR, data_12$TPR, data_13$TPR, data_14$TPR, data_15$TPR, data_16$TPR, data_17$TPR, data_18$TPR, data_19$TPR, data_20$TPR,
                            data_21$TPR, data_22$TPR, data_23$TPR, data_24$TPR, data_25$TPR, data_26$TPR, data_27$TPR, data_28$TPR, data_29$TPR, data_30$TPR,
                            data_31$TPR, data_32$TPR, data_33$TPR, data_34$TPR, data_35$TPR, data_36$TPR, data_37$TPR, data_38$TPR, data_39$TPR, data_40$TPR,
                            data_41$TPR, data_42$TPR, data_43$TPR, data_44$TPR, data_45$TPR, data_46$TPR, data_47$TPR, data_48$TPR, data_49$TPR, data_50$TPR,
                            data_51$TPR, data_52$TPR, data_53$TPR, data_54$TPR, data_55$TPR, data_56$TPR, data_57$TPR, data_58$TPR, data_59$TPR, data_60$TPR,
                            data_61$TPR, data_62$TPR, data_63$TPR, data_64$TPR, data_65$TPR, data_66$TPR, data_67$TPR, data_68$TPR, data_69$TPR, data_70$TPR,
                            data_71$TPR, data_72$TPR, data_73$TPR, data_74$TPR, data_75$TPR, data_76$TPR, data_77$TPR, data_78$TPR, data_79$TPR, data_80$TPR,
                            data_81$TPR, data_82$TPR, data_83$TPR, data_84$TPR, data_85$TPR, data_86$TPR, data_87$TPR, data_88$TPR, data_89$TPR, data_90$TPR,
                            data_91$TPR, data_92$TPR, data_93$TPR, data_94$TPR, data_95$TPR, data_96$TPR, data_97$TPR, data_98$TPR, data_99$TPR, data_100$TPR)
data_frame_tps = data.frame(data_1$X, data_1$TP, data_2$TP, data_3$TP, data_4$TP, data_5$TP, data_6$TP, data_7$TP, data_8$TP, data_9$TP, data_10$TP,
                            data_11$TP, data_12$TP, data_13$TP, data_14$TP, data_15$TP, data_16$TP, data_17$TP, data_18$TP, data_19$TP, data_20$TP,
                            data_21$TP, data_22$TP, data_23$TP, data_24$TP, data_25$TP, data_26$TP, data_27$TP, data_28$TP, data_29$TP, data_30$TP,
                            data_31$TP, data_32$TP, data_33$TP, data_34$TP, data_35$TP, data_36$TP, data_37$TP, data_38$TP, data_39$TP, data_40$TP,
                            data_41$TP, data_42$TP, data_43$TP, data_44$TP, data_45$TP, data_46$TP, data_47$TP, data_48$TP, data_49$TP, data_50$TP,
                            data_51$TP, data_52$TP, data_53$TP, data_54$TP, data_55$TP, data_56$TP, data_57$TP, data_58$TP, data_59$TP, data_60$TP,
                            data_61$TP, data_62$TP, data_63$TP, data_64$TP, data_65$TP, data_66$TP, data_67$TP, data_68$TP, data_69$TP, data_70$TP,
                            data_71$TP, data_72$TP, data_73$TP, data_74$TP, data_75$TP, data_76$TP, data_77$TP, data_78$TP, data_79$TP, data_80$TP,
                            data_81$TP, data_82$TP, data_83$TP, data_84$TP, data_85$TP, data_86$TP, data_87$TP, data_88$TP, data_89$TP, data_90$TP,
                            data_91$TP, data_92$TP, data_93$TP, data_94$TP, data_95$TP, data_96$TP, data_97$TP, data_98$TP, data_99$TP, data_100$TP)

# Visualize results of the first 8 simulations

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

# FNs
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

# FPs
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

# TNs
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

# TPs
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

# Visualize Averages
average_data_frame = data.frame(Tools = tools, Av_ACC = rowMeans(data_frame_accuracy[2:i]), Av_DEGs = rowMeans(data_frame_degs[2:i]), Av_FDR = rowMeans(data_frame_fdr[2:i]), Av_FNR = rowMeans(data_frame_fnr[2:i]), Av_FNs = rowMeans(data_frame_fns[2:i]), Av_FPR = rowMeans(data_frame_fpr[2:i]), Av_FPs = rowMeans(data_frame_fps[2:i]), Av_PRE = rowMeans(data_frame_pre[2:i]), Av_TNR = rowMeans(data_frame_tnr[2:i]), Av_TNs = rowMeans(data_frame_tns[2:i]), Av_TPR = rowMeans(data_frame_tpr[2:i]), Av_TPs = rowMeans(data_frame_tps[2:i]))

# Accuray
png(filename = "benchmarking_diagram_av_acc.png", width = 20, height = 20, units = "cm", res = 600)
plot(average_data_frame$Av_ACC, type = "o", lty = 1, pch = 15, col = "blue", axes = FALSE, xlab = "Tools", ylab = "Av_Acc", lwd = 2.5)
axis(1, 1:7, lab = tools)
axis(2)
box()
dev.off()

# Precision
png(filename = "benchmarking_diagram_av_pre.png", width = 20, height = 20, units = "cm", res = 600)
plot(average_data_frame$Av_PRE, type = "o", lty = 1, pch = 15, col = "red", axes = FALSE, xlab = "Tools", ylab = "Av_Pre", lwd = 2.5)
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