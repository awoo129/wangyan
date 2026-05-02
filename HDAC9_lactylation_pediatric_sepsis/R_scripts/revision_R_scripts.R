# =============================================================================
# BMC Bioinformatics Revision - Required R Analyses
# =============================================================================
# 本脚本包含修稿所需的3项新分析：
#   1. M3: K=3/K=4 聚类 + LARS连续风险评分
#   2. M5: GSEA FDR<0.05 严格阈值结果
#   3. m2: 单细胞CD14+ HDAC9梯度供体水平统计检验
# =============================================================================

library(ConsensusClusterPlus)
library(limma)
library(clusterProfiler)
library(survival)
library(survminer)
library(glmnet)
library(lme4)
library(ggplot2)

# =============================================================================
# 通用数据加载（请根据实际路径调整）
# =============================================================================
# 假设以下对象已加载：
# expr_train: GSE66099表达矩阵（行=基因，列=样本）
# meta_train: GSE66099临床数据，含 survival_time, event, HDAC9等
# sc_meta: GSE167363单细胞元数据，含 cell_type, donor_id, group, HDAC9_score
# gene_list: 预排序基因列表（log2FC * -log10(pvalue)）

# =============================================================================
# M3: K=3/K=4 聚类 + LARS连续风险评分
# =============================================================================

cat("\n===== M3: K=3/K=4 Consensus Clustering =====\n")

# --- K=3 和 K=4 聚类 ---
# 重新运行ConsensusClusterPlus，保存K=3和K=4结果
# 注意：需要原始表达矩阵和距离矩阵

# 假设 cc_result 是之前的ConsensusClusterPlus结果
# 如果需要重新运行：
# cc_result <- ConsensusClusterPlus(
#   as.matrix(expr_train),
#   maxK = 6,
#   reps = 1000,
#   pItem = 0.8,
#   pFeature = 1,
#   clusterAlg = "ward.D2",
#   distance = "pearson",
#   seed = 1234,
#   plot = "pdf",
#   writeTable = TRUE
# )

# 提取K=2, K=3, K=4的聚类结果
# cluster_k2 <- cc_result[[2]]$consensusClass
# cluster_k3 <- cc_result[[3]]$consensusClass
# cluster_k4 <- cc_result[[4]]$consensusClass

# --- K=3 生存分析 ---
# meta_train$cluster_k3 <- factor(cluster_k3)
# fit_k3 <- survfit(Surv(survival_time, event) ~ cluster_k3, data = meta_train)
# print(ggsurvplot(fit_k3, pval = TRUE, risk.table = TRUE,
#                  title = "K=3 Consensus Clustering Survival"))

# --- K=4 生存分析 ---
# meta_train$cluster_k4 <- factor(cluster_k4)
# fit_k4 <- survfit(Surv(survival_time, event) ~ cluster_k4, data = meta_train)
# print(ggsurvplot(fit_k4, pval = TRUE, risk.table = TRUE,
#                  title = "K=4 Consensus Clustering Survival"))

# --- 报告各聚类大小 ---
# cat("\nK=2 cluster sizes:\n")
# print(table(cluster_k2))
# cat("\nK=3 cluster sizes:\n")
# print(table(cluster_k3))
# cat("\nK=4 cluster sizes:\n")
# print(table(cluster_k4))

# --- LARS: LASSO-Cox 连续风险评分 ---
cat("\n===== M3: LARS (LASSO-Cox) Continuous Risk Score =====\n")

# 使用乳酸化相关基因构建LARS
# lars_genes <- c("HDAC9", "LDHA", "EP300", "PKM", "GAPDH", "PKM2", "ENO1")
# expr_lars <- t(expr_train[rownames(expr_train) %in% lars_genes, ])

# LASSO-Cox回归
# set.seed(42)
# x_lars <- as.matrix(expr_lars)
# y_lars <- Surv(meta_train$survival_time, meta_train$event)
#
# cv_fit <- cv.glmnet(x_lars, y_lars, family = "cox", alpha = 1, nfolds = 10)
# best_lambda <- cv_fit$lambda.min
# cat("Best lambda:", best_lambda, "\n")
#
# lars_model <- glmnet(x_lars, y_lars, family = "cox", alpha = 1, lambda = best_lambda)
# lars_coef <- coef(lars_model)
# cat("LARS coefficients:\n")
# print(lars_coef)

# 计算LARS评分
# lars_score <- as.numeric(x_lars %*% lars_coef)
# meta_train$LARS <- lars_score

# LARS C-index
# cindex_lars <- survConcordance(Surv(survival_time, event) ~ LARS, data = meta_train)
# cat("LARS C-index:", cindex_lars$concordance, "\n")
# cat("LARS C-index 95%CI:", cindex_lars$concordance - 1.96 * cindex_lars$std.err,
#     "-", cindex_lars$concordance + 1.96 * cindex_lars$std.err, "\n")

# 二分类聚类 C-index 对比
# cindex_binary <- survConcordance(Surv(survival_time, event) ~ as.numeric(cluster_k2)-1, data = meta_train)
# cat("Binary cluster C-index:", cindex_binary$concordance, "\n")

# LARS连续分组（三等分）的生存分析
# meta_train$LARS_tertile <- cut(meta_train$LARS,
#                                 breaks = quantile(meta_train$LARS, c(0, 1/3, 2/3, 1)),
#                                 include.lowest = TRUE,
#                                 labels = c("Low", "Medium", "High"))
# fit_lars <- survfit(Surv(survival_time, event) ~ LARS_tertile, data = meta_train)
# print(ggsurvplot(fit_lars, pval = TRUE, risk.table = TRUE,
#                  title = "LARS Tertile Survival"))

# 外部验证（如果有GSE25504的生存数据）
# ... 同样计算外部C-index

# =============================================================================
# M5: GSEA FDR<0.05 严格阈值
# =============================================================================

cat("\n===== M5: GSEA with FDR<0.05 =====\n")

# 重新运行GSEA，同时报告FDR<0.25和FDR<0.05
# gsea_result <- GSEA(
#   geneList = gene_list,
#   TERM2GENE = msigdb_c2,  # MSigDB C2 (KEGG + Reactome)
#   pvalueCutoff = 0.25,
#   pAdjustMethod = "BH",
#   minGSSize = 10,
#   maxGSSize = 500,
#   seed = 42
# )

# 筛选FDR<0.05的通路
# gsea_strict <- gsea_result[gsea_result$qvalues < 0.05, ]
# cat("\nPathways at FDR<0.25:", nrow(gsea_result@geneSetSizes), "\n")
# cat("Pathways at FDR<0.05:", nrow(gsea_strict), "\n")

# 关键通路在两个阈值下的稳健性
# key_pathways <- c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
#                    "KEGG_HIF1_SIGNALING_PATHWAY",
#                    "KEGG_PYRUVATE_METABOLISM",
#                    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#                    "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY")
#
# cat("\nKey pathway robustness check:\n")
# for (pw in key_pathways) {
#   idx <- which(gsea_result$Description == pw)
#   if (length(idx) > 0) {
#     cat(pw, ": NES =", gsea_result$NES[idx],
#         ", FDR =", gsea_result$qvalues[idx], "\n")
#   }
# }

# 输出Supplementary Table S1
# write.csv(
#   data.frame(
#     Pathway = gsea_result$Description,
#     NES = gsea_result$NES,
#     pvalue = gsea_result$pvalue,
#     FDR = gsea_result$qvalues,
#     Significant_at_FDR0.05 = gsea_result$qvalues < 0.05
#   ),
#   file = "Supplementary_Table_S1_GSEA_comparison.csv",
#   row.names = FALSE
# )

# =============================================================================
# m2: CD14+ HDAC9梯度供体水平统计检验
# =============================================================================

cat("\n===== m2: CD14+ HDAC9 Donor-Level Statistics =====\n")

# 从单细胞数据提取供体水平HDAC9均值
# 假设 sc_meta 包含: cell_id, donor_id, group (HC/Surv/NonSurv), HDAC9_score
#
# donor_hdac9 <- aggregate(HDAC9_score ~ donor_id + group,
#                          data = sc_meta[sc_meta$cell_type == "CD14+ Monocytes", ],
#                          FUN = mean)
# cat("\nDonor-level HDAC9 means:\n")
# print(donor_hdac9)

# 线性混合效应模型
# lmm_result <- lmer(HDAC9_score ~ group + (1|donor_id),
#                    data = sc_meta[sc_meta$cell_type == "CD14+ Monocytes", ])
# cat("\nLinear mixed-effects model:\n")
# print(summary(lmm_result))

# 趋势检验（group作为有序因子）
# donor_hdac9$group_ordered <- factor(donor_hdac9$group,
#                                      levels = c("HC", "Surv", "NonSurv"),
#                                      ordered = TRUE)
# trend_test <- cor.test(as.numeric(donor_hdac9$group_ordered),
#                        donor_hdac9$HDAC9_score,
#                        method = "spearman")
# cat("\nSpearman trend test (donor-level):\n")
# cat("  rho =", trend_test$estimate, ", p =", trend_test$p.value, "\n")

# 配对Wilcoxon检验（供体水平，n=4 per group）
# hc_vals <- donor_hdac9$HDAC9_score[donor_hdac9$group == "HC"]
# surv_vals <- donor_hdac9$HDAC9_score[donor_hdac9$group == "Surv"]
# nonsurv_vals <- donor_hdac9$HDAC9_score[donor_hdac9$group == "NonSurv"]
#
# cat("\nPairwise Wilcoxon (donor-level, n=4 each):\n")
# cat("HC vs Non-survivor:",
#     wilcox.test(hc_vals, nonsurv_vals, exact = FALSE)$p.value, "\n")
# cat("Survivor vs Non-survivor:",
#     wilcox.test(surv_vals, nonsurv_vals, exact = FALSE)$p.value, "\n")
# cat("HC vs Survivor:",
#     wilcox.test(hc_vals, surv_vals, exact = FALSE)$p.value, "\n")

# 同样对LDHA做相同分析
# donor_ldha <- aggregate(LDHA_score ~ donor_id + group,
#                         data = sc_meta[sc_meta$cell_type == "CD14+ Monocytes", ],
#                         FUN = mean)
# ... 同上

# 输出Supplementary Table S2
# write.csv(
#   data.frame(
#     Test = c("LMM_group_effect", "Spearman_trend",
#              "Wilcox_HC_vs_NonSurv", "Wilcox_Surv_vs_NonSurv", "Wilcox_HC_vs_Surv"),
#     Statistic = c(...),
#     P_value = c(...),
#     Sample_size = c(...)
#   ),
#   file = "Supplementary_Table_S2_scRNA_stats.csv",
#   row.names = FALSE
# )

cat("\n===== All analyses complete =====\n")
cat("Please fill in the [TO BE ADDED AFTER R ANALYSIS] placeholders in Manuscript.md\n")
cat("with the results from this script.\n")
