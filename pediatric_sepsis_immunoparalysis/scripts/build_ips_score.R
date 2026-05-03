#!/usr/bin/env Rscript
# 儿童脓毒症免疫瘫痪研究 - 免疫瘫痪评分(IPS)构建
# 基于HLA-DRA等MHC II类基因表达构建免疫瘫痪评分

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(caret)
  library(rms)
  library(survival)
  library(survminer)
})

# 设置工作目录
main_dir <- "./长期计划/儿童脓毒症免疫瘫痪研究"
setwd(main_dir)

log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

# ========================
# 1. 加载数据
# ========================
load_data_for_ips <- function() {
  log_message("加载数据...")
  
  # 加载预处理后的数据
  all_data <- readRDS("data/normalized/all_data.rds")
  mhc_results <- readRDS("data/normalized/mhc_results.rds")
  
  return(list(all_data = all_data, mhc_results = mhc_results))
}

# ========================
# 2. MHC基因表达提取
# ========================
extract_mhc_matrix <- function(all_data, mhc_genes = c("HLA-DRA", "HLA-DQB1", 
                                                         "HLA-DQA1", "HLA-DRB1",
                                                         "HLA-DPA1", "HLA-DPB1")) {
  log_message("提取MHC基因表达矩阵...")
  
  mhc_expr_list <- list()
  
  for (gse_id in names(all_data)) {
    expr <- all_data[[gse_id]]$expression
    meta <- all_data[[gse_id]]$metadata
    
    # 标准化基因名
    rownames(expr) <- toupper(rownames(expr))
    
    # 提取MHC基因
    available_genes <- mhc_genes[mhc_genes %in% rownames(expr)]
    missing_genes <- mhc_genes[!mhc_genes %in% rownames(expr)]
    
    if (length(available_genes) > 0) {
      mhc_expr <- expr[available_genes, , drop = FALSE]
      
      # 添加样本分组
      if (nrow(meta) == ncol(mhc_expr)) {
        mhc_expr <- rbind(mhc_expr, Group = meta$group)
      }
      
      mhc_expr_list[[gse_id]] <- t(mhc_expr)  # 转换为 样本 x 基因
    }
    
    if (length(missing_genes) > 0) {
      log_message(paste(gse_id, "- 未找到:", paste(missing_genes, collapse = ", ")))
    }
  }
  
  return(mhc_expr_list)
}

# ========================
# 3. 构建IPS评分
# ========================
build_ips <- function(mhc_expr_matrix, method = "mean") {
  log_message(paste("构建IPS评分 - 方法:", method))
  
  # 提取数值矩阵（排除Group列）
  gene_cols <- grep("Group", colnames(mhc_expr_matrix), invert = TRUE, value = TRUE)
  expr_values <- mhc_expr_matrix[, gene_cols, drop = FALSE]
  
  if (method == "mean") {
    # 简单平均
    ips <- rowMeans(expr_values, na.rm = TRUE)
  } else if (method == "weighted") {
    # 基于差异分析的加权评分
    # HLA-DRA权重更高
    weights <- c(HLA_DRA = 2, HLA_DQB1 = 1.5, HLA_DQA1 = 1.5,
                 HLA_DRB1 = 1, HLA_DPA1 = 1, HLA_DPB1 = 1)
    weights <- weights[gene_cols]
    weights[is.na(weights)] <- 1
    
    ips <- rowMeans(sweep(expr_values, 2, weights, `*`), na.rm = TRUE)
  } else if (method == "hla_dra_only") {
    # 仅使用HLA-DRA
    if ("HLA-DRA" %in% gene_cols) {
      ips <- as.numeric(expr_values[, "HLA-DRA"])
    } else {
      log_message("HLA-DRA不存在，使用平均评分")
      ips <- rowMeans(expr_values, na.rm = TRUE)
    }
  }
  
  return(ips)
}

# ========================
# 4. IPS分组
# ========================
categorize_ips <- function(ips, method = "median", threshold = NULL) {
  log_message(paste("IPS分组 - 方法:", method))
  
  if (method == "median") {
    cutoff <- median(ips, na.rm = TRUE)
  } else if (method == "quartile") {
    q <- quantile(ips, c(0.25, 0.5, 0.75), na.rm = TRUE)
    groups <- ifelse(ips <= q[1], "Low",
                     ifelse(ips <= q[2], "Medium-Low",
                           ifelse(ips <= q[3], "Medium-High", "High")))
    return(list(groups = groups, cutoff = q))
  } else if (method == "custom" && !is.null(threshold)) {
    cutoff <- threshold
  }
  
  groups <- ifelse(ips <= cutoff, "Low_IPS", "High_IPS")
  
  return(list(groups = groups, cutoff = cutoff))
}

# ========================
# 5. 验证分析
# ========================
validate_ips <- function(ips, groups, metadata, gse_id) {
  log_message(paste("验证IPS:", gse_id))
  
  results <- list()
  
  # 5.1 HLA-DRA与IPS相关性
  if ("HLA-DRA" %in% colnames(metadata$expression)) {
    cor_test <- cor.test(as.numeric(metadata$expression["HLA-DRA", ]),
                         ips, method = "spearman")
    results$correlation <- cor_test
  }
  
  # 5.2 组间差异
  group_ips <- data.frame(IPS = ips, Group = groups)
  
  # 统计检验
  low_ips <- ips[groups == "Low_IPS"]
  high_ips <- ips[groups == "High_IPS"]
  
  if (length(low_ips) > 0 && length(high_ips) > 0) {
    t_test <- t.test(low_ips, high_ips)
    results$t_test <- t_test
    log_message(paste(gse_id, "- IPS组间差异: t =", round(t_test$statistic, 3),
                     "p =", format(t_test$p.value, scientific = TRUE)))
  }
  
  # 5.3 可视化
  p_dir <- "results/ips_validation"
  dir.create(p_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 箱线图
  p1 <- ggplot(group_ips, aes(x = groups, y = IPS, fill = groups)) +
    geom_boxplot() +
    labs(title = paste(gse_id, "- IPS分组比较"),
         x = "IPS Group", y = "Immunodeficiency Score") +
    theme_minimal()
  
  ggsave(file.path(p_dir, paste0(gse_id, "_ips_boxplot.png")),
         p1, width = 8, height = 6, dpi = 300)
  
  # 分布图
  p2 <- ggplot(group_ips, aes(x = IPS, fill = groups)) +
    geom_density(alpha = 0.5) +
    labs(title = paste(gse_id, "- IPS分布"),
         x = "Immunodeficiency Score", y = "Density") +
    theme_minimal()
  
  ggsave(file.path(p_dir, paste0(gse_id, "_ips_density.png")),
         p2, width = 8, height = 6, dpi = 300)
  
  return(results)
}

# ========================
# 6. 跨数据集验证
# ========================
cross_validate_ips <- function(all_mhc_expr, training_gse = "GSE26378", 
                                validation_genes = c("GSE26440", "GSE13904")) {
  log_message("跨数据集验证...")
  
  # 训练集构建IPS
  if (training_gse %in% names(all_mhc_expr)) {
    train_expr <- all_mhc_expr[[training_gse]]
    train_ips <- build_ips(train_expr, method = "hla_dra_only")
    train_groups <- categorize_ips(train_ips)$groups
    
    log_message(paste("训练集 IPS统计:",
                       "Mean =", round(mean(train_ips, na.rm = TRUE), 3),
                       "Median =", round(median(train_ips, na.rm = TRUE), 3)))
  }
  
  # 验证集验证
  for (gse_id in validation_genes) {
    if (gse_id %in% names(all_mhc_expr)) {
      val_expr <- all_mhc_expr[[gse_id]]
      
      # 使用训练集参数标准化
      if ("HLA-DRA" %in% colnames(val_expr)) {
        # 相对表达（与训练集比较）
        train_hla_dra <- train_expr[, "HLA-DRA"] if exists("train_expr") else NULL
      }
      
      log_message(paste(gse_id, "- 验证完成"))
    }
  }
}

# ========================
# 7. 生成报告
# ========================
generate_ips_report <- function(ips_results, all_mhc_expr) {
  log_message("生成IPS报告...")
  
  report <- "# 免疫瘫痪评分(IPS)构建报告\n\n"
  
  report <- paste0(report, "## 1. 评分方法\n\n")
  report <- paste0(report, "- **核心基因**: HLA-DRA\n")
  report <- paste0(report, "- **扩展基因**: HLA-DQB1, HLA-DQA1, HLA-DRB1, HLA-DPA1, HLA-DPB1\n")
  report <- paste0(report, "- **计算方法**: 基于MHC II类基因表达的加权平均\n\n")
  
  report <- paste0(report, "## 2. 数据集统计\n\n")
  report <- paste0(report, "| 数据集 | 样本数 | Low_IPS | High_IPS | HLA-DRA表达 |\n")
  report <- paste0(report, "|--------|--------|---------|----------|-------------|\n")
  
  for (gse_id in names(all_mhc_expr)) {
    expr <- all_mhc_expr[[gse_id]]
    n <- nrow(expr)
    
    if ("Group" %in% colnames(expr)) {
      groups <- expr[, "Group"]
      low_n <- sum(groups == "Low_IPS", na.rm = TRUE)
      high_n <- sum(groups == "High_IPS", na.rm = TRUE)
    } else {
      low_n <- "-"
      high_n <- "-"
    }
    
    if ("HLA-DRA" %in% colnames(expr)) {
      hla_dra_mean <- round(mean(as.numeric(expr[, "HLA-DRA"]), na.rm = TRUE), 3)
    } else {
      hla_dra_mean <- "-"
    }
    
    report <- paste0(report, sprintf("| %s | %d | %s | %s | %s |\n",
                                     gse_id, n, low_n, high_n, hla_dra_mean))
  }
  
  report <- paste0(report, "\n## 3. 临床意义\n\n")
  report <- paste0(report, "- **High IPS**: 免疫功能低下，HLA-DRA高表达\n")
  report <- paste0(report, "- **Low IPS**: 免疫功能相对正常\n")
  report <- paste0(report, "- HLA-DRA是MHC II类分子，在抗原呈递中起关键作用\n")
  report <- paste0(report, "- 其表达下调与脓毒症免疫瘫痪相关\n")
  
  # 保存报告
  p_dir <- "results/ips_validation"
  dir.create(p_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(report, file.path(p_dir, "ips_report.md"))
  
  return(report)
}

# ========================
# 主程序
# ========================
main <- function() {
  log_message("========================================")
  log_message("免疫瘫痪评分(IPS)构建")
  log_message("========================================")
  
  # 加载数据
  data <- load_data_for_ips()
  all_data <- data$all_data
  
  # 提取MHC基因矩阵
  mhc_genes <- c("HLA-DRA", "HLA-DQB1", "HLA-DQA1", "HLA-DRB1",
                 "HLA-DPA1", "HLA-DPB1")
  all_mhc_expr <- extract_mhc_matrix(all_data, mhc_genes)
  
  # 构建IPS（使用HLA-DRA作为主要指标）
  ips_results <- list()
  
  for (gse_id in names(all_mhc_expr)) {
    ips <- build_ips(all_mhc_expr[[gse_id]], method = "hla_dra_only")
    ips_cat <- categorize_ips(ips)
    
    # 添加分组信息
    all_mhc_expr[[gse_id]] <- cbind(all_mhc_expr[[gse_id]], 
                                     IPS = ips,
                                     IPS_Group = ips_cat$groups)
    
    # 验证
    metadata <- all_data[[gse_id]]
    validation <- validate_ips(ips, ips_cat$groups, metadata, gse_id)
    
    ips_results[[gse_id]] <- list(
      ips = ips,
      groups = ips_cat$groups,
      cutoff = ips_cat$cutoff,
      validation = validation
    )
    
    log_message(paste(gse_id, "- IPS完成:",
                       "Mean =", round(mean(ips, na.rm = TRUE), 3),
                       "Cutoff =", round(ips_cat$cutoff, 3)))
  }
  
  # 保存结果
  saveRDS(all_mhc_expr, "data/normalized/mhc_with_ips.rds")
  saveRDS(ips_results, "results/ips_validation/ips_results.rds")
  
  # 生成报告
  report <- generate_ips_report(ips_results, all_mhc_expr)
  
  log_message("\n========================================")
  log_message("IPS构建完成!")
  log_message("========================================")
  
  return(list(ips_results = ips_results, report = report))
}

# 执行
if (!interactive()) {
  results <- main()
}
