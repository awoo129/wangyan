#!/usr/bin/env Rscript
# 儿童脓毒症免疫瘫痪研究 - GEO数据下载与预处理
# Author: Bioinformatic Analysis Pipeline
# Date: 2024

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(Biobase)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(ggdendro)
  library(factoextra)
  library(SummarizedExperiment)
  library(matrixStats)
})

# 设置工作目录
main_dir <- "./长期计划/儿童脓毒症免疫瘫痪研究"
setwd(main_dir)

# 创建子目录（如果不存在）
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/normalized", showWarnings = FALSE, recursive = TRUE)
dir.create("data/metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("results/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("results/preprocessing", showWarnings = FALSE, recursive = TRUE)

# 日志记录函数
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

# 错误记录函数
log_error <- function(msg) {
  cat(paste0("[", Sys.time(), "] ERROR: ", msg, "\n"))
}

# ========================
# 1. 数据下载
# ========================
download_gse_data <- function(gse_id, dest_dir = "data/raw") {
  log_message(paste("开始下载数据集:", gse_id))
  
  matrix_file <- file.path(dest_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  soft_file <- file.path(dest_dir, paste0(gse_id, "_soft.RData"))
  
  # 检查是否已存在
  if (file.exists(matrix_file) && file.exists(soft_file)) {
    log_message(paste(gse_id, "已存在，跳过下载"))
    return(TRUE)
  }
  
  tryCatch({
    # 下载series matrix文件
    geo_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gse_id, "&format=file")
    
    # 使用GEOquery下载
    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE, destdir = dest_dir)
    
    if (is.null(gse) || length(gse) == 0) {
      log_error(paste(gse_id, "下载失败，返回为空"))
      return(NULL)
    }
    
    log_message(paste(gse_id, "下载成功"))
    return(gse)
    
  }, error = function(e) {
    log_error(paste(gse_id, "下载失败:", e$message))
    return(NULL)
  })
}

# ========================
# 2. 数据预处理
# ========================
preprocess_expression_data <- function(gse_object, gse_id) {
  log_message(paste("开始预处理:", gse_id))
  
  # 获取表达矩阵
  expr_data <- exprs(gse_object[[1]])
  
  # 获取样本信息
  sample_info <- pData(gse_object[[1]])
  
  # 获取基因注释
  feature_info <- fData(gse_object[[1]])
  
  # 获取平台信息
  platform <- annotation(gse_object[[1]])
  
  result <- list(
    expression = expr_data,
    sample_info = sample_info,
    feature_info = feature_info,
    platform = platform,
    gse_id = gse_id
  )
  
  return(result)
}

# ========================
# 3. 探针到基因符号映射
# ========================
map_probes_to_genes <- function(preprocessed_data, gpl_info = NULL) {
  expr_matrix <- preprocessed_data$expression
  feature_info <- preprocessed_data$feature_info
  
  # 查找基因符号列
  gene_symbol_col <- NULL
  
  # 常见列名
  possible_names <- c("Gene Symbol", "gene_symbol", "SYMBOL", "symbol", 
                      "GENE_SYMBOL", "GeneSymbol", "gene")
  
  for (name in possible_names) {
    if (name %in% colnames(feature_info)) {
      gene_symbol_col <- name
      break
    }
  }
  
  if (is.null(gene_symbol_col)) {
    log_message("未找到基因符号列，尝试其他方法...")
    # 返回原始矩阵
    return(expr_matrix)
  }
  
  log_message(paste("使用列进行映射:", gene_symbol_col))
  
  # 创建探针-基因映射
  probe_to_gene <- feature_info[[gene_symbol_col]]
  names(probe_to_gene) <- rownames(feature_info)
  
  # 过滤掉空的基因符号
  valid_probes <- !is.na(probe_to_gene) & probe_to_gene != "" & probe_to_gene != "NA"
  probe_to_gene <- probe_to_gene[valid_probes]
  
  expr_filtered <- expr_matrix[names(probe_to_gene), ]
  
  # 按基因符号聚合（取均值）
  expr_aggregated <- apply(expr_filtered, 2, function(x) {
    tapply(x, probe_to_gene, mean, na.rm = TRUE)
  })
  
  return(expr_aggregated)
}

# ========================
# 4. MHC II类基因提取
# ========================
extract_mhc_genes <- function(expr_matrix, mhc_genes = c("HLA-DRA", "HLA-DQB1", "HLA-DQA1", 
                                                          "CIITA", "CD74", "HLA-DRB1", 
                                                          "HLA-DPA1", "HLA-DPB1", "HLA-DMA", 
                                                          "HLA-DMB")) {
  log_message("提取MHC II类基因...")
  
  available_genes <- rownames(expr_matrix)
  
  found_genes <- mhc_genes[mhc_genes %in% available_genes]
  not_found_genes <- mhc_genes[!mhc_genes %in% available_genes]
  
  if (length(not_found_genes) > 0) {
    log_message(paste("未找到的基因:", paste(not_found_genes, collapse = ", ")))
  }
  
  if (length(found_genes) == 0) {
    log_message("警告: 未找到任何MHC II类基因!")
    return(NULL)
  }
  
  mhc_expr <- expr_matrix[found_genes, , drop = FALSE]
  log_message(paste("成功提取", length(found_genes), "个MHC II类基因"))
  
  return(list(
    expression = mhc_expr,
    found_genes = found_genes,
    not_found_genes = not_found_genes
  ))
}

# ========================
# 5. 数据质控
# ========================
perform_qc <- function(expr_matrix, sample_info, gse_id, output_dir = "results/QC") {
  log_message(paste("开始质控分析:", gse_id))
  
  # 创建质控结果目录
  qc_dir <- file.path(output_dir, gse_id)
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 5.1 样本分布检查 - 箱线图
  pdf(file.path(qc_dir, paste0(gse_id, "_boxplot.pdf")))
  boxplot(expr_matrix, outline = FALSE, main = paste(gse_id, "样本表达分布"),
          xlab = "Samples", ylab = "Expression", las = 2)
  dev.off()
  
  # 5.2 密度图
  pdf(file.path(qc_dir, paste0(gse_id, "_density.pdf")))
  plot(density(expr_matrix[, 1]), main = paste(gse_id, "表达密度分布"),
       xlab = "Expression", col = 1)
  for (i in 2:ncol(expr_matrix)) {
    lines(density(expr_matrix[, i]), col = i)
  }
  dev.off()
  
  # 5.3 PCA分析
  pca_result <- tryCatch({
    pca <- prcomp(t(expr_matrix), scale. = TRUE)
    pca
  }, error = function(e) {
    log_error(paste(gse_id, "PCA分析失败:", e$message))
    return(NULL)
  })
  
  if (!is.null(pca_result)) {
    pdf(file.path(qc_dir, paste0(gse_id, "_pca.pdf")))
    print(fviz_pca_ind(pca_result, geom = "point", 
                       title = paste(gse_id, "PCA分析")) +
            theme_minimal())
    dev.off()
    
    # 保存PCA数据
    pca_data <- as.data.frame(pca_result$x)
    pca_data$sample <- rownames(pca_data)
  }
  
  # 5.4 层次聚类
  dist_matrix <- dist(t(expr_matrix), method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  pdf(file.path(qc_dir, paste0(gse_id, "_dendrogram.pdf")))
  plot(hc, main = paste(gse_id, "样本聚类树"), xlab = "Samples")
  dev.off()
  
  # 5.5 质控统计
  qc_stats <- data.frame(
    Metric = c("Total Samples", "Total Genes", "Mean Expression", "Median Expression",
               "Min Expression", "Max Expression", "SD"),
    Value = c(ncol(expr_matrix), nrow(expr_matrix),
              round(mean(expr_matrix), 2), round(median(expr_matrix), 2),
              round(min(expr_matrix), 2), round(max(expr_matrix), 2),
              round(sd(expr_matrix), 2))
  )
  
  write.csv(qc_stats, file.path(qc_dir, paste0(gse_id, "_qc_stats.csv")), row.names = FALSE)
  
  log_message(paste(gse_id, "质控完成"))
  
  return(list(
    pca = pca_result,
    hc = hc,
    stats = qc_stats,
    qc_dir = qc_dir
  ))
}

# ========================
# 6. 批次效应校正（如果需要）
# ========================
combat_correction <- function(expr_matrix, batch) {
  log_message("执行批次效应校正...")
  
  if (length(unique(batch)) < 2) {
    log_message("只有一个批次，跳过校正")
    return(expr_matrix)
  }
  
  library(sva)
  
  corrected <- ComBat(dat = expr_matrix, batch = batch, par.prior = TRUE)
  
  return(corrected)
}

# ========================
# 主程序
# ========================
main <- function() {
  log_message("========================================")
  log_message("儿童脓毒症免疫瘫痪研究 - GEO数据预处理")
  log_message("========================================")
  
  # 定义数据集
  gse_ids <- c("GSE26378", "GSE26440", "GSE13904")
  
  # 存储所有处理后的数据
  all_data <- list()
  all_metadata <- list()
  mhc_results <- list()
  qc_results <- list()
  
  # 处理每个数据集
  for (gse_id in gse_ids) {
    log_message(paste("\n========== 处理数据集:", gse_id, "=========="))
    
    # 下载数据
    gse_data <- download_gse_data(gse_id)
    
    if (is.null(gse_data)) {
      log_error(paste(gse_id, "下载失败，跳过"))
      next
    }
    
    # 预处理
    preprocessed <- preprocess_expression_data(gse_data, gse_id)
    
    # 保存原始数据
    saveRDS(preprocessed$expression, file.path("data/raw", paste0(gse_id, "_raw_expr.rds")))
    saveRDS(preprocessed$sample_info, file.path("data/raw", paste0(gse_id, "_sample_info.rds")))
    
    # 探针到基因映射
    expr_mapped <- map_probes_to_genes(preprocessed)
    
    # 保存映射后数据
    saveRDS(expr_mapped, file.path("data/normalized", paste0(gse_id, "_mapped_expr.rds")))
    
    # 提取样本分组信息
    sample_info <- preprocessed$sample_info
    
    # 查找分组信息
    group_col <- NULL
    possible_cols <- c("characteristics_ch1", "source_name_ch1", "title", 
                       "disease:ch1", "group:ch1")
    
    for (col in possible_cols) {
      if (col %in% colnames(sample_info)) {
        group_col <- col
        break
      }
    }
    
    if (!is.null(group_col)) {
      groups <- as.character(sample_info[[group_col]])
      groups[grepl("normal|healthy|control", tolower(groups))] <- "Normal"
      groups[grepl("sepsis|septic", tolower(groups))] <- "Sepsis"
    } else {
      groups <- rep("Unknown", nrow(sample_info))
    }
    
    # 创建元数据表
    metadata <- data.frame(
      sample_id = rownames(sample_info),
      group = groups,
      gse_id = gse_id,
      platform = preprocessed$platform,
      stringsAsFactors = FALSE
    )
    
    write.csv(metadata, file.path("data/metadata", paste0(gse_id, "_metadata.csv")), row.names = FALSE)
    
    # 提取MHC II类基因
    mhc_data <- extract_mhc_genes(expr_mapped)
    
    if (!is.null(mhc_data)) {
      mhc_results[[gse_id]] <- mhc_data
      
      # 保存MHC表达矩阵
      saveRDS(mhc_data$expression, file.path("data/normalized", paste0(gse_id, "_mhc_expr.rds")))
    }
    
    # 执行质控
    qc_result <- perform_qc(expr_mapped, metadata, gse_id)
    qc_results[[gse_id]] <- qc_result
    
    # 保存处理后的数据
    all_data[[gse_id]] <- list(
      expression = expr_mapped,
      metadata = metadata,
      platform = preprocessed$platform
    )
    
    # 打印样本信息
    log_message(paste(gse_id, "- 样本数量:", ncol(expr_mapped)))
    log_message(paste(gse_id, "- 基因数量:", nrow(expr_mapped)))
    log_message(paste(gse_id, "- 分组:", paste(table(groups), collapse = ", ")))
    
    # 清理内存
    rm(gse_data)
    gc()
  }
  
  # ========================
  # 生成汇总报告
  # ========================
  log_message("\n========== 生成汇总报告 ==========")
  
  # 创建汇总表
  summary_table <- data.frame(
    GSE_ID = gse_ids,
    Platform = sapply(gse_ids, function(x) {
      if (!is.null(all_data[[x]])) all_data[[x]]$platform else NA
    }),
    Total_Samples = sapply(gse_ids, function(x) {
      if (!is.null(all_data[[x]])) ncol(all_data[[x]]$expression) else NA
    }),
    Normal_Samples = sapply(gse_ids, function(x) {
      if (!is.null(all_data[[x]])) sum(all_data[[x]]$metadata$group == "Normal") else NA
    }),
    Sepsis_Samples = sapply(gse_ids, function(x) {
      if (!is.null(all_data[[x]])) sum(all_data[[x]]$metadata$group == "Sepsis") else NA
    }),
    MHC_Genes_Found = sapply(gse_ids, function(x) {
      if (!is.null(mhc_results[[x]])) length(mhc_results[[x]]$found_genes) else NA
    }),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_table, "results/preprocessing/dataset_summary.csv", row.names = FALSE)
  
  # MHC基因检测汇总
  all_found_genes <- unique(unlist(lapply(mhc_results, function(x) x$found_genes)))
  all_missing_genes <- unique(unlist(lapply(mhc_results, function(x) x$not_found_genes)))
  
  mhc_summary <- data.frame(
    Gene = c(all_found_genes, all_missing_genes),
    Status = c(rep("Found", length(all_found_genes)), rep("Not Found", length(all_missing_genes))),
    stringsAsFactors = FALSE
  )
  
  write.csv(mhc_summary, "results/preprocessing/mhc_gene_detection.csv", row.names = FALSE)
  
  # HLA-DRA验证
  hla_dra_present <- sapply(gse_ids, function(x) {
    if (!is.null(mhc_results[[x]])) "HLA-DRA" %in% mhc_results[[x]]$found_genes else FALSE
  })
  
  log_message("\n========== HLA-DRA检测结果 ==========")
  for (i in seq_along(gse_ids)) {
    status <- ifelse(hla_dra_present[i], "✓ 存在", "✗ 不存在")
    log_message(paste(gse_ids[i], ":", status))
  }
  
  # 保存完整结果
  saveRDS(all_data, "data/normalized/all_data.rds")
  saveRDS(mhc_results, "data/normalized/mhc_results.rds")
  saveRDS(qc_results, "results/QC/qc_results.rds")
  
  log_message("\n========================================")
  log_message("数据预处理完成!")
  log_message("========================================")
  
  # 返回结果
  return(list(
    summary = summary_table,
    mhc_summary = mhc_summary,
    hla_dra_status = hla_dra_present
  ))
}

# 执行主程序
results <- main()

# 打印最终结果
print(results$summary)
print(results$mhc_summary)
