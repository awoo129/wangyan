#!/usr/bin/env Rscript
# 儿童脓毒症免疫瘫痪研究 - 差异表达分析
# 运行条件: 需要先运行 download_and_preprocess.py 完成数据预处理

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# 设置工作目录
main_dir <- "./长期计划/儿童脓毒症免疫瘫痪研究"
setwd(main_dir)

# 日志函数
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

# ========================
# 1. 加载预处理数据
# ========================
load_preprocessed_data <- function() {
  log_message("加载预处理数据...")
  
  # 加载所有数据
  all_data <- readRDS("data/normalized/all_data.rds")
  
  return(all_data)
}

# ========================
# 2. 差异表达分析
# ========================
perform_dea <- function(expr_matrix, metadata, contrast_name = "Normal_vs_Sepsis") {
  log_message(paste("差异表达分析:", contrast_name))
  
  # 确保样本顺序一致
  expr_matrix <- expr_matrix[, metadata$sample_id]
  
  # 创建设计矩阵
  group <- factor(metadata$group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # 拟合模型
  fit <- lmFit(expr_matrix, design)
  
  # 定义对比
  if ("Normal" %in% colnames(design) && "Sepsis" %in% colnames(design)) {
    contrast.matrix <- makeContrasts(Normal - Sepsis, levels = design)
  } else {
    log_message("警告: 分组信息不完整，跳过差异分析")
    return(NULL)
  }
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 获取差异基因
  results <- topTable(fit2, adjust.method = "BH", number = Inf)
  results$gene <- rownames(results)
  
  # 分类
  results$regulation <- "Not Significant"
  results$regulation[results$logFC > 0.5 & results$adj.P.Val < 0.05] <- "Up"
  results$regulation[results$logFC < -0.5 & results$adj.P.Val < 0.05] <- "Down"
  
  log_message(paste("差异基因数量:",
                    "Up =", sum(results$regulation == "Up"),
                    "Down =", sum(results$regulation == "Down")))
  
  return(results)
}

# ========================
# 3. MHC基因差异分析
# ========================
analyze_mhc_dea <- function(dea_results, mhc_genes) {
  log_message("MHC基因差异分析...")
  
  mhc_dea <- dea_results[dea_results$gene %in% mhc_genes, ]
  
  if (nrow(mhc_dea) == 0) {
    log_message("警告: MHC基因在差异分析结果中未找到")
    return(NULL)
  }
  
  # HLA-DRA 重点关注
  hla_dra <- mhc_dea[mhc_dea$gene == "HLA-DRA", ]
  if (nrow(hla_dra) > 0) {
    log_message(paste("HLA-DRA:",
                       "logFC =", round(hla_dra$logFC, 3),
                       "adj.P.Val =", format(hla_dra$adj.P.Val, scientific = TRUE)))
  }
  
  return(mhc_dea)
}

# ========================
# 4. 可视化
# ========================
visualize_dea <- function(dea_results, gse_id, output_dir = "results/differential_expression") {
  log_message(paste("生成差异分析图表:", gse_id))
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 4.1 火山图
  dea_results$neg_log10_pval <- -log10(dea_results$adj.P.Val)
  
  volcano_plot <- ggplot(dea_results, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", 
                                   "Not Significant" = "grey80")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
    labs(title = paste(gse_id, "- 差异表达分析火山图"),
         x = "log2 Fold Change", y = "-log10(adjusted P-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 标注关键基因
  key_genes <- c("HLA-DRA", "HLA-DQB1", "HLA-DQA1", "CIITA", "CD74")
  key_data <- dea_results[dea_results$gene %in% key_genes, ]
  
  if (nrow(key_data) > 0) {
    volcano_plot <- volcano_plot + 
      geom_text_repel(data = key_data, aes(label = gene), size = 3)
  }
  
  ggsave(file.path(output_dir, paste0(gse_id, "_volcano.png")), 
         volcano_plot, width = 10, height = 8, dpi = 300)
  
  # 4.2 热图 - MHC基因
  log_message("生成MHC基因热图...")
  
  # 读取表达数据
  expr_file <- file.path("data/normalized", paste0(gse_id, "_expression.csv"))
  if (file.exists(expr_file)) {
    expr_data <- read.csv(expr_file, row.names = 1)
    
    # 提取MHC基因
    mhc_genes <- c("HLA-DRA", "HLA-DQB1", "HLA-DQA1", "HLA-DRB1", 
                   "HLA-DPA1", "HLA-DPB1", "CIITA", "CD74")
    mhc_available <- mhc_genes[mhc_genes %in% rownames(expr_data)]
    
    if (length(mhc_available) > 0) {
      mhc_expr <- expr_data[mhc_available, ]
      
      # Z-score标准化
      mhc_zscore <- t(scale(t(mhc_expr)))
      
      # 读取元数据
      meta_file <- file.path("data/metadata", paste0(gse_id, "_metadata.csv"))
      if (file.exists(meta_file)) {
        metadata <- read.csv(meta_file)
        annotation_col <- data.frame(Group = metadata$group)
        rownames(annotation_col) <- metadata$sample_id
        
        pheatmap(mhc_zscore, 
                 annotation_col = annotation_col,
                 show_colnames = FALSE,
                 main = paste(gse_id, "- MHC II类基因表达热图"),
                 filename = file.path(output_dir, paste0(gse_id, "_mhc_heatmap.png")),
                 width = 12, height = 6)
      }
    }
  }
  
  # 4.3 MHC基因箱线图
  log_message("生成HLA-DRA箱线图...")
  
  hla_dra_plot <- ggplot(dea_results[dea_results$gene %in% key_genes, ], 
                         aes(x = gene, y = logFC, fill = regulation)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste(gse_id, "- MHC II类基因差异表达"),
         x = "Gene", y = "log2 Fold Change (Normal vs Sepsis)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0(gse_id, "_mhc_barplot.png")), 
         hla_dra_plot, width = 8, height = 6, dpi = 300)
}

# ========================
# 5. GO/KEGG富集分析
# ========================
perform_enrichment <- function(dea_results, output_dir = "results/enrichment") {
  log_message("执行富集分析...")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 差异基因
  deg_genes <- dea_results$gene[dea_results$regulation != "Not Significant"]
  
  if (length(deg_genes) < 10) {
    log_message("差异基因数量不足，跳过富集分析")
    return(NULL)
  }
  
  # GO分析
  tryCatch({
    go_results <- enrichGO(gene = deg_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)
    
    if (!is.null(go_results) && nrow(go_results) > 0) {
      write.csv(as.data.frame(go_results), 
                file.path(output_dir, "go_enrichment.csv"), row.names = FALSE)
      
      # 可视化
      dotplot(go_results, showCategory = 20) +
        ggtitle("GO Enrichment - Biological Process") +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(output_dir, "go_enrichment.png"), 
             width = 12, height = 10, dpi = 300)
    }
  }, error = function(e) {
    log_message(paste("GO分析失败:", e$message))
  })
  
  # KEGG分析
  tryCatch({
    kegg_results <- enrichKEGG(gene = deg_genes,
                               organism = 'hsa',
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
    
    if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
      write.csv(as.data.frame(kegg_results), 
                file.path(output_dir, "kegg_enrichment.csv"), row.names = FALSE)
    }
  }, error = function(e) {
    log_message(paste("KEGG分析失败:", e$message))
  })
}

# ========================
# 主程序
# ========================
main <- function() {
  log_message("========================================")
  log_message("差异表达分析")
  log_message("========================================")
  
  # 加载数据
  all_data <- load_preprocessed_data()
  
  # 数据集列表
  gse_ids <- names(all_data)
  
  # 存储结果
  all_dea_results <- list()
  
  # 处理每个数据集
  for (gse_id in gse_ids) {
    log_message(paste("\n========== 处理:", gse_id, "=========="))
    
    data <- all_data[[gse_id]]
    expr_matrix <- data$expression
    metadata <- data$metadata
    
    # 差异表达分析
    dea_results <- perform_dea(expr_matrix, metadata, paste0(gse_id, "_Normal_vs_Sepsis"))
    
    if (!is.null(dea_results)) {
      # MHC基因分析
      mhc_genes <- c("HLA-DRA", "HLA-DQB1", "HLA-DQA1", "HLA-DRB1", 
                     "HLA-DPA1", "HLA-DPB1", "CIITA", "CD74")
      mhc_dea <- analyze_mhc_dea(dea_results, mhc_genes)
      
      # 保存结果
      all_dea_results[[gse_id]] <- list(
        full_results = dea_results,
        mhc_results = mhc_dea
      )
      
      # 保存CSV
      write.csv(dea_results, 
                file.path("results", "differential_expression", 
                          paste0(gse_id, "_dea_results.csv")))
      
      # 可视化
      visualize_dea(dea_results, gse_id)
    }
  }
  
  # 合并所有数据集的差异基因
  log_message("\n========== 合并差异分析结果 ==========")
  
  combined_dea <- lapply(all_dea_results, function(x) x$full_results)
  
  # MHC基因跨数据集比较
  log_message("\n========== MHC基因跨数据集比较 ==========")
  
  for (gene in c("HLA-DRA", "HLA-DQB1", "HLA-DQA1")) {
    log_message(paste("\n基因:", gene))
    for (gse_id in names(combined_dea)) {
      gene_data <- combined_dea[[gse_id]]
      if (gene %in% gene_data$gene) {
        info <- gene_data[gene_data$gene == gene, ]
        log_message(paste("  ", gse_id, ":",
                           "logFC =", round(info$logFC, 3),
                           "P =", format(info$P.Value, scientific = TRUE),
                           "Status:", info$regulation))
      }
    }
  }
  
  # 保存所有结果
  saveRDS(all_dea_results, "results/differential_expression/all_dea_results.rds")
  
  log_message("\n========================================")
  log_message("差异表达分析完成!")
  log_message("========================================")
  
  return(all_dea_results)
}

# 执行
if (interactive()) {
  results <- main()
} else {
  results <- main()
}
