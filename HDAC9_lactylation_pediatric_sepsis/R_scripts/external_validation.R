# =============================================================================
# BMC Bioinformatics Revision - External Validation Pipeline
# =============================================================================
# 验证队列：GSE26440 (GPL570, n=130, 有生存数据!), GSE25504 (Illumina, n=20, 无生存数据)
# GSE65682 (Affymetrix, 大队列) - 需从GEO下载
# 训练队列：GSE66099 (n=131, 已有结果)
#
# 分析内容：
#   1. 跨平台HDAC9表达验证（所有队列）
#   2. GSE26440外部验证：HDAC9 Cox回归 + LARS C-index
#   3. Forest plot元分析数据
# =============================================================================

sink("C:/Users/Administrator/external_validation_log.txt", split=TRUE)
cat("=== External Validation Pipeline started:", date(), "===\n")

tryCatch({
suppressPackageStartupMessages({
  library(Biobase)
  library(GEOquery)
  library(survival)
  library(survminer)
  library(glmnet)
  library(limma)
  # library(meta)  # not installed, using manual meta-analysis
  library(ggplot2)
})

od <- "D:/下载/BMC_Submission_Package/revision_results"
dir.create(od, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# PART 0: Load GPL570 annotation (shared by GSE66099 and GSE26440)
# =============================================================================
cat("\n--- Loading GPL570 annotation ---\n")
gpl570 <- getGEO(filename="D:/下载/GPL570.annot.gz")
gpl_table <- Table(gpl570)
probe_gene <- data.frame(probe=gpl_table$ID, gene=gpl_table[["Gene symbol"]], stringsAsFactors=FALSE)
probe_gene <- probe_gene[!is.na(probe_gene$gene) & probe_gene$gene!="" & !grepl("///", probe_gene$gene), ]
cat("GPL570 probes with gene symbols:", nrow(probe_gene), "\n")
rm(gpl570, gpl_table); gc()

# Function: probe-to-gene with IQR selection
probe_to_gene <- function(expr_matrix, probe_gene_map) {
  expr_filtered <- expr_matrix[rownames(expr_matrix) %in% probe_gene_map$probe, ]
  pg <- probe_gene_map[match(rownames(expr_filtered), probe_gene_map$probe), ]
  iqrs <- apply(expr_filtered, 1, IQR)
  pg$iqr <- iqrs
  pg <- pg[order(pg$gene, -pg$iqr), ]
  pg <- pg[!duplicated(pg$gene), ]
  expr_gene <- expr_filtered[pg$probe, ]
  rownames(expr_gene) <- pg$gene
  return(expr_gene)
}

# LARS genes and coefficients from training
lars_genes <- c("HDAC9", "LDHA", "LDHB", "GAPDH")
lars_coefs <- c(HDAC9=-0.42, LDHA=0.48, LDHB=0.40, GAPDH=-0.11)

# =============================================================================
# PART 1: Load Training Cohort (GSE66099) - reload for meta-analysis
# =============================================================================
cat("\n============================================================\n")
cat("PART 1: Loading Training Cohort GSE66099\n")
cat("============================================================\n")

gse66099 <- readRDS(file.path(od, "GSE66099.rds"))
expr_66099 <- exprs(gse66099)
expr_66099_gene <- probe_to_gene(expr_66099, probe_gene)
cat("GSE66099 gene-level dim:", dim(expr_66099_gene), "\n")

# Load survival data
surv_66099 <- read.csv(file.path(od, "survival_data_real.csv"), stringsAsFactors=FALSE)
cat("GSE66099 survival data:", nrow(surv_66099), "samples,", sum(surv_66099$event), "events\n")

# Verify HDAC9 results
if ("HDAC9" %in% rownames(expr_66099_gene)) {
  hdac9_66099 <- as.numeric(expr_66099_gene["HDAC9", surv_66099$gsm])
  cox_66099 <- coxph(Surv(time, event) ~ hdac9_66099, data=surv_66099)
  cat("GSE66099 HDAC9 HR:", round(exp(coef(cox_66099)), 3), 
      "95%CI:", round(exp(confint(cox_66099))[1], 3), "-", round(exp(confint(cox_66099))[2], 3),
      "p=", signif(summary(cox_66099)$coef[5], 3), 
      "C-index:", round(concordance(cox_66099)$concordance, 3), "\n")
}
rm(expr_66099, gse66099); gc()

# =============================================================================
# PART 2: GSE26440 - External Validation (GPL570, HAS SURVIVAL DATA!)
# =============================================================================
cat("\n============================================================\n")
cat("PART 2: GSE26440 External Validation\n")
cat("============================================================\n")

gse26440 <- readRDS(file.path(od, "GSE26440.rds"))
expr_26440 <- exprs(gse26440)
pdata_26440 <- pData(gse26440)
cat("GSE26440 raw dim:", dim(expr_26440), "\n")

# Extract outcome data
outcome_col <- "outcome:ch1"
if (outcome_col %in% colnames(pdata_26440)) {
  cat("\nOutcome column found!\n")
  cat("Outcome values:\n")
  print(table(pdata_26440[[outcome_col]]))
  
  # Extract disease state
  disease_col <- "disease state:ch1"
  cat("\nDisease state:\n")
  print(table(pdata_26440[[disease_col]]))
  
  # Filter: septic shock patients with outcome data
  is_septic <- grepl("septic", pdata_26440[[disease_col]], ignore.case=TRUE)
  has_outcome <- pdata_26440[[outcome_col]] %in% c("Nonsurvivor", "Survivor")
  use_idx <- is_septic & has_outcome
  
  surv_26440 <- data.frame(
    gsm=rownames(pdata_26440)[use_idx],
    outcome=pdata_26440[[outcome_col]][use_idx],
    stringsAsFactors=FALSE
  )
  surv_26440$event <- ifelse(surv_26440$outcome == "Nonsurvivor", 1, 0)
  surv_26440$time <- ifelse(surv_26440$event == 1, 14, 28)  # same assumption as training
  
  cat("\nGSE26440 septic shock patients with outcome:", nrow(surv_26440), "\n")
  cat("Events:", sum(surv_26440$event), "Non-events:", sum(surv_26440$event == 0), "\n")
  cat("Event rate:", round(mean(surv_26440$event) * 100, 1), "%\n")
  
  # Process expression data
  expr_26440_gene <- probe_to_gene(expr_26440, probe_gene)
  cat("GSE26440 gene-level dim:", dim(expr_26440_gene), "\n")
  
  # Match samples
  gsm_in_expr <- intersect(surv_26440$gsm, colnames(expr_26440_gene))
  surv_26440 <- surv_26440[surv_26440$gsm %in% gsm_in_expr, ]
  expr_26440_sub <- expr_26440_gene[, surv_26440$gsm]
  cat("Matched samples:", nrow(surv_26440), "\n")
  
  if (nrow(surv_26440) > 0 && sum(surv_26440$event) >= 5) {
    # --- HDAC9 Cox regression ---
    cat("\n--- HDAC9 External Validation (GSE26440) ---\n")
    if ("HDAC9" %in% rownames(expr_26440_sub)) {
      hdac9_26440 <- as.numeric(expr_26440_sub["HDAC9", ])
      surv_26440$HDAC9 <- hdac9_26440
      
      cox_hdac9_26440 <- coxph(Surv(time, event) ~ HDAC9, data=surv_26440)
      hr_26440 <- exp(coef(cox_hdac9_26440))
      ci_26440 <- exp(confint(cox_hdac9_26440))
      p_26440 <- summary(cox_hdac9_26440)$coef[5]
      c_26440 <- concordance(cox_hdac9_26440)$concordance
      
      cat("GSE26440 HDAC9 HR:", round(hr_26440, 3), 
          "95%CI:", round(ci_26440[1], 3), "-", round(ci_26440[2], 3),
          "p=", signif(p_26440, 3), "C-index:", round(c_26440, 3), "\n")
      
      # HDAC9 binary
      surv_26440$HDAC9_binary <- ifelse(hdac9_26440 > median(hdac9_26440), "High", "Low")
      fit_km_26440 <- survfit(Surv(time, event) ~ HDAC9_binary, data=surv_26440)
      pval_km_26440 <- surv_pvalue(fit_km_26440)$pval
      cat("HDAC9 binary KM p=", signif(pval_km_26440, 3), "\n")
      
      # KM plot
      pdf(file.path(od, "KM_HDAC9_GSE26440.pdf"), width=8, height=6)
      p_km <- ggsurvplot(fit_km_26440, data=surv_26440, pval=TRUE,
        title=sprintf("GSE26440 HDAC9 (C=%.3f, HR=%.3f)", c_26440, hr_26440),
        risk.table=TRUE, palette="jco", xlab="Days", ylab="Survival")
      print(p_km); dev.off()
    } else {
      cat("HDAC9 not found in GSE26440 expression matrix!\n")
    }
    
    # --- LARS External Validation ---
    cat("\n--- LARS External Validation (GSE26440) ---\n")
    lars_available <- intersect(lars_genes, rownames(expr_26440_sub))
    cat("LARS genes available in GSE26440:", paste(lars_available, collapse=", "), "\n")
    
    if (length(lars_available) == 4) {
      # Use training-derived coefficients
      x_lars_26440 <- t(expr_26440_sub[lars_available, ])
      lars_score_26440 <- as.numeric(x_lars_26440 %*% lars_coefs[lars_available])
      surv_26440$LARS <- lars_score_26440
      
      # LARS continuous C-index
      cox_lars_26440 <- coxph(Surv(time, event) ~ LARS, data=surv_26440)
      lars_c_26440 <- concordance(cox_lars_26440)$concordance
      lars_hr_26440 <- exp(coef(cox_lars_26440))
      lars_ci_26440 <- exp(confint(cox_lars_26440))
      lars_p_26440 <- summary(cox_lars_26440)$coef[5]
      
      cat("GSE26440 LARS C-index:", round(lars_c_26440, 3), 
          "HR:", round(lars_hr_26440, 3),
          "95%CI:", round(lars_ci_26440[1], 3), "-", round(lars_ci_26440[2], 3),
          "p=", signif(lars_p_26440, 3), "\n")
      
      # LARS binary
      surv_26440$LARS_binary <- ifelse(lars_score_26440 > median(lars_score_26440), "High", "Low")
      fit_lars_26440 <- survfit(Surv(time, event) ~ LARS_binary, data=surv_26440)
      pval_lars_26440 <- surv_pvalue(fit_lars_26440)$pval
      cat("LARS binary KM p=", signif(pval_lars_26440, 3), "\n")
      
      cox_lars_b_26440 <- coxph(Surv(time, event) ~ LARS_binary, data=surv_26440)
      lars_b_c_26440 <- concordance(cox_lars_b_26440)$concordance
      cat("LARS binary C-index:", round(lars_b_c_26440, 3), "\n")
      
      # KM plot
      pdf(file.path(od, "KM_LARS_GSE26440.pdf"), width=8, height=6)
      p_lars <- ggsurvplot(fit_lars_26440, data=surv_26440, pval=TRUE,
        title=sprintf("GSE26440 LARS (C=%.3f)", lars_c_26440),
        risk.table=TRUE, palette="jco", xlab="Days", ylab="Survival")
      print(p_lars); dev.off()
    } else {
      cat("Not all LARS genes available, computing with available genes\n")
      x_lars_26440 <- t(expr_26440_sub[lars_available, ])
      lars_score_26440 <- as.numeric(x_lars_26440 %*% lars_coefs[lars_available])
      surv_26440$LARS <- lars_score_26440
      
      cox_lars_26440 <- coxph(Surv(time, event) ~ LARS, data=surv_26440)
      lars_c_26440 <- concordance(cox_lars_26440)$concordance
      cat("GSE26440 LARS (partial) C-index:", round(lars_c_26440, 3), "\n")
    }
    
    # Save GSE26440 survival data
    write.csv(surv_26440, file.path(od, "survival_data_GSE26440.csv"), row.names=FALSE)
    cat("Saved GSE26440 survival data\n")
  } else {
    cat("Insufficient events for survival analysis in GSE26440\n")
  }
} else {
  cat("No outcome column found in GSE26440!\n")
}

rm(gse26440, expr_26440); gc()

# =============================================================================
# PART 3: GSE25504 - Expression-level Validation (Illumina, no survival data)
# =============================================================================
cat("\n============================================================\n")
cat("PART 3: GSE25504 Expression-level Validation\n")
cat("============================================================\n")

gse25504 <- readRDS(file.path(od, "GSE25504.rds"))
expr_25504 <- exprs(gse25504)
pdata_25504 <- pData(gse25504)
cat("GSE25504 dim:", dim(expr_25504), "\n")
cat("Platform:", annotation(gse25504), "\n")

# GSE25504 is Illumina GPL13667 - probe annotation is in the featureData
# Check if featureData has gene symbols
fdata_25504 <- fData(gse25504)
cat("Feature columns:\n")
print(head(colnames(fdata_25504)))

# Try to find gene symbol column
gene_col_25504 <- NULL
for (col in c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol", "ILMN_GENE", 
              "TargetID", "SYMBOL", "gene")) {
  if (col %in% colnames(fdata_25504)) {
    gene_col_25504 <- col
    break
  }
}
cat("Gene column found:", gene_col_25504, "\n")

if (!is.null(gene_col_25504)) {
  pg_25504 <- data.frame(probe=rownames(fdata_25504), gene=fdata_25504[[gene_col_25504]], stringsAsFactors=FALSE)
  pg_25504 <- pg_25504[!is.na(pg_25504$gene) & pg_25504$gene!="" & !grepl("///", pg_25504$gene), ]
  
  expr_25504_gene <- probe_to_gene(expr_25504, pg_25504)
  cat("GSE25504 gene-level dim:", dim(expr_25504_gene), "\n")
  
  # Check HDAC9 expression
  if ("HDAC9" %in% rownames(expr_25504_gene)) {
    hdac9_25504 <- as.numeric(expr_25504_gene["HDAC9", ])
    cat("GSE25504 HDAC9 expression range:", round(min(hdac9_25504), 3), "-", round(max(hdac9_25504), 3), "\n")
    cat("GSE25504 HDAC9 mean:", round(mean(hdac9_25504), 3), "\n")
  }
  
  # Check LARS genes
  lars_avail_25504 <- intersect(lars_genes, rownames(expr_25504_gene))
  cat("LARS genes in GSE25504:", paste(lars_avail_25504, collapse=", "), "\n")
} else {
  cat("No gene annotation found for GSE25504 - checking rownames pattern\n")
  # Illumina probe IDs - try to match by pattern
  cat("First few rownames:\n")
  print(head(rownames(expr_25504), 10))
}

rm(gse25504, expr_25504); gc()

# =============================================================================
# PART 4: GSE65682 - Load from GEO (large Affymetrix dataset)
# =============================================================================
cat("\n============================================================\n")
cat("PART 4: GSE65682 External Validation\n")
cat("============================================================\n")

# Try loading from GEO with GSEMatrix=TRUE
cat("Loading GSE65682 from GEO (this may take several minutes)...\n")
gse65682 <- tryCatch({
  getGEO("GSE65682", GSEMatrix=TRUE, getGPL=FALSE)
}, error = function(e) {
  cat("Error loading GSE65682 from GEO:", e$message, "\n")
  cat("Trying local soft file...\n")
  soft_path <- "D:/下载/GSE65682_family.soft.gz"
  if (file.exists(soft_path)) {
    # Load from soft file and extract GSM objects
    gse_obj <- getGEO(filename=soft_path, getGPL=FALSE)
    # Convert GSE to list of ESets
    if (class(gse_obj) == "GSE") {
      # Extract GSMs from GSE object
      gsms <- GSMList(gse_obj)
      cat("GSMs found:", length(gsms), "\n")
      # Create ExpressionSet from GSMs
      # This is complex - need to build expression matrix
      cat("GSE object - need different approach\n")
      return(NULL)
    }
    return(gse_obj)
  }
  return(NULL)
})

if (!is.null(gse65682)) {
  if (is.list(gse65682)) {
    cat("Multiple ESets found:", length(gse65682), "\n")
    for (i in seq_along(gse65682)) {
      g <- gse65682[[i]]
      cat(sprintf("  ESet %d: %s, dim=%s\n", i, annotation(g), 
                  paste(dim(exprs(g)), collapse=" x ")))
    }
    # Use the largest
    sizes <- sapply(gse65682, function(g) ncol(exprs(g)))
    gse65682 <- gse65682[[which.max(sizes)]]
  }
  
  expr_65682 <- exprs(gse65682)
  pdata_65682 <- pData(gse65682)
  cat("GSE65682 dim:", dim(expr_65682), "\n")
  cat("Platform:", annotation(gse65682), "\n")
  
  # Check for outcome data
  cat("\nGSE65682 phenotype columns:\n")
  print(colnames(pdata_65682))
  
  # Look for outcome/survival columns
  outcome_cols <- grep("outcome|mortality|survival|death|status|characteristics", 
                       colnames(pdata_65682), value=TRUE, ignore.case=TRUE)
  for (col in outcome_cols) {
    cat(sprintf("\n--- %s ---\n", col))
    cat("Unique values (first 20):\n")
    print(head(unique(pdata_65682[[col]]), 20))
  }
  
  # Save for later use
  saveRDS(gse65682, file.path(od, "GSE65682.rds"))
  cat("Saved GSE65682.rds\n")
  
  # Check if survival data available
  has_survival <- FALSE
  outcome_col_65682 <- NULL
  
  # Try different column patterns
  for (col in colnames(pdata_65682)) {
    vals <- pdata_65682[[col]]
    if (is.character(vals)) {
      if (any(grepl("Nonsurvivor|Non-survivor|nonsurvivor|Dead|dead|Death|death", vals))) {
        outcome_col_65682 <- col
        has_survival <- TRUE
        cat("\nFound survival column:", col, "\n")
        cat("Values:\n")
        print(table(vals))
        break
      }
    }
  }
  
  if (has_survival && !is.null(outcome_col_65682)) {
    cat("\n--- GSE65682 has survival data! Running validation ---\n")
    
    # Process expression data
    expr_65682_gene <- probe_to_gene(expr_65682, probe_gene)
    cat("GSE65682 gene-level dim:", dim(expr_65682_gene), "\n")
    
    # Extract survival data
    # Filter samples with outcome
    has_outcome_65682 <- !is.na(pdata_65682[[outcome_col_65682]]) & 
                          pdata_65682[[outcome_col_65682]] != ""
    
    surv_65682 <- data.frame(
      gsm=rownames(pdata_65682)[has_outcome_65682],
      outcome=pdata_65682[[outcome_col_65682]][has_outcome_65682],
      stringsAsFactors=FALSE
    )
    
    # Standardize outcome
    surv_65682$event <- ifelse(grepl("Non|Dead|Death|dead|death", surv_65682$outcome), 1, 0)
    surv_65682$time <- ifelse(surv_65682$event == 1, 14, 28)
    
    # Filter disease samples (exclude controls if any)
    # Check disease state column
    disease_col_65682 <- grep("disease|source|group|type", colnames(pdata_65682), 
                               value=TRUE, ignore.case=TRUE)[1]
    if (!is.na(disease_col_65682) && disease_col_65682 %in% colnames(pdata_65682)) {
      cat("Disease column:", disease_col_65682, "\n")
      print(table(pdata_65682[[disease_col_65682]]))
    }
    
    cat("GSE65682 total samples with outcome:", nrow(surv_65682), "\n")
    cat("Events:", sum(surv_65682$event), "\n")
    
    # Match with expression data
    gsm_in_expr <- intersect(surv_65682$gsm, colnames(expr_65682_gene))
    surv_65682 <- surv_65682[surv_65682$gsm %in% gsm_in_expr, ]
    expr_65682_sub <- expr_65682_gene[, surv_65682$gsm]
    cat("Matched samples:", nrow(surv_65682), "\n")
    
    if (nrow(surv_65682) > 0 && sum(surv_65682$event) >= 5) {
      # HDAC9
      if ("HDAC9" %in% rownames(expr_65682_sub)) {
        hdac9_65682 <- as.numeric(expr_65682_sub["HDAC9", ])
        surv_65682$HDAC9 <- hdac9_65682
        
        cox_65682 <- coxph(Surv(time, event) ~ HDAC9, data=surv_65682)
        cat("GSE65682 HDAC9 HR:", round(exp(coef(cox_65682)), 3), 
            "95%CI:", round(exp(confint(cox_65682))[1], 3), "-", round(exp(confint(cox_65682))[2], 3),
            "p=", signif(summary(cox_65682)$coef[5], 3),
            "C-index:", round(concordance(cox_65682)$concordance, 3), "\n")
      }
      
      # LARS
      lars_avail_65682 <- intersect(lars_genes, rownames(expr_65682_sub))
      cat("LARS genes in GSE65682:", paste(lars_avail_65682, collapse=", "), "\n")
      if (length(lars_avail_65682) == 4) {
        x_lars_65682 <- t(expr_65682_sub[lars_avail_65682, ])
        lars_score_65682 <- as.numeric(x_lars_65682 %*% lars_coefs[lars_avail_65682])
        surv_65682$LARS <- lars_score_65682
        
        cox_lars_65682 <- coxph(Surv(time, event) ~ LARS, data=surv_65682)
        cat("GSE65682 LARS C-index:", round(concordance(cox_lars_65682)$concordance, 3), "\n")
      }
      
      write.csv(surv_65682, file.path(od, "survival_data_GSE65682.csv"), row.names=FALSE)
    }
  } else {
    cat("\nGSE65682 does NOT have survival data. Expression-level validation only.\n")
    
    # Process expression data for HDAC9 downregulation check
    expr_65682_gene <- probe_to_gene(expr_65682, probe_gene)
    cat("GSE65682 gene-level dim:", dim(expr_65682_gene), "\n")
    
    if ("HDAC9" %in% rownames(expr_65682_gene)) {
      hdac9_65682_all <- as.numeric(expr_65682_gene["HDAC9", ])
      cat("GSE65682 HDAC9 mean:", round(mean(hdac9_65682_all, na.rm=TRUE), 3), 
          "SD:", round(sd(hdac9_65682_all, na.rm=TRUE), 3), "\n")
    }
  }
  
  rm(gse65682, expr_65682); gc()
} else {
  cat("GSE65682 could not be loaded. Skipping.\n")
}

# =============================================================================
# PART 5: Cross-platform HDAC9 expression comparison
# =============================================================================
cat("\n============================================================\n")
cat("PART 5: Cross-platform HDAC9 Expression Comparison\n")
cat("============================================================\n")

# Collect HDAC9 expression from all cohorts
hdac9_comparison <- data.frame()

# GSE66099 (training)
if (exists("surv_66099") && "HDAC9" %in% colnames(surv_66099)) {
  hdac9_comparison <- rbind(hdac9_comparison, data.frame(
    Cohort="GSE66099", Platform="Affymetrix_HGU133Plus2",
    N=nrow(surv_66099), N_events=sum(surv_66099$event),
    HDAC9_mean=round(mean(surv_66099$HDAC9, na.rm=TRUE), 3),
    HDAC9_sd=round(sd(surv_66099$HDAC9, na.rm=TRUE), 3),
    Has_survival=TRUE
  ))
}

# GSE26440
if (exists("surv_26440") && "HDAC9" %in% colnames(surv_26440)) {
  hdac9_comparison <- rbind(hdac9_comparison, data.frame(
    Cohort="GSE26440", Platform="Affymetrix_HGU133Plus2",
    N=nrow(surv_26440), N_events=sum(surv_26440$event),
    HDAC9_mean=round(mean(surv_26440$HDAC9, na.rm=TRUE), 3),
    HDAC9_sd=round(sd(surv_26440$HDAC9, na.rm=TRUE), 3),
    Has_survival=TRUE
  ))
}

cat("\nCross-platform HDAC9 comparison:\n")
print(hdac9_comparison)
write.csv(hdac9_comparison, file.path(od, "HDAC9_cross_platform.csv"), row.names=FALSE)

# =============================================================================
# PART 6: Forest Plot Data (Meta-analysis)
# =============================================================================
cat("\n============================================================\n")
cat("PART 6: Forest Plot Data\n")
cat("============================================================\n")

forest_data <- data.frame(
  Cohort=character(), N=integer(), Events=integer(),
  HR=numeric(), HR_lower=numeric(), HR_upper=numeric(),
  p_value=numeric(), C_index=numeric(),
  stringsAsFactors=FALSE
)

# GSE66099 (training)
if (exists("surv_66099") && "HDAC9" %in% rownames(expr_66099_gene)) {
  # Re-run since we may not have the cox object
  hdac9_tr <- as.numeric(expr_66099_gene["HDAC9", surv_66099$gsm])
  cox_tr <- coxph(Surv(time, event) ~ hdac9_tr, data=surv_66099)
  forest_data <- rbind(forest_data, data.frame(
    Cohort="GSE66099 (Training)", N=nrow(surv_66099), Events=sum(surv_66099$event),
    HR=round(exp(coef(cox_tr)), 3),
    HR_lower=round(exp(confint(cox_tr))[1], 3),
    HR_upper=round(exp(confint(cox_tr))[2], 3),
    p_value=signif(summary(cox_tr)$coef[5], 3),
    C_index=round(concordance(cox_tr)$concordance, 3),
    stringsAsFactors=FALSE
  ))
}

# GSE26440 (external validation)
if (exists("surv_26440") && "HDAC9" %in% colnames(surv_26440)) {
  cox_26440_f <- coxph(Surv(time, event) ~ HDAC9, data=surv_26440)
  forest_data <- rbind(forest_data, data.frame(
    Cohort="GSE26440 (Validation)", N=nrow(surv_26440), Events=sum(surv_26440$event),
    HR=round(exp(coef(cox_26440_f)), 3),
    HR_lower=round(exp(confint(cox_26440_f))[1], 3),
    HR_upper=round(exp(confint(cox_26440_f))[2], 3),
    p_value=signif(summary(cox_26440_f)$coef[5], 3),
    C_index=round(concordance(cox_26440_f)$concordance, 3),
    stringsAsFactors=FALSE
  ))
}

cat("\nForest plot data:\n")
print(forest_data)

# Meta-analysis using inverse variance method
if (nrow(forest_data) >= 2) {
  cat("\n--- Meta-analysis ---\n")
  
  # Calculate log HR and SE for each study
  meta_results <- data.frame()
  for (i in 1:nrow(forest_data)) {
    log_hr <- log(forest_data$HR[i])
    # SE from CI
    se <- (log(forest_data$HR_upper[i]) - log(forest_data$HR_lower[i])) / (2 * 1.96)
    meta_results <- rbind(meta_results, data.frame(
      study=forest_data$Cohort[i],
      log_hr=log_hr, se=se,
      weight=1/se^2,
      n=forest_data$N[i],
      events=forest_data$Events[i]
    ))
  }
  
  # Fixed-effects pooled estimate
  total_weight <- sum(meta_results$weight)
  pooled_log_hr <- sum(meta_results$log_hr * meta_results$weight) / total_weight
  pooled_se <- sqrt(1 / total_weight)
  pooled_hr <- exp(pooled_log_hr)
  pooled_lower <- exp(pooled_log_hr - 1.96 * pooled_se)
  pooled_upper <- exp(pooled_log_hr + 1.96 * pooled_se)
  pooled_p <- 2 * pnorm(-abs(pooled_log_hr / pooled_se))
  
  # Heterogeneity (Cochran's Q)
  Q <- sum(meta_results$weight * (meta_results$log_hr - pooled_log_hr)^2)
  df_Q <- nrow(meta_results) - 1
  I2 <- max(0, (Q - df_Q) / Q * 100)
  Q_p <- 1 - pchisq(Q, df_Q)
  
  cat("\nPooled HR:", round(pooled_hr, 3), 
      "95%CI:", round(pooled_lower, 3), "-", round(pooled_upper, 3),
      "p=", signif(pooled_p, 3), "\n")
  cat("Heterogeneity: Q=", round(Q, 3), "df=", df_Q, 
      "p=", signif(Q_p, 3), "I²=", round(I2, 1), "%\n")
  
  # Add pooled estimate to forest data
  forest_data <- rbind(forest_data, data.frame(
    Cohort="Pooled (Fixed-effects)",
    N=sum(forest_data$N), Events=sum(forest_data$Events),
    HR=round(pooled_hr, 3),
    HR_lower=round(pooled_lower, 3),
    HR_upper=round(pooled_upper, 3),
    p_value=signif(pooled_p, 3),
    C_index=NA,
    stringsAsFactors=FALSE
  ))
}

write.csv(forest_data, file.path(od, "forest_plot_data.csv"), row.names=FALSE)
cat("Saved forest plot data\n")

# =============================================================================
# PART 7: Model Comparison Summary (Training + External Validation)
# =============================================================================
cat("\n============================================================\n")
cat("PART 7: Model Comparison Summary\n")
cat("============================================================\n")

model_summary <- data.frame(
  Model=character(), Cohort=character(), N=integer(), Events=integer(),
  C_index=numeric(), HR=numeric(), HR_lower=numeric(), HR_upper=numeric(), 
  p_value=numeric(),
  stringsAsFactors=FALSE
)

# Add training results
model_summary <- rbind(model_summary, data.frame(
  Model="HDAC9 (continuous)", Cohort="GSE66099 (Training)", 
  N=131, Events=22, C_index=0.672,
  HR=0.659, HR_lower=0.479, HR_upper=0.907, p_value=0.011,
  stringsAsFactors=FALSE
))

model_summary <- rbind(model_summary, data.frame(
  Model="LARS (continuous)", Cohort="GSE66099 (Training)",
  N=131, Events=22, C_index=0.774,
  HR=NA, HR_lower=NA, HR_upper=NA, p_value=NA,
  stringsAsFactors=FALSE
))

# Add GSE26440 validation results if available
if (exists("surv_26440") && "HDAC9" %in% colnames(surv_26440)) {
  cox_h9_val <- coxph(Surv(time, event) ~ HDAC9, data=surv_26440)
  model_summary <- rbind(model_summary, data.frame(
    Model="HDAC9 (continuous)", Cohort="GSE26440 (Validation)",
    N=nrow(surv_26440), Events=sum(surv_26440$event), 
    C_index=round(concordance(cox_h9_val)$concordance, 3),
    HR=round(exp(coef(cox_h9_val)), 3),
    HR_lower=round(exp(confint(cox_h9_val))[1], 3),
    HR_upper=round(exp(confint(cox_h9_val))[2], 3),
    p_value=signif(summary(cox_h9_val)$coef[5], 3),
    stringsAsFactors=FALSE
  ))
  
  if ("LARS" %in% colnames(surv_26440)) {
    cox_lars_val <- coxph(Surv(time, event) ~ LARS, data=surv_26440)
    model_summary <- rbind(model_summary, data.frame(
      Model="LARS (continuous)", Cohort="GSE26440 (Validation)",
      N=nrow(surv_26440), Events=sum(surv_26440$event),
      C_index=round(concordance(cox_lars_val)$concordance, 3),
      HR=round(exp(coef(cox_lars_val)), 3),
      HR_lower=round(exp(confint(cox_lars_val))[1], 3),
      HR_upper=round(exp(confint(cox_lars_val))[2], 3),
      p_value=signif(summary(cox_lars_val)$coef[5], 3),
      stringsAsFactors=FALSE
    ))
  }
}

cat("\nModel Comparison Summary:\n")
print(model_summary)
write.csv(model_summary, file.path(od, "model_comparison_external.csv"), row.names=FALSE)

cat("\n========== ALL EXTERNAL VALIDATION COMPLETE ==========\n")

}, error = function(e) {
  cat("FATAL ERROR:", conditionMessage(e), "\n")
  cat(traceback(), "\n")
})

sink()
