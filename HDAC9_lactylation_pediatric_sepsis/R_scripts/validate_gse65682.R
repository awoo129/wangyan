# =============================================================================
# GSE65682 - Install annotation package and do probe-to-gene mapping
# =============================================================================

library(GEOquery)
library(survival)

od <- "D:/下载/BMC_Submission_Package/revision_results"

# Install annotation package if needed
if (!requireNamespace("hgu219.db", quietly=TRUE)) {
  cat("Installing hgu219.db...\n")
  BiocManager::install("hgu219.db", ask=FALSE, update=FALSE)
}

library(hgu219.db)

# Load expression data
cat("Loading GSE65682 RDS...\n")
gse65682 <- readRDS(file.path(od, "GSE65682.rds"))
expr_mat <- exprs(gse65682)
pdata <- pData(gse65682)

cat("Expression dim:", dim(expr_mat), "\n")

# === Probe to gene mapping using hgu219.db ===
cat("\n=== Probe-to-Gene Mapping ===\n")

# Get all probe-to-symbol mappings
probe2symbol <- toTable(hgu219SYMBOL)
cat("Total probe-symbol mappings:", nrow(probe2symbol), "\n")
cat("Unique genes:", length(unique(probe2symbol$symbol)), "\n")

# Check HDAC9 probes
hdac9_probes <- probe2symbol$symbol == "HDAC9"
cat("\nHDAC9 probes in hgu219.db:\n")
print(probe2symbol[hdac9_probes, ])

# Check LARS genes
lars_genes <- c("HDAC9", "LDHA", "LDHB", "GAPDH")
for (g in lars_genes) {
  pg <- probe2symbol$symbol == g
  avail <- intersect(probe2symbol$probe_id[pg], rownames(expr_mat))
  cat(sprintf("\n%s: %d probes in DB, %d in expression matrix\n", g, sum(pg), length(avail)))
  if (length(avail) > 0) {
    cat("  Available probes:", paste(avail, collapse=", "), "\n")
  }
}

# === Build gene-level expression ===
cat("\n=== Building Gene-Level Expression ===\n")

# Map all probes to genes
mapped_probes <- probe2symbol$probe_id %in% rownames(expr_mat)
cat("Mapped probes:", sum(mapped_probes), "of", nrow(probe2symbol), "\n")

probe2symbol_mapped <- probe2symbol[mapped_probes, ]

# For genes with multiple probes, use IQR-based selection (same as training cohort)
unique_genes <- unique(probe2symbol_mapped$symbol)
cat("Unique genes:", length(unique_genes), "\n")

# Build gene-level matrix
expr_gene_list <- list()
for (gene in unique_genes) {
  probes <- probe2symbol_mapped$probe_id[probe2symbol_mapped$symbol == gene]
  probes <- intersect(probes, rownames(expr_mat))
  if (length(probes) == 0) next
  if (length(probes) == 1) {
    expr_gene_list[[gene]] <- as.numeric(expr_mat[probes, ])
  } else {
    # IQR-based selection: pick probe with largest IQR
    iqrs <- apply(expr_mat[probes, , drop=FALSE], 1, IQR, na.rm=TRUE)
    best <- probes[which.max(iqrs)]
    expr_gene_list[[gene]] <- as.numeric(expr_mat[best, ])
  }
}

expr_gene_mat <- do.call(rbind, expr_gene_list)
colnames(expr_gene_mat) <- colnames(expr_mat)
cat("Gene-level expression dim:", dim(expr_gene_mat), "\n")

# Verify LARS genes
for (g in lars_genes) {
  if (g %in% rownames(expr_gene_mat)) {
    vals <- expr_gene_mat[g, ]
    cat(sprintf("%s: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n", 
                g, mean(vals, na.rm=T), sd(vals, na.rm=T), 
                min(vals, na.rm=T), max(vals, na.rm=T)))
  } else {
    cat(g, ": NOT FOUND\n")
  }
}

# === Prepare survival data ===
cat("\n=== Preparing Survival Data ===\n")

mort_col <- "mortality_event_28days:ch1"
time_col <- "time_to_event_28days:ch1"

mort_clean <- gsub("mortality_event_28days: ", "", as.character(pdata[[mort_col]]))
time_clean <- gsub("time_to_event_28days: ", "", as.character(pdata[[time_col]]))

mort_num <- as.numeric(mort_clean)
time_num <- as.numeric(time_clean)

# ICU patients with survival data
icu_samples <- grepl("intensive-care", pdata$title, ignore.case=TRUE)
has_surv <- !is.na(mort_num) & !is.na(time_num)
valid <- icu_samples & has_surv

cat("Valid ICU patients with survival data:", sum(valid), "\n")
cat("Events:", sum(mort_num[valid]==1), "\n")

# Subset
expr_valid <- expr_gene_mat[, valid]
mort_valid <- mort_num[valid]
time_valid <- time_num[valid]

cat("Subset expression dim:", dim(expr_valid), "\n")

# === HDAC9 Cox Regression ===
cat("\n=== HDAC9 Cox Regression ===\n")

hdac9_expr <- as.numeric(expr_valid["HDAC9", ])
cox_hdac9 <- coxph(Surv(time_valid, mort_valid) ~ hdac9_expr)
cat("HDAC9 HR:", round(exp(coef(cox_hdac9)), 4), "\n")
cat("HDAC9 95%CI:", round(exp(confint(cox_hdac9)), 4), "\n")
cat("HDAC9 p-value:", round(summary(cox_hdac9)$coef[5], 6), "\n")

# C-index
if (!requireNamespace("survcomp", quietly=TRUE)) {
  cat("Installing survcomp...\n")
  BiocManager::install("survcomp", ask=FALSE, update=FALSE)
}
library(survcomp)
cidx_hdac9 <- concordance.index(hdac9_expr, time_valid, mort_valid, method="noether")
cat("HDAC9 C-index:", round(cidx_hdac9$c.index, 4), "\n")
cat("HDAC9 C-index se:", round(cidx_hdac9$se, 4), "\n")

# === HDAC9 Binary Stratification ===
cat("\n=== HDAC9 Binary Stratification ===\n")
hdac9_med <- median(hdac9_expr, na.rm=TRUE)
hdac9_binary <- ifelse(hdac9_expr > hdac9_med, "High", "Low")
cat("Median HDAC9:", round(hdac9_med, 3), "\n")
cat("High:", sum(hdac9_binary=="High"), "Low:", sum(hdac9_binary=="Low"), "\n")

km_fit <- survfit(Surv(time_valid, mort_valid) ~ hdac9_binary)
lr_test <- survdiff(Surv(time_valid, mort_valid) ~ hdac9_binary)
lr_p <- 1 - pchisq(lr_test$chisq, df=1)
cat("HDAC9 binary log-rank p:", round(lr_p, 6), "\n")

# C-index for binary
hdac9_bin_num <- ifelse(hdac9_binary=="High", 1, 0)
cidx_hdac9_bin <- concordance.index(hdac9_bin_num, time_valid, mort_valid, method="noether")
cat("HDAC9 binary C-index:", round(cidx_hdac9_bin$c.index, 4), "\n")

# === LARS Risk Score ===
cat("\n=== LARS Risk Score ===\n")

# Training-derived coefficients: HDAC9=-0.42, LDHA=0.48, LDHB=0.40, GAPDH=-0.11
lars_coef <- c(HDAC9=-0.42, LDHA=0.48, LDHB=0.40, GAPDH=-0.11)
lars_avail <- intersect(names(lars_coef), rownames(expr_valid))
cat("LARS genes available:", paste(lars_avail, collapse=", "), "\n")

if (length(lars_avail) == 4) {
  lars_score <- as.numeric(t(lars_coef[lars_avail] %*% expr_valid[lars_avail, ]))
  
  # LARS Cox
  cox_lars <- coxph(Surv(time_valid, mort_valid) ~ lars_score)
  cat("LARS HR:", round(exp(coef(cox_lars)), 4), "\n")
  cat("LARS 95%CI:", round(exp(confint(cox_lars)), 4), "\n")
  cat("LARS p-value:", round(summary(cox_lars)$coef[5], 6), "\n")
  
  # LARS C-index
  cidx_lars <- concordance.index(lars_score, time_valid, mort_valid, method="noether")
  cat("LARS C-index:", round(cidx_lars$c.index, 4), "\n")
  cat("LARS C-index se:", round(cidx_lars$se, 4), "\n")
  
  # LARS binary
  lars_med <- median(lars_score, na.rm=TRUE)
  lars_binary <- ifelse(lars_score > lars_med, "High", "Low")
  cat("LARS High:", sum(lars_binary=="High"), "Low:", sum(lars_binary=="Low"), "\n")
  
  km_lars <- survfit(Surv(time_valid, mort_valid) ~ lars_binary)
  lr_lars <- survdiff(Surv(time_valid, mort_valid) ~ lars_binary)
  lr_lars_p <- 1 - pchisq(lr_lars$chisq, df=1)
  cat("LARS binary log-rank p:", round(lr_lars_p, 6), "\n")
  
  lars_bin_num <- ifelse(lars_binary=="High", 1, 0)
  cidx_lars_bin <- concordance.index(lars_bin_num, time_valid, mort_valid, method="noether")
  cat("LARS binary C-index:", round(cidx_lars_bin$c.index, 4), "\n")
} else {
  cat("WARNING: Not all LARS genes available!\n")
}

# === Three-cohort Meta-analysis ===
cat("\n=== Three-Cohort Meta-analysis (HDAC9) ===\n")

# GSE66099: HR=0.659, 95%CI: 0.479-0.907, p=0.011
# GSE26440: HR=0.366, 95%CI: 0.044-3.027, p=0.351
# GSE65682: (from above)
# Fixed-effects meta-analysis using log(HR) and SE

# GSE66099
hr1 <- 0.659; lo1 <- 0.479; hi1 <- 0.907
se1 <- (log(hi1) - log(lo1)) / (2 * 1.96)

# GSE26440
hr2 <- 0.366; lo2 <- 0.044; hi2 <- 3.027
se2 <- (log(hi2) - log(lo2)) / (2 * 1.96)

# GSE65682
hr3 <- as.numeric(exp(coef(cox_hdac9)))
lo3 <- as.numeric(exp(confint(cox_hdac9)[1]))
hi3 <- as.numeric(exp(confint(cox_hdac9)[2]))
se3 <- (log(hi3) - log(lo3)) / (2 * 1.96)

cat(sprintf("GSE65682: HR=%.3f, 95%%CI: %.3f-%.3f, SE=%.4f\n", hr3, lo3, hi3, se3))

# Fixed-effects inverse-variance meta-analysis
log_hrs <- c(log(hr1), log(hr2), log(hr3))
ses <- c(se1, se2, se3)
weights <- 1 / ses^2
pooled_loghr <- sum(weights * log_hrs) / sum(weights)
pooled_se <- sqrt(1 / sum(weights))
pooled_hr <- exp(pooled_loghr)
pooled_lo <- exp(pooled_loghr - 1.96 * pooled_se)
pooled_hi <- exp(pooled_loghr + 1.96 * pooled_se)
pooled_z <- pooled_loghr / pooled_se
pooled_p <- 2 * pnorm(-abs(pooled_z))

# Heterogeneity
Q <- sum(weights * (log_hrs - pooled_loghr)^2)
df_Q <- length(log_hrs) - 1
Q_p <- 1 - pchisq(Q, df_Q)
I2 <- max(0, (Q - df_Q) / Q * 100)

cat(sprintf("\nPooled HR (3-cohort): %.3f (95%%CI: %.3f-%.3f)\n", pooled_hr, pooled_lo, pooled_hi))
cat(sprintf("Pooled p: %.6f\n", pooled_p))
cat(sprintf("Cochran Q: %.3f, df=%d, p=%.4f\n", Q, df_Q, Q_p))
cat(sprintf("I²: %.1f%%\n", I2))

# Also do 2-cohort (GSE66099 + GSE65682) excluding GSE26440
log_hrs_2 <- c(log(hr1), log(hr3))
ses_2 <- c(se1, se3)
weights_2 <- 1 / ses_2^2
pooled_loghr_2 <- sum(weights_2 * log_hrs_2) / sum(weights_2)
pooled_se_2 <- sqrt(1 / sum(weights_2))
pooled_hr_2 <- exp(pooled_loghr_2)
pooled_lo_2 <- exp(pooled_loghr_2 - 1.96 * pooled_se_2)
pooled_hi_2 <- exp(pooled_loghr_2 + 1.96 * pooled_se_2)
pooled_z_2 <- pooled_loghr_2 / pooled_se_2
pooled_p_2 <- 2 * pnorm(-abs(pooled_z_2))

Q_2 <- sum(weights_2 * (log_hrs_2 - pooled_loghr_2)^2)
df_Q_2 <- 1
Q_p_2 <- 1 - pchisq(Q_2, df_Q_2)
I2_2 <- max(0, (Q_2 - df_Q_2) / Q_2 * 100)

cat(sprintf("\n2-cohort (GSE66099+GSE65682) Pooled HR: %.3f (95%%CI: %.3f-%.3f)\n", 
            pooled_hr_2, pooled_lo_2, pooled_hi_2))
cat(sprintf("2-cohort p: %.6f\n", pooled_p_2))
cat(sprintf("2-cohort I²: %.1f%%\n", I2_2))

cat("\n========== GSE65682 VALIDATION COMPLETE ==========\n")
