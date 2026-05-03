# Figure Legends and Captions

---

## Figure 1. Study Design and Analytical Pipeline

![Study Flowchart](figures/fig1_study_flowchart.png)

**Figure 1. Study Design and Analytical Pipeline**

**(A)** Overview of the four-route analytical framework. Six independent cohorts were analyzed across pediatric (GSE26440, GSE26378, GSE13904) and adult (GSE65682 discovery, GSE65682 validation, GAinS) populations. **(B)** Route A: Pediatric IPS model construction using the discovery cohort (GSE26440, n=130) with 10-fold cross-validation and independent validation (GSE26378, n=103). **(C)** Route B: Adult IPS model using the MARS cohort (GSE65682, n=479) with discovery-validation split. **(D)** Route C: External validation of the pediatric signature in an independent cohort (GSE13904, n=139). **(E)** Route D: Mechanistic analysis of age-dependent CIITA expression patterns and IFN-γ-STAT1-IRF1 axis regulation.

**Sample sizes**: Pediatric discovery: n=130 (Subclass A: 28, B: 45, C: 25, NA: 32); Pediatric validation: n=103; Pediatric external: n=139; Adult discovery: n=99; Adult validation: n=380; Adult external (GAinS): n=106.

---

## Figure 2. CIITA Expression Comparison Between Pediatric and Adult Sepsis

![CIITA Expression Comparison](figures/fig2_ciita_comparison.png)

**Figure 2. Age-Dependent Discordance of CIITA Expression in Sepsis-Induced Immunoparalysis**

**(A)** Pediatric sepsis (GSE26440): CIITA expression was significantly downregulated in immunoparalysis (Subclass A, red) compared to non-paralysis (Subclass B+C, blue). Mean expression: 7.23 vs 7.81, log2FC = -0.576, p = 7.97×10⁻⁹. **(B)** Adult sepsis (GSE65682): CIITA expression was significantly upregulated in immunoparalysis (Mars1, red) compared to non-paralysis (Mars2/3/4, blue). Mean expression: 2.44 vs 2.32, log2FC = +0.122, p = 0.0034.

Box plots display median (center line), interquartile range (box), and 1.5×IQR (whiskers). Each point represents an individual sample. Statistical significance assessed by two-sided t-test with Benjamini-Hochberg FDR correction.

**Key finding**: CIITA demonstrates opposite directional changes in pediatric versus adult immunoparalysis.

---

## Figure 3. ROC Curves for Immunoparalysis Classification

![ROC Curves](figures/fig3_roc_curves.png)

**Figure 3. Diagnostic Performance of IPS Models for Immunoparalysis Classification**

**(A)** Pediatric IPS model (Route A): ROC curve demonstrating excellent discrimination between immunoparalysis and non-paralysis in the discovery cohort (GSE26440, n=98). Training AUC = 0.870 (95% CI: 0.804–0.936), optimal threshold = 0.31. **(B)** Confusion matrix showing model classification performance: 60 true negatives, 25 true positives, 10 false positives, and 3 false negatives.

**(C)** Adult IPS model (Route B): ROC curve for the MARS cohort discovery set (GSE65682, n=99). Training AUC = 0.830 (95% CI: 0.756–0.904). **(D)** Validation performance: AUC = 0.696 (95% CI: 0.641–0.751) in the validation cohort (n=380).

Models were constructed using L2-regularized logistic regression with 10-fold stratified cross-validation. External validation AUC = 0.550 in the GAinS cohort (n=106) due to platform effects.

---

## Figure 4. Mechanistic Hypothesis: Age-Dependent Immunoparalysis Patterns

![Mechanism Hypothesis](figures/fig4_mechanism_hypothesis.png)

**Figure 4. Proposed Mechanistic Framework for Age-Dependent CIITA Expression Discordance**

**(A)** Pediatric "Pure Suppression" Pattern: In children, sepsis induces coordinated downregulation of the entire IFN-γ-STAT1-IRF1-CIITA axis. IFN-γ production decreases, leading to reduced STAT1 phosphorylation and IRF1 transcription factor activity. This results in CIITA silencing and subsequent MHC class II molecule downregulation, creating a "complete shutdown" of antigen presentation capacity.

**(B)** Adult "Compensatory Inflammation" Pattern: In adults, despite MHC II protein downregulation (as measured by mHLA-DR flow cytometry), CIITA mRNA is paradoxically upregulated. This suggests compensatory activation of the IFN-γ pathway attempting to restore MHC II expression, but failing due to post-transcriptional or post-translational regulatory mechanisms (potential involvement of miR-155 or protein trafficking defects).

Arrows indicate direction of regulatory changes (↓ down, ↑ up). Dashed arrow indicates attempted but incomplete compensation.

---

## Figure 5. MHC Class II Gene Expression Heatmap

![MHC II Heatmap](figures/fig5_mhc_heatmap.png)

**Figure 5. Coordinated Downregulation of MHC Class II Antigen Presentation Pathway in Pediatric Sepsis**

Heatmap showing expression levels of 9 MHC class II-related genes across pediatric sepsis subgroups (GSE26440). Genes are arranged by functional category: CIITA (master regulator), MHC II alpha/beta chains (HLA-DRA, HLA-DQA1, HLA-DQB1, HLA-DRB1, HLA-DPB1), peptide loading complex (HLA-DMA, HLA-DMB), and invariant chain (CD74).

Samples are clustered by immunoparalysis status: Subclass A (immunoparalysis, red bar) versus Subclass B+C (non-paralysis, blue bar). Color scale represents z-score normalized expression (blue: low, white: intermediate, red: high).

The heatmap reveals coordinated downregulation of CIITA and downstream MHC II genes in immunoparalysis, with CIITA (log2FC = -0.576, FDR = 6.37×10⁻⁸) showing the most significant change.

---

## Figure S1. Volcano Plots Showing Differential Expression

![Volcano Plots](figures/figS1_volcano_plot.png)

**Figure S1. Differential Expression Analysis in Pediatric Sepsis Discovery Cohort**

**(A)** Volcano plot displaying all genes in the discovery cohort (GSE26440) comparing immunoparalysis (Subclass A) versus non-paralysis (Subclass B+C). X-axis: log2 fold change; Y-axis: -log10(p-value). Horizontal dashed line indicates FDR < 0.05 threshold. **(B)** Validation cohort (GSE26378) volcano plot confirming similar expression patterns.

Genes highlighted: CIITA (red), HLA-DRA (orange), HLA-DPB1 (green). These three genes constitute the minimal pediatric IPS signature.

---

## Figure S2. IPS Score Distribution by Survival Outcome

![IPS Distribution](figures/figS2_ips_distribution.png)

**Figure S2. IPS Scores Stratify Mortality Risk in Pediatric Sepsis**

**(A)** Distribution of IPS scores in the pediatric discovery cohort (GSE26440) stratified by 28-day mortality. High IPS (≥0.31) associated with numerically higher mortality (20.0% vs 15.9%), though not statistically significant (p = 0.592, Fisher's exact test). **(B)** Validation cohort (GSE26378): High IPS group showed 2.7-fold higher mortality (27.3% vs 10.0%).

Dashed line indicates optimal threshold (0.31) determined by Youden's index. Box plots display median (center line), interquartile range (box), and individual data points.

---

## Figure S3. Adult Cohort Model Performance

![Adult ROC](figures/figS3_adult_roc.png)

**Figure S3. Adult IPS Model Performance Across Cohorts**

**(A)** ROC curves comparing adult IPS model performance: Discovery cohort (blue, AUC = 0.830), internal validation (orange, AUC = 0.696), and external validation in GAinS (gray, AUC = 0.550). **(B)** Model coefficient visualization showing contribution of each MHC II gene to the adult IPS score. CIITA (green) contributes positively despite being a marker of immunoparalysis, reflecting the compensatory upregulation pattern.

The reduced performance in external validation likely reflects cross-platform technical variance rather than biological differences.

---

## Figure S4. IFN-γ-STAT1-IRF1 Axis Expression Patterns

![Key Findings](figures/figS4_axis_summary.png)

**Figure S4. Summary of Key Mechanistic Findings**

**(A)** Comparison of IFN-γ, STAT1, and IRF1 expression between pediatric and adult immunoparalysis. Pediatric patients show coordinated downregulation of all pathway components, while adults show maintenance or upregulation of these genes. **(B)** Proposed model explaining the age-dependent differences in CIITA regulation and immune response to sepsis.

This figure summarizes the core finding that immunoparalysis mechanisms differ fundamentally between age groups, with important implications for biomarker development and therapeutic targeting.

---

*Figure legends prepared for: "Age-dependent discordance of CIITA expression in sepsis-induced immunoparalysis: Mechanisms and clinical implications"*
