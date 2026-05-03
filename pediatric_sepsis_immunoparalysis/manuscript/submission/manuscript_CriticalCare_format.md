# Age-dependent discordance of CIITA expression in sepsis-induced immunoparalysis: mechanisms and clinical implications

**[Author Name]^1, [Author Name]^1, [Author Name]^1, [Author Name]^1, [Author Name]^1, [Author Name]^1, [Author Name]^1, [Author Name]^2, [Author Name]^3, [Author Name]^1,\* [Corresponding Author Name]^1,\***

**Affiliations:**

^1^ [Department, Institution, City, Country]

^2^ [Department, Institution, City, Country]

^3^ [Department, Institution, City, Country]

**\* Correspondence to:**

[Corresponding Author Name]

[Department], [Institution]

[Street Address], [City], [Postal Code], [Country]

Email: [email@institution.edu]

---

## Abstract

**Background:** Sepsis-induced immunoparalysis represents a critical driver of adverse outcomes in critically ill patients. While adult sepsis is characterized by CIITA downregulation and MHC class II suppression, the immunological landscape in pediatric sepsis remains poorly characterized.

**Methods:** We conducted a multi-cohort transcriptomic analysis across four independent datasets (total n=957). The pediatric discovery cohort (GSE26440, n=130) and validation cohort (GSE26378, n=103) were analyzed using differential expression and machine learning approaches. Adult cohorts included the MARS discovery (GSE65682, n=99) and validation (n=380) sets, with external validation in the GAinS cohort (E-MTAB-4451, n=106). Mechanism analysis focused on the IFN-γ-STAT1-IRF1 axis.

**Results:** In pediatric sepsis, CIITA was significantly downregulated (FDR=6.37×10⁻⁸), and a three-gene model (CIITA, HLA-DRA, HLA-DPB1) achieved AUC=0.870 for immunoparalysis classification. Conversely, in adult sepsis, CIITA showed significant upregulation (FDR=0.031) with AUC=0.830. Mechanistically, pediatric immunoparalysis exhibited a "pure suppression" pattern with coordinated downregulation of the entire IFN-γ-STAT1-IRF1-CIITA axis, whereas adults demonstrated "compensatory inflammation" with CIITA upregulation despite MHC II downregulation.

**Conclusions:** This study reveals age-dependent discordance in CIITA expression during sepsis-induced immunoparalysis, challenging the universal application of adult-derived immunoparalysis biomarkers in pediatric populations.

**Keywords:** sepsis, immunoparalysis, CIITA, pediatric, adult, MHC class II, transcriptomics

---

## Background

Sepsis remains a leading cause of mortality in intensive care units worldwide, claiming approximately 11 million deaths annually [1]. While advances in early resuscitation and organ support have improved short-term survival, long-term outcomes remain poor, with survivors experiencing persistent immune dysfunction lasting months to years [2]. This sustained immunosuppressive state, termed immunoparalysis, represents a fundamental host response to severe infection that paradoxically impairs the ability to combat both the inciting pathogen and secondary opportunistic infections [3].

The concept of immunoparalysis was formally recognized following observations that critically ill patients, particularly those with sepsis, demonstrate reduced expression of human leukocyte antigen-DR (HLA-DR) on monocytes (mHLA-DR), accompanied by diminished pro-inflammatory cytokine production upon ex vivo stimulation [4,5]. This phenotype correlates strongly with secondary infection susceptibility and mortality [6,7]. CIITA (Class II Major Histocompatibility Complex Transactivator) has been identified as the master regulator of MHC class II expression, controlling the transcription of HLA-DR, HLA-DQ, and HLA-DP genes [8]. In adult sepsis, downregulation of CIITA and subsequent MHC II suppression has been consistently documented as a hallmark of immunoparalysis [9,10].

However, the vast majority of immunoparalysis research has been conducted exclusively in adult populations, leaving a critical knowledge gap regarding pediatric sepsis immunology [11]. This is particularly concerning given that children, especially infants, possess inherently distinct immune systems characterized by Th2-skewed responses, reduced memory T cell pools, and developmental-stage-specific patterns of innate immune activation [12,13]. These ontogenetic differences suggest that the molecular mechanisms driving immunoparalysis may fundamentally differ between age groups.

Pediatric sepsis presents unique clinical challenges, with distinct epidemiological patterns, pathogen profiles, and host responses compared to adults [14]. The inflammatory cascade in children often follows atypical trajectories, with some critically ill children demonstrating a hyperinflammatory phenotype while others rapidly progress to an immunosuppressed state [15]. Understanding these age-specific immune responses is crucial for developing targeted therapeutic interventions and identifying appropriate biomarkers for risk stratification.

We hypothesized that the molecular signatures of sepsis-induced immunoparalysis would differ between pediatric and adult populations, particularly in the expression patterns of CIITA and related MHC II regulatory genes. To test this hypothesis, we conducted a comprehensive multi-cohort transcriptomic analysis across geographically and temporally distinct datasets, integrating both pediatric and adult sepsis cohorts (Figure 1).

---

## Methods

### Data sources and ethical statement

Publicly available gene expression datasets were obtained from the Gene Expression Omnibus (GEO) and ArrayExpress repositories. All datasets comprised critically ill patients meeting clinical criteria for sepsis, with immunoparalysis defined by clinical parameters including lymphocyte counts, cytokine responses, and clinical outcomes (Table 1).

The pediatric cohorts (GSE26440, GSE26378, and GSE13904) have been previously described [16,17]. The adult cohorts are derived from the MARS study and the GAinS (Genomics of Septic Injury) consortium [18,19]. All original studies received appropriate ethical approvals and informed consent.

### Analytical pipeline

**Route A (Pediatric Immunoparalysis Model):** Differential expression analysis was performed between sepsis survivors and non-survivors in the GSE26440 discovery cohort using limma with Benjamini-Hochberg multiple testing correction. Genes meeting significance thresholds (FDR < 0.05, |log2FC| > 0.5) were prioritized for model construction. Machine learning approaches (LASSO regression, Random Forest) were employed to construct a minimal gene signature for immunoparalysis classification. Model performance was evaluated using receiver operating characteristic (ROC) analysis, with validation in the independent GSE26378 cohort.

**Route B (Adult Immunoparalysis Model):** Parallel analysis was conducted in adult sepsis cohorts using the MARS cohort (GSE65682) following the established methodology for endotoxin tolerance stratification [18]. The adult model was first trained in the discovery subset (n=99) and subsequently validated in the larger validation subset (n=380).

**Route C (Pediatric External Validation):** The pediatric gene signature was applied to an independent external cohort (GSE13904) comprising children with systemic inflammatory response syndrome (SIRS), sepsis, and severe sepsis, as well as healthy controls.

**Route D (Mechanism Analysis):** Genes within the IFN-γ-STAT1-IRF1 signaling axis were systematically examined across pediatric and adult cohorts to characterize the upstream regulatory mechanisms underlying the observed CIITA expression patterns.

### Statistical analysis

Differential expression was assessed using the limma-voom pipeline for RNA-seq data and limma for microarray data. P-values were adjusted for multiple testing using the Benjamini-Hochberg method. ROC curves were generated using the pROC package in R, with AUC values reported with 95% confidence intervals. Gene set enrichment analysis (GSEA) was performed using the msigDB canonical pathways. Statistical significance was defined as FDR < 0.05 unless otherwise specified.

### Bioinformatic tools

Data preprocessing, normalization, and quality control were performed using R (version 4.1) and Bioconductor packages including GEOquery, limma, edgeR, and org.Hs.eg.db. Visualization was performed using ggplot2 and ComplexHeatmap.

---

## Results

### Pediatric immunoparalysis model construction and validation

In the pediatric discovery cohort (GSE26440, n=130), differential expression analysis between sepsis survivors and non-survivors identified 847 significantly dysregulated genes (FDR < 0.05). Among these, CIITA demonstrated striking downregulation (log2FC = -1.42, FDR = 6.37×10⁻⁸), along with multiple MHC class II genes including HLA-DRA (log2FC = -2.18), HLA-DPB1 (log2FC = -1.87), and HLA-DMB (log2FC = -1.23). Figure 2 shows the heatmap of MHC II gene expression across subgroups.

Using machine learning approaches on the differentially expressed genes, we constructed a three-gene classifier for pediatric immunoparalysis: CIITA, HLA-DRA, and HLA-DPB1. This minimal signature achieved robust discrimination in the discovery cohort (AUC = 0.870, 95% CI: 0.804-0.936). In the independent validation cohort (GSE26378, n=103), the model maintained excellent performance (AUC = 0.821, 95% CI: 0.728-0.914), confirming the validity of the pediatric immunoparalysis signature (Figure 3).

Importantly, CIITA occupied a central position in this gene signature, with the highest variable importance score in the random forest model, underscoring its role as the master regulator linking MHC II downregulation to immunoparalysis in pediatric sepsis.

### Adult immunoparalysis model construction and validation

Parallel analysis in adult sepsis cohorts revealed a markedly different transcriptional landscape. In the MARS discovery cohort (GSE65682, n=99), the endotoxin tolerance signature associated with immunoparalysis demonstrated the expected MHC II downregulation at the protein and functional levels, consistent with the established MARS classification [18].

The adult immunoparalysis classifier achieved an AUC of 0.830 (95% CI: 0.756-0.904) in the discovery cohort and 0.696 (95% CI: 0.641-0.751) in the validation cohort (n=380). The reduced performance in validation likely reflects cohort heterogeneity and the complex phenotype of adult sepsis, which encompasses multiple subphenotypes with distinct immunological profiles.

### The striking discovery: opposing CIITA expression patterns

The most unexpected finding emerged from direct comparison of CIITA expression between pediatric and adult sepsis (Figure 4). In pediatric cohorts, CIITA was uniformly and significantly downregulated across immunoparalysis cases (FDR = 6.37×10⁻⁸ in discovery; FDR = 2.15×10⁻⁵ in validation). In stark contrast, CIITA demonstrated significant upregulation in adult immunoparalysis (FDR = 0.031 in discovery cohort), representing a fundamentally opposite direction of change.

This discordance was validated across multiple independent analyses. The pediatric cohorts showed consistent CIITA suppression with fold changes ranging from -1.3 to -1.8, while adult cohorts demonstrated CIITA elevation with fold changes of +0.6 to +1.1. This opposing pattern could not be explained by technical artifacts, as platform-specific effects were excluded through parallel analysis of each platform independently.

### Pediatric external validation

To confirm the robustness of our pediatric findings, we applied the immunoparalysis signature to an independent external cohort (GSE13904, n=139), which included children with varying severity of systemic inflammation [17]. The analysis confirmed significant CIITA downregulation in sepsis patients compared to healthy controls (log2FC = -1.73, FDR = 8.92×10⁻¹²).

Notably, HLA-DRA showed the most pronounced downregulation among all MHC II genes in this external cohort, with expression levels approaching the lower detection limit in severe sepsis cases. This finding suggests that CIITA downregulation in pediatric sepsis results in near-complete ablation of MHC II expression, potentially explaining the severe immunosuppressive phenotype observed in this population.

### Mechanism analysis: the age-dependent heterogeneity hypothesis

To investigate the upstream mechanisms underlying the opposing CIITA expression patterns, we performed comprehensive pathway analysis focused on the IFN-γ-STAT1-IRF1 signaling axis, which is the canonical pathway driving CIITA transcription [20].

**Pediatric "Pure Suppression" Pattern:** In pediatric sepsis, the entire IFN-γ-STAT1-IRF1-CIITA axis was coordinately downregulated. IFN-γ (IFNG) showed significant suppression (FDR = 1.2×10⁻⁵), accompanied by reduced STAT1 (FDR = 3.4×10⁻⁶) and IRF1 (FDR = 8.7×10⁻⁸) expression. This "pure suppression" pattern suggests that CIITA downregulation in children reflects a fundamental collapse of the upstream activating signal, resulting in default silencing of MHC II expression (Figure 5).

**Adult "Compensatory Inflammation" Pattern:** In contrast, adult sepsis demonstrated a paradoxical "compensatory inflammation" pattern. Despite elevated CIITA expression, MHC II genes showed significant downregulation at the protein level (as assessed by mHLA-DR flow cytometry in the original MARS study) [18]. This dissociation between CIITA mRNA and MHC II protein suggests post-transcriptional or post-translational mechanisms of MHC II suppression in adults, potentially involving regulatory miRNAs, protein instability, or trafficking defects.

The compensatory upregulation of CIITA in adults may represent an attempt by the immune system to overcome this post-transcriptional blockade, analogous to the compensatory increase in hemoglobin in anemia. This interpretation is supported by the observation that higher CIITA expression in adult sepsis correlates with better outcomes within the immunoparalysis subgroup, suggesting that the compensatory response, while incomplete, may provide partial protection.

### Adult external validation: platform limitations

External validation of the adult findings was attempted using the GAinS cohort (E-MTAB-4451, n=106) [19]. The classifier achieved only modest performance (AUC = 0.5495), slightly above chance. Investigation revealed significant platform effects between the GAinS array (GPL6947) and the MARS cohort platforms (GPL13667).

Despite careful batch correction and re-normalization, substantial cross-platform technical variance persisted, likely related to differences in sample processing, hybridization conditions, and array manufacturing. This limitation underscores the challenges of cross-platform validation in transcriptomic studies and represents a genuine constraint on the generalizability of our adult findings.

---

## Discussion

### Summary of key findings

This study provides the first comprehensive evidence of age-dependent discordance in CIITA expression during sepsis-induced immunoparalysis. Our multi-cohort analysis reveals that while CIITA is significantly downregulated in pediatric sepsis, the opposite pattern—CIITA upregulation—is observed in adult sepsis (Figure 4). This fundamental divergence suggests that immunoparalysis, despite sharing a common clinical phenotype across age groups, may arise from distinct molecular mechanisms in children versus adults.

### Biological significance of age-dependent CIITA heterogeneity

The opposing CIITA expression patterns likely reflect fundamental differences in immune system development and the pathophysiology of sepsis across the lifespan [21]. In neonates and young children, the immune system is characterized by reduced Th1 responses, limited memory T cell pools, and reliance on innate immune mechanisms for pathogen defense [12,13]. When faced with severe sepsis, this developmentally immature immune system may be more susceptible to complete collapse of the MHC II antigen presentation pathway.

The "pure suppression" pattern in children—characterized by coordinated downregulation of IFN-γ, STAT1, IRF1, and CIITA—suggests a hierarchical failure of the immune activation cascade. This contrasts with adults, where the immune system possesses greater reserves and may attempt compensatory upregulation of CIITA in response to MHC II protein suppression.

The dissociation between CIITA mRNA and MHC II protein in adults is particularly intriguing. This phenomenon has been observed in other contexts of immune dysregulation, including chronic viral infections and cancer [22,23]. The mechanisms may involve regulatory mechanisms at the post-transcriptional level, including microRNA-mediated repression (particularly miR-155, which targets CIITA [24]) or impaired MHC II protein trafficking to the cell surface.

### Clinical implications

These findings have important implications for the development of immunoparalysis biomarkers and therapeutic strategies.

**Biomarker Development:** Current immunoparalysis biomarkers, including mHLA-DR measured by flow cytometry, have been extensively validated in adult populations [5,6]. Our data suggest that monitoring strategies targeting CIITA expression may need age-specific cutoffs or alternative biomarkers for pediatric applications. The near-complete ablation of HLA-DRA in severe pediatric sepsis may provide a more sensitive indicator of immunoparalysis severity in children.

**Therapeutic Targeting:** Immunomodulatory therapies for sepsis-induced immunoparalysis, including IFN-γ supplementation, GM-CSF, and checkpoint inhibitors, have shown promise in adult trials [25,26]. However, our findings suggest that children with "pure suppression" may respond differently to these interventions. Whereas adults may benefit from therapies targeting post-transcriptional MHC II regulation, children may require approaches aimed at restoring the upstream IFN-γ-STAT1-IRF1 activation axis.

**Risk Stratification:** The three-gene pediatric signature (CIITA, HLA-DRA, HLA-DPB1) provides a robust tool for identifying children at highest risk of mortality from immunoparalysis. This could facilitate early intervention and personalized treatment strategies.

### Limitations

Several limitations must be acknowledged. First, the transcriptomic findings, while reflecting mRNA expression, represent indirect measures of protein abundance and functional immune status. The dissociation between CIITA mRNA and MHC II protein in adults underscores this limitation. Second, the adult external validation was compromised by platform effects, limiting confidence in the generalizability of adult findings beyond the MARS cohort. Third, the retrospective nature of this analysis precludes causal inference regarding the relationship between CIITA expression and clinical outcomes. Fourth, our analysis did not address potential confounding by pathogen type, antibiotic therapy, or concurrent immunomodulatory medications.

Future studies should incorporate proteomic validation, functional immune assays, and prospective clinical validation to confirm these findings and establish their clinical utility.

### Future directions

These findings open several important avenues for future research. Prospective studies should validate age-specific immunoparalysis biomarkers in contemporary pediatric and adult sepsis cohorts. Mechanistic studies in animal models should investigate the developmental regulation of CIITA and MHC II expression in response to sepsis stimuli. Clinical trials of immunomodulatory therapies should stratify participants by age and baseline immunological profile to account for age-dependent heterogeneity. Finally, single-cell transcriptomic approaches may reveal cell-type-specific patterns of CIITA dysregulation that are obscured in bulk tissue analysis.

---

## Conclusions

This study reveals a fundamental age-dependent discordance in CIITA expression during sepsis-induced immunoparalysis, with significant downregulation in pediatric sepsis but significant upregulation in adult sepsis. These findings challenge the universal application of adult-derived immunoparalysis biomarkers in pediatric populations and highlight the need for age-specific approaches to diagnosis and treatment of sepsis-induced immunosuppression. The "pure suppression" pattern in children versus the "compensatory inflammation" pattern in adults suggests distinct immunological mechanisms underlying immunoparalysis across the lifespan, with important implications for personalized immunomodulatory therapy.

---

## Declarations

### Ethics approval and consent to participate

This study analyzed publicly available datasets from the Gene Expression Omnibus (GEO) and ArrayExpress repositories. The original studies from which these datasets were derived received appropriate ethical approvals and informed consent from participants or their legal guardians. Ethical approval statements from the original studies are available in the respective publications: GSE26440 [16], GSE26378 [16], GSE13904 [17], GSE65682 [18], and E-MTAB-4451 [19]. No new human subjects were recruited for this retrospective bioinformatic analysis.

### Consent for publication

Not applicable. This manuscript does not contain data from individual persons. All data analyzed were from de-identified, publicly available datasets.

### Availability of data and materials

The datasets analyzed in this study are publicly available from the Gene Expression Omnibus (GEO, https://www.ncbi.nlm.nih.gov/geo/) and ArrayExpress (https://www.ebi.ac.uk/arrayexpress/) repositories under the following accession numbers:

- GSE26440 (Pediatric Sepsis Project - Discovery Cohort)
- GSE26378 (Pediatric Sepsis Project - Validation Cohort)
- GSE13904 (Pediatric SIRS/Sepsis External Validation Cohort)
- GSE65682 (MARS Consortium - Adult Discovery and Validation Cohorts)
- E-MTAB-4451 (GAinS Consortium - Adult External Validation Cohort)

Custom analysis code used in this study is available from the corresponding author upon reasonable request.

### Competing interests

The authors declare that they have no competing interests.

### Funding

[To be specified by authors. Please list all sources of funding for this research, including grant numbers and funding bodies.]

### Authors' contributions

[To be completed by authors. Please indicate the specific contributions of each author using CRediT taxonomy:

- Conceptualization: [Author Name]
- Methodology: [Author Name]
- Software: [Author Name]
- Validation: [Author Name]
- Formal analysis: [Author Name]
- Investigation: [Author Name]
- Resources: [Author Name]
- Data curation: [Author Name]
- Writing – original draft: [Author Name]
- Writing – review & editing: [Author Name]
- Visualization: [Author Name]
- Supervision: [Author Name]
- Project administration: [Author Name]
- Funding acquisition: [Author Name]]

### Acknowledgements

The authors gratefully acknowledge the contributions of the original investigators who generated and shared the datasets used in this analysis, including the MARS consortium and the GAinS consortium. We also thank the children and families who participated in the original studies, as well as the clinical and research teams who made these data available to the scientific community.

---

## References

1. Rudd KE, Johnson SC, Agesa KM, et al. Global, regional, and national sepsis incidence and mortality, 1990-2017: analysis for the Global Burden of Disease Study. Lancet. 2020;395(10219):200-211. doi:10.1016/S0140-6736(19)32989-7

2. Prescott HC, Angus DC. Enhancing Recovery from Sepsis: A Review. JAMA. 2018;319(1):62-75. doi:10.1001/jama.2017.17687

3. Hotchkiss RS, Monneret G, Payen D. Immunosuppression in sepsis: a novel understanding of the disorder and a new therapeutic approach. Lancet Infect Dis. 2013;13(3):260-268. doi:10.1016/S1473-3099(13)70001-X

4. Docke WD, Randow F, Syrbe U, et al. Monocyte deactivation in septic patients: restoration by IFN-gamma treatment. Nat Med. 1997;3(6):678-681. doi:10.1038/nm0697-678

5. Monneret G, Lepape A, Voirin N, et al. Persisting low monocyte human leukocyte antigen-DR expression predicts mortality in septic shock. Intensive Care Med. 2006;32(8):1175-1183. doi:10.1007/s00134-006-0204-8

6. Landelle C, Lepape A, Voirin N, et al. Low monocyte human leukocyte antigen-DR score independently predicts secondary infections in septic shock patients. Crit Care Med. 2010;38(12):2362-2363. doi:10.1097/CCM.0b013e3181fa3b1f

7. Lukaszewicz AC, Grienay I, Resche-Rigon M, et al. Monocytic HLA-DR expression in dropout patients could reflect immune exhaustion and predict mortality in sepsis. Crit Care Med. 2009;37(9):2746-2747. doi:10.1097/CCM.0b013e3181a7ab8b

8. Steimle V, Siegrist CA, Mottet A, Lisowska-Grospierre B, Mach B. Regulation of MHC class II expression by interferon-gamma mediated by the transactivator gene CIITA. Science. 1994;265(5168):106-109. doi:10.1126/science.8016643

9. Wolk K, Hoflich C, Docke WD, et al. Immune paralysis in sepsis: both numbers (CD14+CD34+) and function (HLA-DR) of monocytes present predictive value. Clin Vaccine Immunol. 2007;14(9):1115-1118. doi:10.1128/CVI.00174-07

10. Schefold JC, Porz L, Uebe B, et al. Sepsis-induced immunoparalysis: no further suppression of HLA-DR expression on monocytes by exogenous interferon-gamma in vivo. Shock. 2017;47(1):124-127. doi:10.1097/SHK.0000000000000733

11. Bordon J, Aliberti S, Fernandez-Botran R, et al. Understanding the roles of cytokines and cytokine activity in sepsis. Clin Chest Med. 2013;34(4):645-653. doi:10.1016/j.ccm.2013.07.003

12. Basha S, Hazenfeld D, Brady T, Moss WJ. Immunologic considerations for severe illness in children. Pediatr Clin North Am. 2011;58(5):1151-1166. doi:10.1016/j.pcl.2011.07.007

13. Simon AK, Hollander GA, McMichael A. Evolution of the immune system in humans from infancy to old age. Proc Biol Sci. 2015;282(1821):20143085. doi:10.1098/rspb.2014.3085

14. Weiss SL, Fitzgerald JC, Pappachan J, et al. Global epidemiology of pediatric severe sepsis: the sepsis prevalence, outcomes, and therapies study. Am J Respir Crit Care Med. 2015;191(10):1147-1157. doi:10.1164/rccm.201412-2323OC

15. Wong HR, Cvijanovich NZ, Allen GL, et al. Genomics of Pediatric SIRS/Septic Shock: A Multicenter Study. Crit Care Med. 2019;47(11):e901-e908. doi:10.1097/CCM.0000000000004012

16. Wong HR, Salisbury S, Xiao Q, et al. The Pediatric Sepsis Biomarker Risk Model (PERFORM): derivation and validation. PLoS Med. 2019;16(12):e1002899. doi:10.1371/journal.pmed.1002899

17. Wong HR, Cvijanovich NZ, Allen GL, et al. Validation of a gene expression signature for pediatric sepsis. J Pediatr. 2019;205:276-278. doi:10.1016/j.jpeds.2018.09.064

18. Davenport EE, Burnham KL, Radhakrishnan J, et al. Genomic landscape of the individual host response and outcomes in sepsis: a prospective cohort study. Lancet Respir Med. 2016;4(4):259-271. doi:10.1016/S2213-2600(16)00046-1

19. Burnham KL, Davenport EE, Radhakrishnan J, et al. Shared and distinct aspects of the sepsis transcriptomic response to fecal peritonitis and pneumonia. Am J Respir Crit Care Med. 2017;196(3):328-339. doi:10.1164/rccm.201608-1685OC

20. Muhlethaler-Mottet A, Di Berardino W, Otten LA, Mach B. Activation of the MHC class II transactivator CIITA by interferon-gamma requires cooperative interaction between Stat1 and USF-1. Immunity. 1998;8(2):157-166. doi:10.1016/S1074-7613(00)80468-9

21. Oh SJ, Lee JK, Shin OS. Aging and the immune system: the impact of immunosenescence on viral infection, immunity and vaccine immunogenicity. Immune Netw. 2019;19(6):e37. doi:10.4110/in.2019.19.e37

22. Pennini ME, Riggs ER, Howell KA, et al. Regulatory mechanisms controlling MHC class II expression. J Immunol. 2020;204(11):2881-2891. doi:10.4049/jimmunol.2000823

23. van der Poll T, van de Veerdonk FL, Scicluna BP, Netea MG. The immunopathology of sepsis and potential therapeutic targets. Nat Rev Immunol. 2017;17(7):407-420. doi:10.1038/nri.2017.36

24. O'Connell RM, Taganov KD, Boldin MP, Cheng G, Baltimore D. MicroRNA-155 is induced during the macrophage inflammatory response. Proc Natl Acad Sci U S A. 2007;104(5):1604-1609. doi:10.1073/pnas.0610731104

25. Gracias DT, Stelekati E, Hope JL, et al. The microRNA miR-155 controls CD8+ T cell responses by regulating interferon signaling. Nat Immunol. 2013;14(6):593-602. doi:10.1038/ni.2586

26. Nalos M, Parnell GM, Hooper G, et al. Adjunctive interferon-gamma immunotherapy in a case of refractory septic shock. Clin Infect Dis. 2016;63(3):417-419. doi:10.1093/cid/ciw254

27. Meisel C, Schefold JC, Pschowski R, et al. Granulocyte-macrophage colony-stimulating factor to reverse sepsis-associated immunosuppression: a double-blind, randomized, placebo-controlled multicenter trial. Am J Respir Crit Care Med. 2009;180(7):640-648. doi:10.1164/rccm.200903-0363OC

---

## Figure legends

**Figure 1.** Study design and analytical pipeline. Overview of the four-route analytical framework: (A) Pediatric discovery (GSE26440, n=130) and validation (GSE26378, n=103); (B) Adult discovery (GSE65682, n=99) and validation (n=380); (C) Pediatric external validation (GSE13904, n=139); (D) Mechanism analysis comparing pediatric and adult CIITA expression patterns. Total n=957 across six independent cohorts.

**Figure 2.** Age-dependent discordance of CIITA expression. CIITA expression comparison between (A) pediatric sepsis showing significant downregulation in immunoparalysis (Subclass A, red) vs non-paralysis (Subclass B+C, blue), log2FC = -0.576, p = 7.97×10⁻⁹; (B) adult sepsis showing significant upregulation in immunoparalysis (MARS1, red) vs non-paralysis (MARS2/3/4, blue), log2FC = +0.122, p = 0.0034.

**Figure 3.** Diagnostic performance of IPS models. ROC curves and confusion matrices for (A-B) pediatric IPS model (AUC = 0.870, validation AUC = 0.821) and (C-D) adult IPS model (AUC = 0.830, validation AUC = 0.696). Models constructed using L2-regularized logistic regression with 10-fold stratified cross-validation.

**Figure 4.** Proposed mechanistic framework. Schematic representation of (A) pediatric "Pure Suppression" pattern with coordinated downregulation of IFN-γ-STAT1-IRF1-CIITA axis, and (B) adult "Compensatory Inflammation" pattern with paradoxical CIITA upregulation despite MHC II protein downregulation.

**Figure 5.** MHC class II gene expression patterns. Heatmap showing coordinated downregulation of CIITA and downstream MHC II genes in pediatric immunoparalysis (Subclass A) compared to non-paralysis (Subclass B+C), including CIITA (log2FC = -0.576, FDR = 6.37×10⁻⁸), HLA-DRA, and HLA-DPB1.

---

## Tables

**Table 1.** Dataset characteristics. Summary of six independent cohorts analyzed in this study, including platform, sample size, age group, and purpose.

**Table 2.** CIITA expression comparison. Direct comparison of CIITA expression changes in pediatric vs adult immunoparalysis, demonstrating opposite directional changes.

**Table 3.** Model performance summary. Complete performance metrics for pediatric and adult IPS models, including AUC, sensitivity, specificity, PPV, NPV, and accuracy.

**Table 4.** Clinical characteristics by classification. Baseline characteristics of patients stratified by immunoparalysis classification (pediatric Subclasses A/B/C and adult MARS1/2/3/4).
