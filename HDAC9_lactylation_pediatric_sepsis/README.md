# Multi-omics identification of HDAC9 downregulation and a lactate-associated metabolic signature in pediatric septic shock

**Journal:** BMC Bioinformatics

**Authors:** Yan Wang, Lan Huang

## Contents

| File | Description |
|------|-------------|
| `R_scripts/external_validation.R` | External validation analysis using GSE26440 and GSE65682 cohorts |
| `R_scripts/revision_R_scripts.R` | Revision scripts for resubmission (response to reviewers) |
| `R_scripts/validate_gse65682.R` | GSE65682 adult cohort validation analysis |
| `Table1.csv` | Baseline characteristics of study cohorts |
| `Table2.csv` | Univariate and multivariate Cox regression results |
| `Supplementary_Table_S1.csv` | Dataset characteristics |
| `Supplementary_Table_S2.csv` | Full GSEA results |
| `Supplementary_Table_S3.csv` | Metabolomics differential analysis |
| `Supplementary_Table_S4.csv` | Consensus clustering results |
| `Supplementary_Table_S5.csv` | Single-cell RNA sequencing analysis results |

## Reproducibility

1. Clone this repository
2. Run R scripts in order:
   - `R_scripts/external_validation.R`
   - `R_scripts/validate_gse65682.R`
   - `R_scripts/revision_R_scripts.R`
3. Processed expression matrices and sessionInfo are available on Zenodo (DOI: TBD)

## Data Availability

Raw transcriptomic data are available from GEO:
- GSE66099, GSE26440, GSE25504, GSE65682
- Single-cell: GSE167363
- Metabolomics: ST003136, ST003439 (MetaboLights)

## Contact

Yan Wang - 245198067@qq.com
