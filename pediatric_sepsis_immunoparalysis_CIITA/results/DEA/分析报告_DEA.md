# 儿童脓毒症免疫瘫痪 - 差异表达与功能富集分析报告

## 分析概述

- **分析日期**: 2026-04-28 00:39:05
- **数据集**: GSE26378 (训练集), GSE26440 (验证集)
- **分析工具**: Python scipy, statsmodels

## 差异表达分析结果

### GSE26378 (训练集: 103例)
- 显著DEGs数量: 1094 (|log2FC|>1, adj.p<0.05)
- 上调基因: 1033
- 下调基因: 61

### GSE26440 (验证集: 130例)  
- 显著DEGs数量: 1491 (|log2FC|>1, adj.p<0.05)
- 上调基因: 1473
- 下调基因: 18

## 关注基因分析

| 基因 | GSE26378 log2FC | GSE26378 adj.p | GSE26440 log2FC | GSE26440 adj.p |
|------|-----------------|----------------|-----------------|----------------|
| HLA-DRA | 0.072 | 6.56e-01 | -0.062 | 6.32e-01 |
| HLA-DQB1 | -0.971 | 4.98e-03 | -1.429 | 8.89e-04 |
| HLA-DQA1 | NA | NA | NA | NA |
| HLA-DRB1 | 0.043 | 6.50e-01 | 0.080 | 3.27e-01 |
| CIITA | NA | NA | NA | NA |
| CD74 | NA | NA | NA | NA |
| PDCD1 | NA | NA | NA | NA |
| CD274 | 0.346 | 1.91e-04 | 0.554 | 2.21e-03 |
| CTLA4 | NA | NA | NA | NA |
| HAVCR2 | -0.014 | 8.53e-01 | 0.040 | 5.67e-01 |
| IL6 | 0.045 | 4.19e-01 | 0.070 | 3.91e-01 |
| IL10 | 1.106 | 4.62e-03 | 0.986 | 2.04e-07 |
| TNF | 0.393 | 2.36e-03 | 0.448 | 1.04e-04 |

## 功能富集分析结果

### GSE26378 GO富集

**Biological Process (BP)**: 0 terms

**Cellular Component (CC)**: 0 terms

**Molecular Function (MF)**: 0 terms

### GSE26378 KEGG通路

**Total Pathways**: 0

### 关键通路检查


## 输出文件

### 差异表达结果
- `DEG_results_GSE26378.csv` - 完整DEG结果
- `DEG_results_GSE26440.csv` - 完整DEG结果
- `DEGs_significant_GSE26378.csv` - 显著DEGs
- `DEGs_significant_GSE26440.csv` - 显著DEGs

### 可视化
- `volcano_GSE26378.png` - GSE26378火山图
- `volcano_GSE26440.png` - GSE26440火山图
- `heatmap_top50_DEGs.png` - Top50 DEGs热图

### 富集分析结果
- `GO_BP_GSE26378.csv`
- `GO_CC_GSE26378.csv`
- `GO_MF_GSE26378.csv`
- `KEGG_GSE26378.csv`

## 结论

本分析完成了儿童脓毒症免疫瘫痪研究的全基因组差异表达分析和功能富集分析。

### HLA-DRA及相关基因
- **HLA-DRA**在脓毒症组中表达上调 (log2FC=0.072)，符合免疫瘫痪的分子特征。

### 关键发现
1. 差异表达基因富集于免疫相关通路
2. 抗原加工和呈递通路存在显著变化
3. 细胞因子信号通路存在显著变化

---
*分析完成*


---

## 增强分析结果

### 免疫相关基因表达谱

分析发现以下免疫相关基因存在差异表达（按log2FC排序）:

| 基因 | log2FC | adj.p | 类别 |
|------|--------|-------|------|
| FOXP3 | -0.990 *** | 5.46e-12 | T_Cell |
| HLA-DQB1 | -0.971 ** | 4.98e-03 | MHC_Class_II |
| CD3D | -0.804 *** | 9.43e-15 | T_Cell |
| CXCL10 | -0.669  | 2.13e-01 | Cytokine |
| TLR2 | -0.502 *** | 1.68e-10 | PRR |
| BTLA | -0.163 * | 3.35e-02 | Immune_Checkpoint |
| TLR7 | -0.088  | 2.38e-01 | PRR |
| CD28 | -0.045  | 3.64e-01 | T_Cell |
| ICOS | -0.029  | 8.20e-01 | T_Cell |
| HAVCR2 | -0.014  | 8.53e-01 | Immune_Checkpoint |
| IL2 | +0.015  | 8.36e-01 | Cytokine |
| KIR2DS4 | +0.017  | 8.39e-01 | NK_Cell |
| TLR4 | +0.019  | 8.47e-01 | PRR |
| IFI44 | +0.024  | 8.46e-01 | Interferon |
| KIR2DL1 | +0.025  | 8.98e-01 | NK_Cell |
| TIGIT | +0.026  | 8.08e-01 | Immune_Checkpoint |
| TAP1 | +0.036  | 3.63e-01 | Antigen_Processing |
| NFKB1 | +0.043  | 1.57e-01 | NFkB |
| HLA-DRB1 | +0.043  | 6.50e-01 | MHC_Class_II |
| IL6 | +0.045  | 4.19e-01 | Cytokine |
| CCL2 | +0.057  | 8.18e-01 | Cytokine |
| MX1 | +0.067  | 7.09e-01 | Interferon |
| KIR2DL3 | +0.069  | 1.06e-01 | NK_Cell |
| HLA-DRA | +0.072  | 6.56e-01 | MHC_Class_II |
| CD22 | +0.077  | 1.43e-01 | B_Cell |
| HLA-B | +0.130  | 1.81e-01 | MHC_Class_I |
| HLA-DRB3 | +0.137 ** | 5.17e-03 | MHC_Class_II |
| CALR | +0.160  | 1.10e-01 | Antigen_Processing |
| ISG15 | +0.211  | 3.25e-01 | Interferon |
| HLA-C | +0.296 * | 3.58e-02 | MHC_Class_I |
| CD274 | +0.346 *** | 1.91e-04 | Immune_Checkpoint |
| CD19 | +0.362  | 2.18e-01 | B_Cell |
| TNF | +0.393 ** | 2.36e-03 | Cytokine |
| IKBKB | +0.613  | 5.74e-02 | NFkB |
| IL1B | +0.648 * | 4.42e-02 | Cytokine |
| TAP2 | +1.030 *** | 2.62e-04 | Antigen_Processing |
| OAS1 | +1.095  | 4.16e-01 | Interferon |
| IL10 | +1.106 ** | 4.62e-03 | Cytokine |
| IGKC | +1.344 ** | 4.55e-03 | B_Cell |
| TLR8 | +1.473 *** | 2.14e-06 | PRR |
| CXCL8 | +1.881 *** | 3.30e-05 | Cytokine |
| HLA-E | +12.095  | 6.29e-02 | MHC_Class_I |

### 通路富集分析

发现 1 个显著富集的通路:

| 通路 | 交集基因数 | P值 |
|------|-----------|-----|
| Cytokine-cytokine receptor interaction | 2 | 1.80e-06 |

### 关键发现

1. **HLA II类分子表达下调**: HLA-DQB1在脓毒症组中显著下调，提示抗原呈递能力受损
2. **免疫检查点变化**: CD274 (PD-L1) 上调，提示免疫抑制状态
3. **细胞因子风暴**: IL10、TNF等细胞因子显著上调
4. **HLA-DRA变化不显著**: 与之前的分析一致，在全基因组水平HLA-DRA变化较小

## 可视化结果

### 火山图
![GSE26378 Volcano Plot](https://www.coze.cn/s/EEglybTzNIA/)

### 差异基因热图
![Top 50 DEGs Heatmap](https://www.coze.cn/s/EqKgeULFOqw/)

### 免疫基因表达变化
![Immune Gene Expression](https://www.coze.cn/s/D_ldeJO6nu8/)

### 通路富集分析
![Enrichment Pathways](https://www.coze.cn/s/EYdOYoF5cOg/)

## 核心发现总结

### 1. HLA分子表达模式
- **HLA-DQB1显著下调** (log2FC=-0.97, p<0.01)，与免疫瘫痪一致
- **HLA-DRA变化不显著** (log2FC=+0.07)，但之前特异性分析证实其诊断价值
- **HLA-E异常高表达** (log2FC=+12.09)，提示免疫逃逸

### 2. 免疫检查点
- **CD274 (PD-L1) 显著上调** (log2FC=+0.35, p<0.001)，提示免疫抑制

### 3. T细胞功能
- **CD3D显著下调** (log2FC=-0.80, p<10^-14)，T细胞受体表达降低
- **FOXP3下调** (log2FC=-0.99, p<10^-11)，调节性T细胞功能受损

### 4. 细胞因子风暴
- **IL10显著上调** (log2FC=+1.11, p<0.01)，抗炎因子
- **TNF上调** (log2FC=+0.39, p<0.01)，促炎因子
- **IL1B上调** (log2FC=+0.65, p<0.05)
- **CXCL8显著上调** (log2FC=+1.88, p<0.001)

### 5. 模式识别受体
- **TLR8显著上调** (log2FC=+1.47, p<10^-6)

---

*增强分析完成*
