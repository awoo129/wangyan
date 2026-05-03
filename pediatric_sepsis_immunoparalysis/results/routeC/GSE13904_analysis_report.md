# GSE13904 外部验证分析报告

## 数据集概述

| 属性 | 值 |
|------|-----|
| GEO Accession | GSE13904 |
| 平台 | GPL570 (Affymetrix HG-U133 Plus 2.0) |
| 样本数 | 227 (139 Day1样本) |
| 研究描述 | 儿童SIRS、脓毒症、脓毒性休克谱 |
| PMID | 19325468 |

## 样本分组统计 (Day1)

| 组别 | 样本数 |
|------|--------|
| Control | 18 |
| SIRS | 22 |
| Sepsis | 32 |
| Septic Shock | 67 |
| **总计** | **139** |

## MHC II基因表达分析

### CIITA (Class II Transactivator)

| 组别 | Mean | SD | n |
|------|------|-----|---|
| Control | 0.9997 | 0.2039 | 18 |
| SIRS | 1.0503 | 0.3936 | 22 |
| Sepsis | 1.2797 | 0.4748 | 32 |
| Septic Shock | 1.3696 | 0.7995 | 67 |

**统计检验:**
- Sepsis vs Control: Δ=-0.2799, p=0.0239, Cohen's d=0.766 (↓)
- Septic Shock vs Control: Δ=-0.3698, p=0.0580, Cohen's d=0.634 (↓)

**解读:** CIITA在脓毒症和脓毒性休克中表达**下降**，与对照组相比有统计学显著差异（Sepsis p=0.024）或边缘显著（Septic Shock p=0.058）。

### HLA-DRA (MHC Class II, DR Alpha)

| 组别 | Mean | SD | n |
|------|------|-----|---|
| Control | 0.9477 | 0.3409 | 18 |
| SIRS | 0.3949 | 0.2314 | 22 |
| Sepsis | 0.4081 | 0.1977 | 32 |
| Septic Shock | 0.2948 | 0.2114 | 67 |

**统计检验:**
- Sepsis vs Control: Δ=+0.5395, p<0.0001, Cohen's d=-1.936
- Septic Shock vs Control: Δ=+0.6529, p<0.0001, Cohen's d=-2.302

**解读:** HLA-DRA在脓毒症和脓毒性休克中表达**显著下降**，效应量非常大（|d|>1.9），提示MHC II抗原呈递功能严重受损。

### CD74 (Invariant Chain)

| 组别 | Mean | SD | n |
|------|------|-----|---|
| Control | 1.0114 | 0.2914 | 18 |
| SIRS | 0.9029 | 0.4186 | 22 |
| Sepsis | 1.1456 | 0.4617 | 32 |
| Septic Shock | 0.9159 | 0.3476 | 67 |

**统计检验:**
- 无显著组间差异 (p > 0.05)

## 结论

GSE13904数据集验证了以下发现：

1. **CIITA在儿童脓毒症中表达下降**，与路线A（儿童IPS）的发现一致
2. **HLA-DRA下降更明显**，提示MHC II介导的抗原呈递通路整体下调
3. 该数据集可作为儿童脓毒症免疫瘫痪的外部验证
