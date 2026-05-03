# 儿童脓毒症免疫瘫痪研究 - GEO数据分析项目

## 项目概述

本研究旨在通过GEO公共数据库筛选儿童脓毒症免疫瘫痪的分子标志物，采用纯生信分析方法构建免疫瘫痪评分(IPS)。

## 目录结构

```
儿童脓毒症免疫瘫痪研究/
├── data/
│   ├── raw/                    # 原始GEO数据
│   │   ├── GSE26378_series_matrix.txt.gz
│   │   ├── GSE26440_series_matrix.txt.gz
│   │   └── GSE13904_series_matrix.txt.gz
│   ├── normalized/             # 标准化后数据
│   │   ├── GSE26378_expression.csv
│   │   ├── GSE26440_expression.csv
│   │   ├── GSE13904_expression.csv
│   │   └── all_data.pkl
│   └── metadata/               # 元数据
│       ├── GSE26378_metadata.csv
│       ├── GSE26440_metadata.csv
│       └── GSE13904_metadata.csv
├── results/
│   ├── QC/                     # 质控报告
│   │   ├── GSE26378/
│   │   ├── GSE26440/
│   │   └── GSE13904/
│   └── preprocessing/          # 预处理结果
│       ├── dataset_summary.csv
│       ├── mhc_gene_detection.csv
│       └── preprocessing_report.md
└── scripts/                    # 分析脚本
    ├── download_and_preprocess.py
    ├── differential_expression.R
    ├── build_ips_score.R
    └── survival_analysis.R
```

## 研究策略（v2更新）

采用**路线A + 路线B并行**策略，详见 `docs/研究路线图_v2.md`

### 路线A：儿童脓毒症IPS评分研究
- **Derivation**: GSE26440 (130样本，子类A/B/C分型)
- **Validation**: GSE26378 (103样本，子类A/B/C分型)
- **目标**: 子类A（免疫瘫痪型）vs B/C的IPS评分构建

### 路线B：成人脓毒症IPS扩展验证
- **Discovery**: GSE65682 Discovery队列 (263样本，Mars1/2/3/4分型)
- **Validation**: GSE65682 Validation队列 (216样本)
- **目标**: Mars1（免疫瘫痪型，34%死亡率）vs 其他类型的验证

### 核心数据集

| 数据集 | 样本量 | 平台 | 分型标签 | 用途 |
|--------|--------|------|----------|------|
| GSE26440 | 130 | GPL570 | 子类A/B/C | 路线A训练 |
| GSE26378 | 103 | GPL570 | 子类A/B/C | 路线A验证 |
| GSE65682 | 479(有标签) | GPL13667 | Mars1/2/3/4 | 路线B |
| GSE13904 | 227 | GPL570 | - | 备用 |
| GSE66099 | 276 | GPL570 | - | 备用 |

## 关键基因

### MHC II类基因（免疫瘫痪核心）
- **HLA-DRA** (主要目标) - MHC II类DRα链
- HLA-DQB1 - DQβ链
- HLA-DQA1 - DQα链
- HLA-DRB1 - DRβ链
- HLA-DPA1, HLA-DPB1 - DP链
- CIITA - MHC II类转录激活因子
- CD74 - MHC II类不变链

### 免疫相关基因
- 炎症因子: IL6, TNF, IL1B
- 免疫检查点: PDL1, CTLA4
- 免疫细胞标记: CD4, CD8, CD14

## 分析流程

### Phase 1: 数据下载与预处理
1. 从GEO下载series matrix文件
2. 探针注释转换（GPL570 → 基因符号）
3. 多探针合并（取均值）
4. 缺失值处理与标准化
5. 数据质控（PCA、聚类、箱线图）

### Phase 2: 差异表达分析
1. 正常 vs 脓毒症 差异分析
2. 脓毒症休克 vs 非休克 差异分析
3. HLA-DRA表达验证

### Phase 3: 免疫瘫痪评分(IPS)构建
```
IPS = w1 × HLA-DRA + w2 × HLA-DQB1 + w3 × CIITA + ...
```

### Phase 4: 临床验证
1. IPS与临床预后关联
2. 独立数据集验证
3. ROC曲线分析

## 使用方法

### 1. 数据下载与预处理
```bash
cd scripts
python3 download_and_preprocess.py
```

### 2. 后续分析
```R
# 在RStudio中运行
source("differential_expression.R")
source("build_ips_score.R")
```

## 预期结果

### MHC基因检测
| 数据集 | HLA-DRA | HLA-DQB1 | HLA-DQA1 | CIITA |
|--------|---------|----------|----------|-------|
| GSE26378 | ✓ | ✓ | ✓ | ✓ |
| GSE26440 | ✓ | ✓ | ✓ | ✓ |
| GSE13904 | ✓ | ✓ | ✓ | ✓ |

### 质控指标
- 样本离群检测
- 数据分布检验
- 批次效应评估

## 注意事项

1. **网络问题**: GEO数据下载可能较慢，建议使用代理或分批下载
2. **GPL570注释**: 该平台使用HG-U133_Plus_2芯片，需使用最新注释文件
3. **样本量**: 正常对照样本较少，可能影响差异分析效能

## 参考资料

- GEO: https://www.ncbi.nlm.nih.gov/geo/
- Platform GPL570: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
- IMMPORT: https://www.immport.org/

## 更新日志

- 2024-XX-XX: 项目初始化，完成目录结构和下载脚本
