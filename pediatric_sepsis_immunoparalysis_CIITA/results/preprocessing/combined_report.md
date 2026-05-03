# 三数据集验证分析报告

## 分析日期
2026-04-27

## GSE13904 数据集处理结果

### 样本信息
- 总样本数: 227
- Normal组: 18
- Sepsis组: 209

### MHC II类基因差异表达 (GSE13904)

| 基因 | Normal均值 | Sepsis均值 | Log2FC | P值 |
|------|-----------|------------|--------|-----|
| HLA-DRA | 1.048 | 1.203 | 0.155 ↑ | 9.70e-02 ns |
| HLA-DRB1 | 1.310 | 0.672 | -0.638 ↓ | 1.84e-04 *** |
| HLA-DQB1 | 1.599 | 0.585 | -1.014 ↓ | 2.06e-06 *** |
| HLA-DQA1 | 1.015 | 0.756 | -0.259 ↓ | 2.37e-04 *** |
| HLA-DPA1 | 0.957 | 1.077 | 0.121 ↑ | 2.83e-02 * |
| HLA-DPB1 | 1.018 | 1.101 | 0.083 ↑ | 9.09e-02 ns |
| CIITA | 1.231 | 1.057 | -0.174 ↓ | 1.98e-01 ns |
| CD74 | 1.011 | 1.049 | 0.037 ↑ | 7.16e-01 ns |

### HLA-DRA 关键发现
- **Log2FC**: 0.155
- **P值**: 9.70e-02
- **方向**: 上调 (不符合免疫瘫痪特征)

## 三数据集一致性验证

### 汇总对比表

| 数据集 | 样本量 | Normal | Sepsis | HLA-DRA Log2FC | P值 | 效应方向 |
|--------|--------|--------|--------|----------------|-----|----------|
| GSE26378 | 103 | 21 | 82 | -1.214 | 5.79e-03** | ↓下调 |
| GSE26440 | 130 | 32 | 98 | -1.735 | 1.44e-07*** | ↓下调 |
| GSE13904 | 227 | 18 | 209 | 0.155 | 9.70e-02ns | ↑上调 |

### 一致性检验
- **所有数据集效应方向一致**: ✗ 不一致
- **统计显著性**: 部分显著

## Meta分析结果

### 合并效应量
- **Pooled Log2FC**: -0.034
- **Pooled Z-score**: -0.390
- **Pooled P-value**: 6.96e-01

### 异质性评估
- **I²统计量**: 94.7%
- **Q统计量**: 37.886
- **研究数**: 3

### I²解读
- I² < 25%: 低异质性
- I² = 25-50%: 中等异质性  
- I² > 50%: 高异质性

**当前I² = 94.7%**: 高异质性，解读需谨慎

## 重要发现与分析

### GSE13904数据集的独特性

GSE13904数据集表现出与其他两个数据集显著不同的模式：

1. **HLA-DRA表达差异**: 
   - GSE13904中HLA-DRA在脓毒症组**上调** (Log2FC = +0.155)
   - 而GSE26378和GSE26440中均**下调**

2. **可能的原因**:
   - **样本构成差异**: GSE13904包含SIRS（全身炎症反应综合征）患者，而GSE26378/GSE26440仅包含脓毒症休克患者
   - **数据处理差异**: 不同数据集的预处理流程可能导致标准化差异
   - **生物学异质性**: 儿童脓毒症谱系的临床表现具有高度异质性

3. **MHC II类基因整体模式**:
   - HLA-DRB1, HLA-DQB1, HLA-DQA1在GSE13904中仍呈下调趋势
   - 这表明HLA-DRA可能是该数据集中的异常值

### 结论

1. **GSE13904验证结果**: HLA-DRA在脓毒症组上调(Log2FC = 0.15, p = 9.70e-02)

2. **三数据集一致性**: ✗ 不一致

3. **Meta分析**: 合并效应量不显著(pooled Log2FC = -0.034, p = 6.96e-01)

4. **高异质性**: I² = 94.7%，提示数据集间存在显著差异

## 生成的文件

### 数据文件
- `data/normalized/GSE13904_validation_expression.csv`
- `data/normalized/GSE13904_validation_metadata.csv`
- `data/normalized/GSE13904_validation_mhc_expression.csv`
- `data/normalized/GSE13904_validation_ips_scores.csv`
- `data/normalized/GSE13904_validation_de_results.csv`
- `data/normalized/three_datasets_summary.csv`
- `data/normalized/meta_analysis_results.csv`

### 可视化图表
- `results/visualization/three_datasets_validation.png`
