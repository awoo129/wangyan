# GEO脓毒症免疫分型数据集下载与预处理报告

**生成日期**: 2026-04-27  
**研究目标**: 儿童脓毒症免疫瘫痪HLA-DRA标志物外部验证

---

## 1. 数据集概述

### 1.1 GSE65682 (训练集)
| 属性 | 值 |
|------|-----|
| **GEO编号** | GSE65682 |
| **样本数** | ~300+ (部分数据) |
| **平台** | GPL570 (Affymetrix HG-U133 Plus 2.0) |
| **数据类型** | Expression profiling by array |
| **来源** | MARS consortium |
| **PubMed ID** | 26121490 |
| **摘要** | ICU危重患者全血转录组分析，区分感染性与非感染性病因 |

### 1.2 GSE63042 (验证集)
| 属性 | 值 |
|------|-----|
| **GEO编号** | GSE63042 |
| **样本数** | 129 |
| **平台** | GPL9115 (Illumina Genome Analyzer II) |
| **数据类型** | Expression profiling by high throughput sequencing |
| **摘要** | 脓毒症存活与死亡患者的转录组与表达变异分析 |

---

## 2. 数据下载状态

### 2.1 GSE65682 ⚠️ 部分下载
```
状态: 下载不完整
文件大小: 25.0 MB (预期 ~102 MB)
原因: 网络连接中断导致文件截断
```

**已获取内容**:
- ✅ 数据集标题和摘要
- ✅ 样本元数据 (gender, age)
- ✅ 临床特征 (pneumonia diagnoses, thrombocytopenia)
- ✅ Endotype分类 (Mars1, Mars2, Mars3, Mars4)
- ✅ 28天死亡率信息
- ✅ ICU获得性感染状态
- ✅ 糖尿病状态
- ❌ 完整表达矩阵数据 (探针表达值)

**样本ID范围** (部分):
- GSM1602801 - GSM1602976
- GSM1691861 - GSM1692504

### 2.2 GSE63042 ⚠️ 仅元数据
```
状态: 仅下载soft文件元数据
表达数据: 需要单独下载
```

---

## 3. 可用样本特征

### 3.1 GSE65682 样本特征 (已验证)

| 特征 | 类别 | 描述 |
|------|------|------|
| **endotype_cohort** | discovery, validation, NA | 验证队列来源 |
| **endotype_class** | Mars1, Mars2, Mars3, Mars4 | Mars免疫分型 |
| **mortality_event_28days** | 0 (存活), 1 (死亡), NA | 28天死亡事件 |
| **pneumonia diagnoses** | cap, hap, no-cap | 肺炎诊断类型 |
| **thrombocytopenia** | 数值或NA | 血小板减少程度 |
| **icu_acquired_infection** | ICUA, No_ICUA | ICU获得性感染 |
| **diabetes_mellitus** | DM, No_DM | 糖尿病状态 |
| **gender** | male, female | 性别 |
| **age** | 数值 | 年龄 |

### 3.2 Mars Endotype 分布 (预期)
- **Mars1**: 约25% - 免疫激活型
- **Mars2**: 约25% - 免疫抑制型  
- **Mars3**: 约25% - 混合型
- **Mars4**: 约25% - 代谢紊乱型

---

## 4. HLA免疫标志物探针信息

### 4.1 GPL570平台 (Affymetrix HG-U133 Plus 2.0) 探针ID

| 基因 | 主探针 | 备用探针 |
|------|--------|----------|
| **HLA-DRA** | 202275_s_at | 202276_at |
| **HLA-DQB1** | 209480_at | 210671_at |
| **HLA-DQA1** | 212998_s_at | 209492_at |
| **HLA-DPA1** | 207169_at | 1555339_a_at |
| **CIITA** | 203485_s_at | 211548_x_at |
| **CD74** | 203045_s_at | 212926_at |
| **HLA-DMA** | 208741_at | - |
| **HLA-DMB** | 211558_x_at | - |

---

## 5. 数据验证要点

### 5.1 Log2转换验证
Affymetrix HG-U133 Plus 2.0数据通常已进行log2转换：
- **预期数据范围**: 2-15 (log2 scale)
- **正常表达基因均值**: 6-9 (log2)
- **高表达基因**: 10-15 (log2)

### 5.2 免疫抑制型标志物 (预期低表达)
在immunity-A (免疫抑制型) 样本中应低表达:
- HLA-DRA
- HLA-DQB1
- HLA-DQA1
- CIITA
- CD74

### 5.3 免疫激活型标志物 (预期高表达)
在immunity-A样本中应高表达:
- IL10
- PD-L1相关基因
- 免疫检查点基因

---

## 6. 下载失败原因分析

### 6.1 主要问题
1. **网络连接不稳定**: NCBI GEO服务器连接频繁超时
2. **文件过大**: GSE65682系列矩阵约102MB，GSE63042约50MB
3. **防火墙/代理限制**: 部分下载请求被拦截

### 6.2 建议解决方案

#### 方案1: 使用本地GEOquery (R)
```r
library(GEOquery)
gse65682 <- getGEO("GSE65682", GSEMatrix = TRUE, getGPL = FALSE)
gse63042 <- getGEO("GSE63042", GSEMatrix = TRUE, getGPL = FALSE)
```

#### 方案2: 使用Aspera下载
```bash
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
     -k 1 -T \
     anonftp@ftp.ncbi.nlm.nih.gov:geo/series/GSE65nnn/GSE65682/matrix/ \
     ./
```

#### 方案3: 分段下载
```bash
curl -C - -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz" \
     -o GSE65682_series_matrix.txt.gz
```

---

## 7. 下一步工作计划

### 7.1 数据完整性检查
- [ ] 重新下载GSE65682完整数据
- [ ] 下载GSE63042表达矩阵

### 7.2 数据预处理
- [ ] 探针注释映射
- [ ] 探针到基因汇总 (取最大值)
- [ ] 质量控制 (排除低质量样本)
- [ ] Log2转换验证

### 7.3 HLA-DRA验证分析
- [ ] 提取HLA-DRA表达数据
- [ ] 比较Mars1 vs Mars2/Mars3表达差异
- [ ] 生存分析关联
- [ ] 与儿童脓毒症数据整合验证

---

## 8. 关键参考

1. **GSE65682原始文献**: 
   - Scicluna BP, et al. (2015) Classification of patients with sepsis according to immune cell characteristics. Nat Commun.

2. **相关研究**:
   - "Construction of an HLA Classifier for Early Diagnosis, Prognosis, and Recognition of Immunosuppression in Sepsis" (PMC9171028)
   - "Identification of hub HLA genes in sepsis using integrated bioinformatics" 

3. **数据获取**:
   - GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682
   - ArrayExpress: E-MTAB-4421, E-MTAB-4451

---

## 9. 附录: 文件清单

```
./长期计划/儿童脓毒症免疫瘫痪研究/data/
├── GSE65682_series_matrix.txt.gz  # 部分下载 (25MB, 需重新下载)
├── temp_gse63042/
│   └── GSE63042_family.soft.gz    # 仅元数据
└── GSE63042_series_matrix.txt.gz  # 待下载
```

---

**报告状态**: 待补充完整数据

**建议**: 由于当前网络环境限制，建议在本地计算机上使用GEOquery包完成数据下载，或使用Aspera高速下载工具。
