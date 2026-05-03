# 模块3: 免疫细胞亚群推断分析框架

## 3.1 分析方法

### 推荐工具
1. **CIBERSORT** - 22种免疫细胞比例估算
2. **xCell** - 64种细胞类型富集分析
3. **MCP-counter** - 8种免疫细胞+2种基质细胞

## 3.2 关键细胞类型

| 细胞类型 | 功能 | 儿童Subclass A | 成人Mars1 |
|----------|------|----------------|-----------|
| 树突状细胞 | MHC II呈递 | ↓↓ | ↓/代偿 |
| CD4+ T细胞 | 辅助免疫 | ↓↓ | ↓ |
| 巨噬细胞M1 | 促炎 | ↓ | ↑ |
| 巨噬细胞M2 | 抗炎/修复 | ↑ | ↑ |
| Treg细胞 | 免疫抑制 | ↑ | ↑ |

## 3.3 预期结果

### 儿童 (Subclass A)
- DC细胞显著减少
- CD4+ T细胞减少
- M2极化增加
- 与CIITA下调相关

### 成人 (Mars1)
- DC细胞可能部分保留
- 免疫细胞异质性增加
- 存在炎症细胞浸润
- CIITA上调可能与DC激活相关

## 3.4 分析代码框架

```r
# CIBERSORT分析
library(CIBERSORT)
result <- CIBERSORT(LM22, expr_matrix, QN=TRUE)

# 比较组间差异
immunoparalysis <- result[immunoparalysis_samples, ]
control <- result[control_samples, ]
```

---

*此框架待完整数据后执行*
