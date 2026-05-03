# 模块4: GSEA通路富集分析框架

## 4.1 分析方法

### 推荐工具
1. **clusterProfiler** - GSEA+ORA分析
2. **fgsea** - 快速GSEA
3. **WebGestalt** - 在线富集分析

## 4.2 基因集选择

### Hallmark基因集
- HALLMARK_INTERFERON_GAMMA_RESPONSE
- HALLMARK_INFLAMMATORY_RESPONSE
- HALLMARK_ALLOGRAFT_REJECTION
- HALLMARK_IL6_JAK_STAT3_SIGNALING
- HALLMARK_COMPLEMENT

### C7免疫基因集
- C7.IMMUNESIGDB.v7.5.1.symbols.gmt

## 4.3 预期结果

### 儿童 (Subclass A vs B+C)
**下调通路**:
- INTERFERON_GAMMA_RESPONSE ↓↓
- INFLAMMATORY_RESPONSE ↓
- T_CELL_RECEPTOR_SIGNALING ↓
- ANTIGEN_PRESENTATION ↓

**上调通路**:
- APOPTOSIS ↑
- P53_PATHWAY ↑

### 成人 (Mars1 vs Mars2/3/4)
**下调通路**:
- MHC_II_ANTIGEN_PRESENTATION ↓

**上调通路**:
- INTERFERON_GAMMA_RESPONSE ↑ (代偿)
- INFLAMMATORY_RESPONSE ↑

## 4.4 CIITA通路网络

```
IFN-γ → STAT1 → IRF1 ──┐
                       ├──→ CIITA → MHC II genes
TLR → NF-κB ───────────┘
```

### 各通路富集预期

| 通路 | 儿童 | 成人 | 解释 |
|------|------|------|------|
| IFN-γ response | ↓↓ | ↑ | 炎症背景相反 |
| MHC II presentation | ↓↓ | ↓ | 功能均受损 |
| T cell activation | ↓ | ↓ | 适应性免疫抑制 |

## 4.5 分析代码框架

```r
library(clusterProfiler)
library(msigdbr)

# 获取Hallmark基因集
msig_h <- msigdbr(species = "Human", category = "H")

# GSEA分析
gsea_result <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = msig_h[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.25
)
```

---

*此框架待完整数据后执行*
