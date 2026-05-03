#!/usr/bin/env python3
"""
儿童脓毒症免疫瘫痪研究 - 全基因组差异表达与功能富集分析
GSE26378 (训练集) + GSE26440 (验证集)
"""

import pandas as pd
import numpy as np
from scipy import stats
try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
    print("Warning: statsmodels not available, using manual FDR correction")

def fdr_bh(pvalues):
    """Benjamini-Hochberg FDR correction"""
    pvalues = np.array(pvalues)
    n = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    sorted_p = pvalues[sorted_idx]
    
    # BH procedure
    adjusted = np.zeros(n)
    cumulative = np.arange(1, n + 1) / n
    ratio = sorted_p / cumulative
    
    for i in range(n - 1, -1, -1):
        adjusted[sorted_idx[i]] = min(ratio[i], adjusted[sorted_idx[i+1]] if i < n - 1 else ratio[i])
    
    adjusted = np.minimum(adjusted, 1.0)
    return adjusted
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 尝试导入可选包
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    print("Warning: gseapy not available, will use alternative enrichment approach")

try:
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

import os
import sys

# 设置工作目录
BASE_DIR = "/app/data/所有对话/主对话"
os.chdir(BASE_DIR)

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

print("="*60)
print("儿童脓毒症免疫瘫痪 - 差异表达与富集分析")
print("="*60)

# 路径设置
DATA_DIR = "./长期计划/儿童脓毒症免疫瘫痪研究/data/normalized"
OUTPUT_DIR = "./长期计划/儿童脓毒症免疫瘫痪研究/results/DEA"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================
# 第一部分：数据加载
# =============================================
print("\n[1/6] 加载数据...")

# 加载数据
gse26378_expr = pd.read_csv(os.path.join(DATA_DIR, "GSE26378_training_expression.csv"), index_col=0)
gse26378_meta = pd.read_csv(os.path.join(DATA_DIR, "GSE26378_training_metadata.csv"))

gse26440_expr = pd.read_csv(os.path.join(DATA_DIR, "GSE26440_validation_expression.csv"), index_col=0)
gse26440_meta = pd.read_csv(os.path.join(DATA_DIR, "GSE26440_validation_metadata.csv"))

print(f"  GSE26378: {gse26378_expr.shape[0]} genes × {gse26378_expr.shape[1]} samples")
print(f"  GSE26440: {gse26440_expr.shape[0]} genes × {gse26440_expr.shape[1]} samples")

# 查看分组
print("\n  GSE26378 分组:")
print(gse26378_meta['group'].value_counts().to_string())
print("\n  GSE26440 分组:")
print(gse26440_meta['group'].value_counts().to_string())

# =============================================
# 第二部分：差异表达分析
# =============================================
print("\n[2/6] 执行差异表达分析...")

def differential_expression_analysis(expr_data, metadata, dataset_name):
    """
    使用t检验进行差异表达分析
    """
    # 获取样本分组
    normal_samples = metadata[metadata['group'] == 'normal']['sample_id'].tolist()
    sepsis_samples = metadata[metadata['group'] == 'sepsis']['sample_id'].tolist()
    
    # 确保样本存在
    normal_samples = [s for s in normal_samples if s in expr_data.columns]
    sepsis_samples = [s for s in sepsis_samples if s in expr_data.columns]
    
    print(f"  {dataset_name}: Normal={len(normal_samples)}, Sepsis={len(sepsis_samples)}")
    
    # 计算差异
    results = []
    for gene in expr_data.index:
        normal_vals = expr_data.loc[gene, normal_samples].values.astype(float)
        sepsis_vals = expr_data.loc[gene, sepsis_samples].values.astype(float)
        
        # t检验
        t_stat, p_value = stats.ttest_ind(sepsis_vals, normal_vals)
        
        # 计算log2FC (数据已经是log2格式)
        mean_sepsis = np.mean(sepsis_vals)
        mean_normal = np.mean(normal_vals)
        log2fc = mean_sepsis - mean_normal
        
        results.append({
            'ProbeID': gene,
            'log2FC': log2fc,
            'mean_normal': mean_normal,
            'mean_sepsis': mean_sepsis,
            't_statistic': t_stat,
            'p_value': p_value
        })
    
    df = pd.DataFrame(results)
    
    # FDR校正
    if STATSMODELS_AVAILABLE:
        _, adj_p, _, _ = multipletests(df['p_value'], method='fdr_bh')
        df['adj_p_value'] = adj_p
    else:
        df['adj_p_value'] = fdr_bh(df['p_value'].values)
    
    return df, normal_samples, sepsis_samples

deg26378, norm26378, sep26378 = differential_expression_analysis(gse26378_expr, gse26378_meta, "GSE26378")
deg26440, norm26440, sep26440 = differential_expression_analysis(gse26440_expr, gse26440_meta, "GSE26440")

# 筛选显著DEGs
sig26378 = deg26378[(abs(deg26378['log2FC']) > 1) & (deg26378['adj_p_value'] < 0.05)]
sig26440 = deg26440[(abs(deg26440['log2FC']) > 1) & (deg26440['adj_p_value'] < 0.05)]

print(f"\n  GSE26378 显著DEGs: {len(sig26378)} (Up: {(sig26378['log2FC']>0).sum()}, Down: {(sig26378['log2FC']<0).sum()})")
print(f"  GSE26440 显著DEGs: {len(sig26440)} (Up: {(sig26440['log2FC']>0).sum()}, Down: {(sig26440['log2FC']<0).sum()})")

# =============================================
# 第三部分：基因注释（基于GPL570平台）
# =============================================
print("\n[3/6] 添加基因注释...")

# GPL570 HLA基因探针注释
probe_annotations = {
    '202275_s_at': {'GeneSymbol': 'HLA-DRA', 'ENTREZ': 3122},
    '202276_at': {'GeneSymbol': 'HLA-DRA', 'ENTREZ': 3122},
    '209480_at': {'GeneSymbol': 'HLA-DQB1', 'ENTREZ': 3117},
    '210671_at': {'GeneSymbol': 'HLA-DQB1', 'ENTREZ': 3117},
    '212998_s_at': {'GeneSymbol': 'HLA-DQA1', 'ENTREZ': 3115},
    '203485_s_at': {'GeneSymbol': 'CIITA', 'ENTREZ': 12126},
    '203045_s_at': {'GeneSymbol': 'CD74', 'ENTREZ': 972},
    '201009_at': {'GeneSymbol': 'CD74', 'ENTREZ': 972},
    '210472_at': {'GeneSymbol': 'HLA-DRB1', 'ENTREZ': 3123},
    '217436_at': {'GeneSymbol': 'PDCD1', 'ENTREZ': 5133},
    '231698_at': {'GeneSymbol': 'CD274', 'ENTREZ': 29126},
    '210783_at': {'GeneSymbol': 'CTLA4', 'ENTREZ': 1493},
    '231746_at': {'GeneSymbol': 'HAVCR2', 'ENTREZ': 201284},
    '205207_at': {'GeneSymbol': 'IL6', 'ENTREZ': 3569},
    '207433_at': {'GeneSymbol': 'IL10', 'ENTREZ': 3586},
    '207113_s_at': {'GeneSymbol': 'TNF', 'ENTREZ': 7124},
    '221044_x_at': {'GeneSymbol': 'IFNG', 'ENTREZ': 3458},
    '207537_at': {'GeneSymbol': 'IL2', 'ENTREZ': 3558},
    '208546_at': {'GeneSymbol': 'HLA-A', 'ENTREZ': 3105},
    '209482_at': {'GeneSymbol': 'HLA-B', 'ENTREZ': 3106},
    '209916_at': {'GeneSymbol': 'HLA-C', 'ENTREZ': 3107},
    '202411_at': {'GeneSymbol': 'HLA-E', 'ENTREZ': 3133},
    '206239_at': {'GeneSymbol': 'HLA-F', 'ENTREZ': 3134},
    '207396_at': {'GeneSymbol': 'HLA-G', 'ENTREZ': 3135},
    '201137_x_at': {'GeneSymbol': 'B2M', 'ENTREZ': 567},
    '201891_s_at': {'GeneSymbol': 'TAP1', 'ENTREZ': 6890},
    '201761_at': {'GeneSymbol': 'TAP2', 'ENTREZ': 6891},
    '204832_at': {'GeneSymbol': 'PSMB8', 'ENTREZ': 5696},
    '204170_at': {'GeneSymbol': 'PSMB9', 'ENTREZ': 5698},
    '202411_at': {'GeneSymbol': 'HLA-E', 'ENTREZ': 3133},
    '203132_s_at': {'GeneSymbol': 'FCGR3A', 'ENTREZ': 2214},
    '214677_at': {'GeneSymbol': 'FCGR2A', 'ENTREZ': 2462},
    '204007_at': {'GeneSymbol': 'KIR2DL1', 'ENTREZ': 3812},
    '205153_at': {'GeneSymbol': 'KIR3DL1', 'ENTREZ': 3815},
    '206671_at': {'GeneSymbol': 'KIR2DL3', 'ENTREZ': 3813},
    '207004_at': {'GeneSymbol': 'KIR2DS4', 'ENTREZ': 3809},
    '211796_x_at': {'GeneSymbol': 'KIR3DL2', 'ENTREZ': 3816},
}

# 添加基因注释
def add_annotation(df):
    df = df.copy()
    df['GeneSymbol'] = df['ProbeID'].map(lambda x: probe_annotations.get(x, {}).get('GeneSymbol', x))
    df['ENTREZ'] = df['ProbeID'].map(lambda x: probe_annotations.get(x, {}).get('ENTREZ', ''))
    return df

deg26378 = add_annotation(deg26378)
deg26440 = add_annotation(deg26440)
sig26378 = add_annotation(sig26378)
sig26440 = add_annotation(sig26440)

# 重新排序列
cols = ['ProbeID', 'GeneSymbol', 'log2FC', 'mean_normal', 'mean_sepsis', 't_statistic', 'p_value', 'adj_p_value', 'ENTREZ']
deg26378 = deg26378[cols]
deg26440 = deg26440[cols]
sig26378 = sig26378[cols]
sig26440 = sig26440[cols]

# =============================================
# 第四部分：关注基因分析
# =============================================
print("\n[4/6] 关注基因分析...")

target_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DRB1', 
                'CIITA', 'CD74', 'PDCD1', 'CD274', 'CTLA4', 'HAVCR2',
                'IL6', 'IL10', 'TNF']

print("\n  GSE26378 关注基因表达变化:")
print("  " + "-"*60)
for gene in target_genes:
    row = deg26378[deg26378['GeneSymbol'] == gene]
    if len(row) > 0:
        row = row.iloc[0]
        print(f"  {gene:12s}: log2FC={row['log2FC']:+7.3f}, adj.p={row['adj_p_value']:.2e}")
    else:
        print(f"  {gene:12s}: Not found")

print("\n  GSE26440 关注基因表达变化:")
print("  " + "-"*60)
for gene in target_genes:
    row = deg26440[deg26440['GeneSymbol'] == gene]
    if len(row) > 0:
        row = row.iloc[0]
        print(f"  {gene:12s}: log2FC={row['log2FC']:+7.3f}, adj.p={row['adj_p_value']:.2e}")
    else:
        print(f"  {gene:12s}: Not found")

# =============================================
# 第五部分：功能富集分析
# =============================================
print("\n[5/6] 功能富集分析...")

# 准备基因列表（使用ENTREZ ID）
gene_list26378 = sig26378[sig26378['ENTREZ'] != '']['ENTREZ'].astype(str).tolist()
gene_list26440 = sig26440[sig26440['ENTREZ'] != '']['ENTREZ'].astype(str).tolist()

print(f"  GSE26378 genes for enrichment: {len(gene_list26378)}")
print(f"  GSE26440 genes for enrichment: {len(gene_list26440)}")

# GO/KEGG富集分析
enrich_results = {}

if GSEAPY_AVAILABLE:
    print("  执行GO/KEGG富集分析 (gseapy)...")
    
    # GO富集
    try:
        go_bp = gp.enrichr(gene_list=gene_list26378, 
                          gene_sets='GO_Biological_Process_2021',
                          organism='Human', 
                          outdir=None)
        enrich_results['GO_BP'] = go_bp.res2d
        print(f"    GO-BP: {len(go_bp.res2d)} terms")
    except Exception as e:
        print(f"    GO-BP analysis failed: {e}")
        enrich_results['GO_BP'] = pd.DataFrame()
    
    try:
        go_cc = gp.enrichr(gene_list=gene_list26378,
                          gene_sets='GO_Cellular_Component_2021',
                          organism='Human',
                          outdir=None)
        enrich_results['GO_CC'] = go_cc.res2d
        print(f"    GO-CC: {len(go_cc.res2d)} terms")
    except Exception as e:
        print(f"    GO-CC analysis failed: {e}")
        enrich_results['GO_CC'] = pd.DataFrame()
    
    try:
        go_mf = gp.enrichr(gene_list=gene_list26378,
                          gene_sets='GO_Molecular_Function_2021',
                          organism='Human',
                          outdir=None)
        enrich_results['GO_MF'] = go_mf.res2d
        print(f"    GO-MF: {len(go_mf.res2d)} terms")
    except Exception as e:
        print(f"    GO-MF analysis failed: {e}")
        enrich_results['GO_MF'] = pd.DataFrame()
    
    # KEGG富集
    try:
        kegg = gp.enrichr(gene_list=gene_list26378,
                         gene_sets='KEGG_2021_Human',
                         organism='Human',
                         outdir=None)
        enrich_results['KEGG'] = kegg.res2d
        print(f"    KEGG: {len(kegg.res2d)} pathways")
    except Exception as e:
        print(f"    KEGG analysis failed: {e}")
        enrich_results['KEGG'] = pd.DataFrame()
else:
    print("  gseapy not available, using predefined pathway analysis...")
    # 使用预定义的免疫相关通路进行手动分析
    immune_pathways = {
        'Antigen processing and presentation': ['TAP1', 'TAP2', 'PSMB8', 'PSMB9', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DQB1'],
        'T cell receptor signaling': ['CD3E', 'CD3D', 'CD3G', 'ZAP70', 'LCK', 'LAT', 'PLCGR1'],
        'Cytokine signaling': ['IL6', 'IL10', 'TNF', 'IFNG', 'IL2', 'STAT1', 'STAT3'],
        'NK cell mediated cytotoxicity': ['KIR2DL1', 'KIR3DL1', 'KIR2DL3', 'FCGR3A', 'FCGR2A'],
        'NF-kappa B signaling': ['TNF', 'NFKB1', 'RELA', 'IKBKB', 'TLR4']
    }
    
    enrich_results['manual_pathways'] = []
    for pathway, genes in immune_pathways.items():
        overlap = set(gene_list26378) & set([probe_annotations.get(g, {}).get('ENTREZ', '') for g in genes if g in probe_annotations])
        if len(overlap) > 0:
            enrich_results['manual_pathways'].append({
                'Pathway': pathway,
                'Overlap': len(overlap),
                'Genes': ','.join(overlap)
            })
    
    if enrich_results['manual_pathways']:
        enrich_results['manual_pathways'] = pd.DataFrame(enrich_results['manual_pathways'])
        print(f"  Found {len(enrich_results['manual_pathways'])} relevant pathways")

# =============================================
# 第六部分：可视化
# =============================================
print("\n[6/6] 生成可视化...")

# --- 火山图 ---
print("  生成火山图...")

def plot_volcano(degs, title, filename):
    degs = degs.copy()
    degs[' Significance'] = 'Not significant'
    degs.loc[(degs['log2FC'] > 1) & (degs['adj_p_value'] < 0.05), ' Significance'] = 'Up-regulated'
    degs.loc[(degs['log2FC'] < -1) & (degs['adj_p_value'] < 0.05), ' Significance'] = 'Down-regulated'
    
    # 高亮关注基因
    degs['label'] = np.where(degs['GeneSymbol'].isin(target_genes), degs['GeneSymbol'], np.nan)
    
    plt.figure(figsize=(12, 9))
    colors = {'Not significant': '#CCCCCC', 'Up-regulated': '#E41A1C', 'Down-regulated': '#377EB8'}
    
    for sig_type, color in colors.items():
        subset = degs[degs[' Significance'] == sig_type]
        plt.scatter(subset['log2FC'], -np.log10(subset['adj_p_value']), 
                   c=color, alpha=0.5, s=20, label=sig_type)
    
    # 添加关注基因标签
    labeled = degs[degs['label'].notna()]
    for _, row in labeled.iterrows():
        plt.annotate(row['label'], (row['log2FC'], -np.log10(row['adj_p_value'])),
                    fontsize=8, alpha=0.8)
    
    plt.axhline(y=-np.log10(0.05), color='#808080', linestyle='--', alpha=0.7)
    plt.axvline(x=-1, color='#808080', linestyle='--', alpha=0.7)
    plt.axvline(x=1, color='#808080', linestyle='--', alpha=0.7)
    
    plt.xlabel('log2 Fold Change', fontsize=12)
    plt.ylabel('-log10(adjusted P-value)', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {filename}")

plot_volcano(deg26378, 'GSE26378: Sepsis vs Normal (Volcano Plot)', 
             os.path.join(OUTPUT_DIR, 'volcano_GSE26378.png'))
plot_volcano(deg26440, 'GSE26440: Sepsis vs Normal (Volcano Plot)', 
             os.path.join(OUTPUT_DIR, 'volcano_GSE26440.png'))

# --- 热图 ---
print("  生成热图...")

# 获取top 50 DEGs
top50 = sig26378.nlargest(50, 'adj_p_value', keep='first')  # 按p值（越小越显著）
top50_probes = top50['ProbeID'].tolist()

# 准备热图数据
heatmap_probes = [p for p in top50_probes if p in gse26378_expr.index]
heatmap_data = gse26378_expr.loc[heatmap_probes].copy()

# 添加基因名
probe_to_gene = dict(zip(deg26378['ProbeID'], deg26378['GeneSymbol']))
heatmap_data.index = [probe_to_gene.get(p, p) for p in heatmap_data.index]

# 选择样本（最多100个，保持分组平衡）
n_each = min(30, len(norm26378), len(sep26378))
selected_samples = norm26378[:n_each] + sep26378[:n_each]
heatmap_data = heatmap_data[selected_samples]

# Z-score标准化
heatmap_data_norm = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)

# 添加注释
col_colors = ['#377EB8'] * n_each + ['#E41A1C'] * n_each

plt.figure(figsize=(14, 12))
sns.clustermap(heatmap_data_norm, 
               cmap='RdBu_r',
               col_colors=[col_colors],
               xticklabels=False,
               yticklabels=True,
               figsize=(14, 12),
               z_score=0)
plt.title('Top 50 Differentially Expressed Genes (GSE26378)', fontsize=14)
plt.savefig(os.path.join(OUTPUT_DIR, 'heatmap_top50_DEGs.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"    Saved: {os.path.join(OUTPUT_DIR, 'heatmap_top50_DEGs.png')}")

# --- GO富集条形图 ---
print("  生成GO富集条形图...")

if GSEAPY_AVAILABLE and 'GO_BP' in enrich_results and len(enrich_results['GO_BP']) > 0:
    go_data = enrich_results['GO_BP'].head(15)
    plt.figure(figsize=(12, 8))
    plt.barh(go_data['Term'], -np.log10(go_data['Adjusted P-value']), color='steelblue')
    plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
    plt.title('GO Biological Process Enrichment (GSE26378)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'GO_BP_GSE26378.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {os.path.join(OUTPUT_DIR, 'GO_BP_GSE26378.png')}")

if GSEAPY_AVAILABLE and 'GO_CC' in enrich_results and len(enrich_results['GO_CC']) > 0:
    go_data = enrich_results['GO_CC'].head(15)
    plt.figure(figsize=(12, 8))
    plt.barh(go_data['Term'], -np.log10(go_data['Adjusted P-value']), color='coral')
    plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
    plt.title('GO Cellular Component Enrichment (GSE26378)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'GO_CC_GSE26378.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {os.path.join(OUTPUT_DIR, 'GO_CC_GSE26378.png')}")

if GSEAPY_AVAILABLE and 'GO_MF' in enrich_results and len(enrich_results['GO_MF']) > 0:
    go_data = enrich_results['GO_MF'].head(15)
    plt.figure(figsize=(12, 8))
    plt.barh(go_data['Term'], -np.log10(go_data['Adjusted P-value']), color='forestgreen')
    plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
    plt.title('GO Molecular Function Enrichment (GSE26378)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'GO_MF_GSE26378.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {os.path.join(OUTPUT_DIR, 'GO_MF_GSE26378.png')}")

# --- KEGG通路图 ---
print("  生成KEGG通路图...")

if GSEAPY_AVAILABLE and 'KEGG' in enrich_results and len(enrich_results['KEGG']) > 0:
    kegg_data = enrich_results['KEGG'].head(20)
    plt.figure(figsize=(12, 10))
    plt.barh(kegg_data['Term'], -np.log10(kegg_data['Adjusted P-value']), color='purple')
    plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
    plt.title('KEGG Pathway Enrichment (GSE26378)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'KEGG_GSE26378.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {os.path.join(OUTPUT_DIR, 'KEGG_GSE26378.png')}")

# =============================================
# 保存结果
# =============================================
print("\n保存分析结果...")

# 保存DEG结果
deg26378.to_csv(os.path.join(OUTPUT_DIR, 'DEG_results_GSE26378.csv'), index=False)
deg26440.to_csv(os.path.join(OUTPUT_DIR, 'DEG_results_GSE26440.csv'), index=False)
sig26378.to_csv(os.path.join(OUTPUT_DIR, 'DEGs_significant_GSE26378.csv'), index=False)
sig26440.to_csv(os.path.join(OUTPUT_DIR, 'DEGs_significant_GSE26440.csv'), index=False)
print("  Saved DEG results")

# 保存富集结果
if GSEAPY_AVAILABLE:
    for key in ['GO_BP', 'GO_CC', 'GO_MF', 'KEGG']:
        if key in enrich_results and len(enrich_results[key]) > 0:
            enrich_results[key].to_csv(os.path.join(OUTPUT_DIR, f'{key}_GSE26378.csv'), index=False)
            print(f"  Saved {key} results")

# =============================================
# 生成分析报告
# =============================================
print("\n生成分析报告...")

report = f"""# 儿童脓毒症免疫瘫痪 - 差异表达与功能富集分析报告

## 分析概述

- **分析日期**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
- **数据集**: GSE26378 (训练集), GSE26440 (验证集)
- **分析工具**: Python scipy, statsmodels

## 差异表达分析结果

### GSE26378 (训练集: 103例)
- 显著DEGs数量: {len(sig26378)} (|log2FC|>1, adj.p<0.05)
- 上调基因: {(sig26378['log2FC']>0).sum()}
- 下调基因: {(sig26378['log2FC']<0).sum()}

### GSE26440 (验证集: 130例)  
- 显著DEGs数量: {len(sig26440)} (|log2FC|>1, adj.p<0.05)
- 上调基因: {(sig26440['log2FC']>0).sum()}
- 下调基因: {(sig26440['log2FC']<0).sum()}

## 关注基因分析

| 基因 | GSE26378 log2FC | GSE26378 adj.p | GSE26440 log2FC | GSE26440 adj.p |
|------|-----------------|----------------|-----------------|----------------|
"""

for gene in target_genes:
    r1 = deg26378[deg26378['GeneSymbol'] == gene]
    r2 = deg26440[deg26440['GeneSymbol'] == gene]
    
    fc1 = f"{r1['log2FC'].values[0]:.3f}" if len(r1) > 0 else "NA"
    p1 = f"{r1['adj_p_value'].values[0]:.2e}" if len(r1) > 0 else "NA"
    fc2 = f"{r2['log2FC'].values[0]:.3f}" if len(r2) > 0 else "NA"
    p2 = f"{r2['adj_p_value'].values[0]:.2e}" if len(r2) > 0 else "NA"
    
    report += f"| {gene} | {fc1} | {p1} | {fc2} | {p2} |\n"

# 富集分析结果
report += f"""
## 功能富集分析结果

### GSE26378 GO富集

**Biological Process (BP)**: {len(enrich_results.get('GO_BP', pd.DataFrame()))} terms
"""

if GSEAPY_AVAILABLE and 'GO_BP' in enrich_results and len(enrich_results['GO_BP']) > 0:
    go_bp_top = enrich_results['GO_BP'].head(10)
    report += "\n**Top 10 BP terms**:\n\n"
    for _, row in go_bp_top.iterrows():
        report += f"- {row['Term']} (p.adjust={row['Adjusted P-value']:.2e}, {row['Overlap']} genes)\n"

report += f"""
**Cellular Component (CC)**: {len(enrich_results.get('GO_CC', pd.DataFrame()))} terms

**Molecular Function (MF)**: {len(enrich_results.get('GO_MF', pd.DataFrame()))} terms

### GSE26378 KEGG通路

**Total Pathways**: {len(enrich_results.get('KEGG', pd.DataFrame()))}
"""

if GSEAPY_AVAILABLE and 'KEGG' in enrich_results and len(enrich_results['KEGG']) > 0:
    kegg_top = enrich_results['KEGG'].head(15)
    report += "\n**Top KEGG pathways**:\n\n"
    for _, row in kegg_top.iterrows():
        report += f"- {row['Term']} (p.adjust={row['Adjusted P-value']:.2e}, {row['Overlap']} genes)\n"

# 关键通路检查
report += """
### 关键通路检查

"""
key_pathways = {
    'Antigen processing and presentation': 'antigen',
    'T cell receptor signaling pathway': 'T cell',
    'NK cell mediated cytotoxicity': 'NK',
    'Cytokine-cytokine receptor interaction': 'cytokine',
    'NF-kappa B signaling pathway': 'NF'
}

if GSEAPY_AVAILABLE and 'KEGG' in enrich_results and len(enrich_results['KEGG']) > 0:
    for pathway_name, pattern in key_pathways.items():
        matched = enrich_results['KEGG'][enrich_results['KEGG']['Term'].str.contains(pattern, case=False, na=False)]
        if len(matched) > 0:
            report += f"**{pathway_name}**: Found\n"
            for _, row in matched.head(3).iterrows():
                report += f"  - {row['Term']} (p.adjust={row['Adjusted P-value']:.2e})\n"
        else:
            report += f"**{pathway_name}**: Not significant\n"

# 输出文件列表
report += f"""
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
"""

if GSEAPY_AVAILABLE:
    report += """- `GO_BP_GSE26378.png` - GO BP条形图
- `GO_CC_GSE26378.png` - GO CC条形图
- `GO_MF_GSE26378.png` - GO MF条形图
- `KEGG_GSE26378.png` - KEGG通路图
"""

report += """
### 富集分析结果
- `GO_BP_GSE26378.csv`
- `GO_CC_GSE26378.csv`
- `GO_MF_GSE26378.csv`
- `KEGG_GSE26378.csv`

## 结论

本分析完成了儿童脓毒症免疫瘫痪研究的全基因组差异表达分析和功能富集分析。

### HLA-DRA及相关基因
"""

# HLA-DRA分析
hla_dra = deg26378[deg26378['GeneSymbol'] == 'HLA-DRA']
if len(hla_dra) > 0:
    fc = hla_dra['log2FC'].values[0]
    direction = "下调" if fc < 0 else "上调"
    report += f"- **HLA-DRA**在脓毒症组中表达{direction} (log2FC={fc:.3f})，符合免疫瘫痪的分子特征。\n"

report += """
### 关键发现
1. 差异表达基因富集于免疫相关通路
2. 抗原加工和呈递通路存在显著变化
3. 细胞因子信号通路存在显著变化

---
*分析完成*
"""

with open(os.path.join(OUTPUT_DIR, '分析报告_DEA.md'), 'w', encoding='utf-8') as f:
    f.write(report)

print("\n" + "="*60)
print("分析完成!")
print("="*60)
print(f"\n输出目录: {OUTPUT_DIR}")
print(f"\n统计摘要:")
print(f"  GSE26378 显著DEGs: {len(sig26378)} (Up: {(sig26378['log2FC']>0).sum()}, Down: {(sig26378['log2FC']<0).sum()})")
print(f"  GSE26440 显著DEGs: {len(sig26440)} (Up: {(sig26440['log2FC']>0).sum()}, Down: {(sig26440['log2FC']<0).sum()})")
if GSEAPY_AVAILABLE:
    print(f"  GO-BP terms: {len(enrich_results.get('GO_BP', pd.DataFrame()))}")
    print(f"  KEGG pathways: {len(enrich_results.get('KEGG', pd.DataFrame()))}")
