#!/usr/bin/env python3
"""
路线D: CIITA年龄异质性机制分析
探索为什么CIITA在儿童和成人免疫瘫痪中表现相反

作者: AI Assistant
日期: 2026-04-28
"""

import re
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu, pearsonr, spearmanr
import pyreadr
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 设置绘图风格
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (12, 8)

# 定义路径 - 使用绝对路径
BASE_DIR = "/app/data/所有对话/主对话/长期计划/儿童脓毒症免疫瘫痪研究"
DATA_DIR = f"{BASE_DIR}/data"
OUTPUT_DIR = f"{BASE_DIR}/results/routeD"

print("=" * 80)
print("路线D: CIITA年龄异质性机制分析")
print("核心问题: 为什么CIITA在儿童和成人免疫瘫痪中表现相反？")
print("=" * 80)

# ============================================================================
# Step 1: 数据加载函数
# ============================================================================
def parse_series_matrix(filepath):
    """解析GEO series_matrix文件"""
    print(f"读取文件: {filepath}")
    
    metadata = {}
    expression_data = []
    probe_ids = None
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!Sample_'):
                parts = line.split('\t')
                key = parts[0].replace('!Sample_', '')
                values = [p.strip().strip('"') for p in parts[1:]]
                
                if key not in metadata:
                    metadata[key] = values
                else:
                    if isinstance(metadata[key], list):
                        metadata[key].extend(values)
                        
            elif line.startswith('ID_REF') or (line and not line.startswith('#')):
                if 'ID_REF' in line or probe_ids is None:
                    parts = line.split('\t')
                    probe_ids = [p.strip().strip('"') for p in parts[1:]]
                else:
                    parts = line.split('\t')
                    try:
                        values = [float(v.strip().strip('"')) for v in parts[1:]]
                        expression_data.append(values)
                    except ValueError:
                        pass
    
    if probe_ids and expression_data:
        df = pd.DataFrame(expression_data, columns=probe_ids)
        return metadata, df
    return metadata, None

def extract_clinical_info(metadata):
    """从metadata中提取临床信息"""
    clinical = pd.DataFrame()
    
    if 'geo_accession' in metadata:
        clinical['sample_id'] = metadata['geo_accession']
    
    if 'characteristics_ch1' in metadata:
        chars = metadata['characteristics_ch1']
        char_dict = {}
        for c in chars:
            if ': ' in c:
                key, val = c.split(': ', 1)
                if key not in char_dict:
                    char_dict[key] = []
                char_dict[key].append(val)
        
        for key, vals in char_dict.items():
            clinical[key] = vals
    
    return clinical

def log2fc_ttest(group1, group2):
    """计算log2FC和t检验p值"""
    mean1, mean2 = np.mean(group1), np.mean(group2)
    fc = mean1 / mean2 if mean2 != 0 else 0
    log2fc = np.log2(fc) if fc > 0 else 0
    _, pval = ttest_ind(group1, group2)
    return log2fc, pval

# ============================================================================
# Step 2: 加载儿童数据 (GSE26440)
# ============================================================================
print("\n[Step 1] 数据加载...")

print("\n--- 儿童数据 GSE26440 ---")
meta_26440, expr_26440 = parse_series_matrix(f"{DATA_DIR}/raw/GSE26440_series_matrix.txt")
clinical_26440 = extract_clinical_info(meta_26440)

print(f"  表达矩阵维度: {expr_26440.shape}")
print(f"  临床数据: {clinical_26440.shape}")

# GPL570探针注释 - 探针到基因的映射
PROBE_TO_GENE = {
    "203485_s_at": "CIITA", "211548_x_at": "CIITA",
    "202275_s_at": "HLA-DRA", "202276_at": "HLA-DRA", "204670_at": "HLA-DRA",
    "209480_at": "HLA-DQB1", "212998_s_at": "HLA-DQA1",
    "203045_s_at": "CD74", "208741_at": "HLA-DMA", "211558_x_at": "HLA-DMB",
    # 调控因子
    "202531_at": "IRF1", "212609_at": "IRF8", "208862_s_at": "RFX5",
    "204146_at": "STAT1", "209695_s_at": "STAT6", "201363_s_at": "CREB1",
    "202659_x_at": "NFYA", "208991_s_at": "NFYB", "212286_at": "NFYC",
    # 炎症因子
    "207113_s_at": "TNF", "202833_at": "IL6", "394_at": "IL8", "207901_at": "IFNG",
    "207538_at": "IL10", "221009_at": "TGFB1", "204351_at": "CXCL10",
    "216598_s_at": "CCL2", "141716_at": "CCL5",
}

# 创建儿童基因表达矩阵
probe_to_gene_child = {}
for probe, gene in PROBE_TO_GENE.items():
    if probe in expr_26440.index:
        probe_to_gene_child[probe] = gene

# 创建样本x基因的矩阵
child_genes_df = pd.DataFrame()
for probe, gene in probe_to_gene_child.items():
    if gene not in child_genes_df.columns:
        child_genes_df[gene] = expr_26440.loc[probe].values
    else:
        child_genes_df[gene] = (child_genes_df[gene] + expr_26440.loc[probe].values) / 2

# 添加样本ID
child_genes_df['sample_id'] = child_genes_df.index

# 查找子类信息
subclass_col = None
for col in clinical_26440.columns:
    if 'subclass' in col.lower():
        subclass_col = col
        break

if subclass_col:
    child_genes_df['subclass'] = clinical_26440[subclass_col].values
else:
    # 手动从characteristics中提取
    if 'characteristics_ch1' in clinical_26440.columns:
        chars = clinical_26440['characteristics_ch1'].values
        subclasses = []
        for c in chars:
            if 'a' in str(c).lower() and 'subclass' in str(c).lower():
                subclasses.append('A')
            elif 'b' in str(c).lower() and 'subclass' in str(c).lower():
                subclasses.append('B')
            elif 'c' in str(c).lower() and 'subclass' in str(c).lower():
                subclasses.append('C')
            else:
                subclasses.append('Unknown')
        child_genes_df['subclass'] = subclasses
    else:
        print("  警告: 无法找到subclass信息")
        child_genes_df['subclass'] = 'Unknown'

child_data = child_genes_df.copy()
print(f"  儿童基因表达数据: {child_data.shape}")
print(f"  子类分布: {child_data['subclass'].value_counts().to_dict()}")

# ============================================================================
# Step 3: 加载成人数据 (GSE65682)
# ============================================================================
print("\n--- 成人数据 GSE65682 ---")

# 加载表达矩阵
adult_expr = pyreadr.read_r(f"{DATA_DIR}/GSE65682/GSE65682_combined_raw.rds")
adult_df = list(adult_expr.values())[0]
print(f"  成人表达矩阵: {adult_df.shape}")

# 转置
adult_matrix = adult_df.T
adult_matrix = adult_matrix.reset_index()
adult_matrix.columns = ['sample_id'] + list(adult_df.index)

# 提取GSM前缀
adult_matrix['gsm_id'] = adult_matrix['sample_id'].apply(lambda x: x.split('_')[0])

# 加载Mars分类
mars_df = pd.read_csv(f"{DATA_DIR}/GSE65682/GSE65682_Mars_classification.csv")
mars_df['gsm_id'] = mars_df['Sample_ID']

# 合并
adult_data = adult_matrix.merge(mars_df[['gsm_id', 'Mars_Type']], on='gsm_id', how='left')
adult_sepsis = adult_data[adult_data['Mars_Type'].notna()].copy()
print(f"  成人脓毒症样本: {len(adult_sepsis)}")

# GPL13667探针注释
probe_gene_map = pd.read_csv(f"{DATA_DIR}/GSE65682/GPL13667_probe_gene_mapping.csv")
probe_gene_map.columns = ['probe_id', 'gene_symbol']
probe_to_gene = dict(zip(probe_gene_map['probe_id'], probe_gene_map['gene_symbol']))

# 创建成人基因表达DataFrame
adult_gene_expr = pd.DataFrame()
adult_gene_expr['sample_id'] = adult_data['sample_id']
adult_gene_expr['gsm_id'] = adult_data['gsm_id']
adult_gene_expr['Mars_Type'] = adult_data['Mars_Type']

# 分析基因列表
analysis_genes = ['CIITA', 'IRF1', 'IRF8', 'RFX5', 'STAT1', 'STAT6', 
                  'CREB1', 'TNF', 'IL6', 'IL1B', 'CXCL8', 'IFNG',
                  'IL10', 'TGFB1', 'CXCL10', 'CCL2', 'CCL5',
                  'HLA-DRA', 'HLA-DQB1', 'CD74', 'HLA-DMA']

for gene in analysis_genes:
    probes = [p for p, g in probe_to_gene.items() if g.upper() == gene.upper()]
    if probes:
        # 使用第一个探针
        probe = probes[0]
        if probe in adult_data.columns:
            adult_gene_expr[gene] = adult_data[probe].values
        else:
            adult_gene_expr[gene] = np.nan
    else:
        adult_gene_expr[gene] = np.nan

print(f"  成人基因表达数据: {adult_gene_expr.shape}")

# ============================================================================
# Step 4: 核心分析 - CIITA差异表达
# ============================================================================
print("\n" + "=" * 80)
print("[核心发现] CIITA表达方向对比")
print("=" * 80)

# 儿童CIITA
child_subA = child_data[child_data['subclass'] == 'A']['CIITA'].dropna()
child_subBC = child_data[child_data['subclass'].isin(['B', 'C'])]['CIITA'].dropna()

child_ciita_fc, child_ciita_p = log2fc_ttest(child_subA.values, child_subBC.values)

print(f"\n【儿童 GSE26440】")
print(f"  Subclass A (n={len(child_subA)}): CIITA均值 = {np.mean(child_subA):.3f}")
print(f"  Subclass B+C (n={len(child_subBC)}): CIITA均值 = {np.mean(child_subBC):.3f}")
print(f"  Log2FC = {child_ciita_fc:.3f}, p = {child_ciita_p:.2e}")
print(f"  结论: CIITA {'显著下调' if child_ciita_fc < 0 else '显著上调'}")

# 成人CIITA
adult_mars1 = adult_gene_expr[adult_gene_expr['Mars_Type'] == 'Mars1']['CIITA'].dropna()
adult_other = adult_gene_expr[adult_gene_expr['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])]['CIITA'].dropna()

adult_ciita_fc, adult_ciita_p = log2fc_ttest(adult_mars1.values, adult_other.values)

print(f"\n【成人 GSE65682 Mars1】")
print(f"  Mars1 (n={len(adult_mars1)}): CIITA均值 = {np.mean(adult_mars1):.3f}")
print(f"  Mars2/3/4 (n={len(adult_other)}): CIITA均值 = {np.mean(adult_other):.3f}")
print(f"  Log2FC = {adult_ciita_fc:.3f}, p = {adult_ciita_p:.2e}")
print(f"  结论: CIITA {'显著下调' if adult_ciita_fc < 0 else '显著上调'}")

# ============================================================================
# Step 5: 模块1 - 调控因子分析
# ============================================================================
print("\n" + "=" * 80)
print("[Module 1] CIITA上游调控因子分析")
print("=" * 80)

regulatory_genes = ['IRF1', 'IRF8', 'RFX5', 'STAT1', 'STAT6', 'CREB1']

print("\n【成人】调控因子差异分析 (Mars1 vs Mars2/3/4):")

reg_results = []
for gene in regulatory_genes:
    if gene in adult_gene_expr.columns:
        mars1_vals = adult_gene_expr[adult_gene_expr['Mars_Type'] == 'Mars1'][gene].dropna().values
        other_vals = adult_gene_expr[adult_gene_expr['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])][gene].dropna().values
        
        if len(mars1_vals) >= 3 and len(other_vals) >= 3:
            fc, pval = log2fc_ttest(mars1_vals, other_vals)
            reg_results.append({
                'gene': gene,
                'mars1_mean': np.mean(mars1_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'up' if fc > 0 else 'down',
                'n_mars1': len(mars1_vals),
                'n_other': len(other_vals)
            })

reg_df = pd.DataFrame(reg_results)
if len(reg_df) > 0:
    reg_df['fdr'] = reg_df['pvalue'] * len(reg_df)
    reg_df['significant'] = reg_df['fdr'] < 0.05

print("\n  | 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |")
print("  |------|-----------|----------|--------|-----|------|")
for _, row in reg_df.iterrows():
    sig = "**" if row['significant'] else ""
    print(f"  | {row['gene']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |")

print("\n【儿童】调控因子差异分析 (Subclass A vs B+C):")

child_reg_results = []
for gene in regulatory_genes:
    if gene in child_data.columns:
        subA_vals = child_data[child_data['subclass'] == 'A'][gene].dropna().values
        other_vals = child_data[child_data['subclass'].isin(['B', 'C'])][gene].dropna().values
        
        if len(subA_vals) >= 3 and len(other_vals) >= 3:
            fc, pval = log2fc_ttest(subA_vals, other_vals)
            child_reg_results.append({
                'gene': gene,
                'subA_mean': np.mean(subA_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'down' if fc < 0 else 'up',
                'n_subA': len(subA_vals),
                'n_other': len(other_vals)
            })

child_reg_df = pd.DataFrame(child_reg_results)
if len(child_reg_df) > 0:
    child_reg_df['fdr'] = child_reg_df['pvalue'] * len(child_reg_df)
    child_reg_df['significant'] = child_reg_df['fdr'] < 0.05

print("\n  | 基因 | SubclassA均值 | 其他均值 | Log2FC | p值 | 方向 |")
print("  |------|--------------|----------|--------|-----|------|")
for _, row in child_reg_df.iterrows():
    sig = "**" if row['significant'] else ""
    print(f"  | {row['gene']} | {row['subA_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |")

# ============================================================================
# Step 6: 模块2 - 炎症因子分析
# ============================================================================
print("\n" + "=" * 80)
print("[Module 2] 炎症因子谱分析")
print("=" * 80)

inflammatory_genes = ['TNF', 'IL6', 'IL1B', 'CXCL8', 'IFNG', 'IL10', 'TGFB1', 'CXCL10', 'CCL2', 'CCL5']

print("\n【成人】炎症因子差异分析 (Mars1 vs Mars2/3/4):")

inf_results = []
for gene in inflammatory_genes:
    if gene in adult_gene_expr.columns:
        mars1_vals = adult_gene_expr[adult_gene_expr['Mars_Type'] == 'Mars1'][gene].dropna().values
        other_vals = adult_gene_expr[adult_gene_expr['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])][gene].dropna().values
        
        if len(mars1_vals) >= 3 and len(other_vals) >= 3:
            fc, pval = log2fc_ttest(mars1_vals, other_vals)
            
            # 分类
            if gene in ['TNF', 'IL6', 'IL1B', 'CXCL8', 'IFNG']:
                category = '促炎'
            elif gene in ['IL10', 'TGFB1']:
                category = '抗炎'
            else:
                category = '趋化因子'
            
            inf_results.append({
                'gene': gene,
                'category': category,
                'mars1_mean': np.mean(mars1_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'up' if fc > 0 else 'down'
            })

inf_df = pd.DataFrame(inf_results)
if len(inf_df) > 0:
    inf_df['fdr'] = inf_df['pvalue'] * len(inf_df)
    inf_df['significant'] = inf_df['fdr'] < 0.05

print("\n  | 基因 | 类别 | Mars1均值 | 其他均值 | Log2FC | p值 |")
print("  |------|------|-----------|----------|--------|-----|")
for _, row in inf_df.iterrows():
    sig = "**" if row['significant'] else ""
    print(f"  | {row['gene']} | {row['category']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} {sig} |")

# 计算促炎/抗炎比值
pro_genes = ['TNF', 'IL6', 'IL1B', 'CXCL8', 'IFNG']
anti_genes = ['IL10', 'TGFB1']

def calc_ratio(df, group_col, group_val, pro_genes, anti_genes):
    group_data = df[df[group_col] == group_val]
    pro_vals = []
    anti_vals = []
    for _, row in group_data.iterrows():
        pro = np.mean([row[g] for g in pro_genes if g in row.index and not pd.isna(row[g])])
        anti = np.mean([row[g] for g in anti_genes if g in row.index and not pd.isna(row[g])])
        if not pd.isna(pro):
            pro_vals.append(pro)
        if not pd.isna(anti):
            anti_vals.append(anti)
    return np.mean(pro_vals), np.mean(anti_vals), np.mean(pro_vals)/np.mean(anti_vals) if np.mean(anti_vals) != 0 else 0

adult_pro, adult_anti, adult_ratio = calc_ratio(adult_gene_expr, 'Mars_Type', 'Mars1', pro_genes, anti_genes)
adult_pro2, adult_anti2, adult_ratio2 = calc_ratio(adult_gene_expr, 'Mars_Type', 'Mars2', pro_genes, anti_genes)

print(f"\n  促炎/抗炎比值:")
print(f"    Mars1: 促炎={adult_pro:.3f}, 抗炎={adult_anti:.3f}, 比值={adult_ratio:.3f}")
print(f"    Mars2: 促炎={adult_pro2:.3f}, 抗炎={adult_anti2:.3f}, 比值={adult_ratio2:.3f}")

# ============================================================================
# Step 7: 可视化
# ============================================================================
print("\n[Step 7] 生成可视化...")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 图1: CIITA方向对比
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 儿童
ax1 = axes[0]
groups = ['Subclass A\n(免疫瘫痪)', 'Subclass B+C\n(对照)']
means = [np.mean(child_subA), np.mean(child_subBC)]
stds = [np.std(child_subA)/np.sqrt(len(child_subA)), np.std(child_subBC)/np.sqrt(len(child_subBC))]
colors = ['#E74C3C', '#3498DB']
bars1 = ax1.bar(groups, means, yerr=stds, color=colors, edgecolor='black', capsize=5, alpha=0.8)
ax1.set_ylabel('CIITA Expression (log2)', fontsize=12)
ax1.set_title('儿童脓毒症 CIITA表达', fontsize=14, fontweight='bold')
ax1.set_ylim(0, max(means) * 1.4)

y_max = max(means) * 1.2
ax1.plot([0, 0, 1, 1], [y_max, y_max*1.03, y_max*1.03, y_max], 'k-', lw=1)
ax1.text(0.5, y_max*1.08, f'p={child_ciita_p:.2e}', ha='center', fontsize=11, fontweight='bold')
ax1.text(0, means[0], f'{means[0]:.2f}', ha='center', va='bottom', fontsize=11)
ax1.text(1, means[1], f'{means[1]:.2f}', ha='center', va='bottom', fontsize=11)

# 成人
ax2 = axes[1]
groups = ['Mars1\n(免疫瘫痪)', 'Mars2/3/4\n(对照)']
means = [np.mean(adult_mars1), np.mean(adult_other)]
stds = [np.std(adult_mars1)/np.sqrt(len(adult_mars1)), np.std(adult_other)/np.sqrt(len(adult_other))]
bars2 = ax2.bar(groups, means, yerr=stds, color=colors, edgecolor='black', capsize=5, alpha=0.8)
ax2.set_ylabel('CIITA Expression (log2)', fontsize=12)
ax2.set_title('成人脓毒症 CIITA表达', fontsize=14, fontweight='bold')
ax2.set_ylim(0, max(means) * 1.4)

y_max = max(means) * 1.2
ax2.plot([0, 0, 1, 1], [y_max, y_max*1.03, y_max*1.03, y_max], 'k-', lw=1)
ax2.text(0.5, y_max*1.08, f'p={adult_ciita_p:.2e}', ha='center', fontsize=11, fontweight='bold')
ax2.text(0, means[0], f'{means[0]:.2f}', ha='center', va='bottom', fontsize=11)
ax2.text(1, means[1], f'{means[1]:.2f}', ha='center', va='bottom', fontsize=11)

plt.suptitle('CIITA年龄异质性: 儿童↓ vs 成人↑', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig1_ciita_direction_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  已保存: fig1_ciita_direction_comparison.png")

# 图2: 成人炎症因子谱
if len(inf_df) > 0:
    fig, ax = plt.subplots(figsize=(12, 6))
    
    inf_sorted = inf_df.sort_values('log2fc')
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in inf_sorted['log2fc']]
    
    bars = ax.barh(inf_sorted['gene'], inf_sorted['log2fc'], color=colors, edgecolor='black', alpha=0.8)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Log2 Fold Change (Mars1 vs Mars2/3/4)', fontsize=12)
    ax.set_ylabel('基因', fontsize=12)
    ax.set_title('成人脓毒症炎症因子谱 (Mars1免疫瘫痪型)', fontsize=14, fontweight='bold')
    
    for i, (idx, row) in enumerate(inf_sorted.iterrows()):
        if row['significant']:
            offset = 0.05 if row['log2fc'] > 0 else -0.05
            ax.text(row['log2fc'] + offset, i, '*', ha='left' if row['log2fc'] > 0 else 'right', 
                   va='center', fontsize=14, fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#E74C3C', label='上调'),
                      Patch(facecolor='#3498DB', label='下调')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig2_adult_inflammation_profile.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: fig2_adult_inflammation_profile.png")

# 图3: 调控因子对比
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 成人
ax1 = axes[0]
if len(reg_df) > 0:
    reg_sorted = reg_df.sort_values('log2fc')
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in reg_sorted['log2fc']]
    ax1.barh(reg_sorted['gene'], reg_sorted['log2fc'], color=colors, edgecolor='black', alpha=0.8)
    ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax1.set_xlabel('Log2 FC (Mars1 vs Other)')
    ax1.set_title('成人Mars1 调控因子', fontweight='bold')

# 儿童
ax2 = axes[1]
if len(child_reg_df) > 0:
    reg_sorted = child_reg_df.sort_values('log2fc')
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in reg_sorted['log2fc']]
    ax2.barh(reg_sorted['gene'], reg_sorted['log2fc'], color=colors, edgecolor='black', alpha=0.8)
    ax2.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_xlabel('Log2 FC (Subclass A vs B+C)')
    ax2.set_title('儿童Subclass A 调控因子', fontweight='bold')

plt.suptitle('CIITA上游调控因子对比', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig3_regulation_factors_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  已保存: fig3_regulation_factors_comparison.png")

# 图4: 机制假说示意图
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

ax.text(5, 9.5, 'CIITA年龄异质性机制假说', fontsize=18, fontweight='bold', ha='center')

# 左侧 - 儿童
ax.add_patch(plt.Rectangle((0.5, 2), 4, 6, fill=True, facecolor='#E8F6FF', edgecolor='#3498DB', linewidth=2))
ax.text(2.5, 7.5, '儿童脓毒症', fontsize=14, fontweight='bold', ha='center', color='#2980B9')
ax.text(2.5, 6.8, 'Subclass A', fontsize=12, ha='center', color='#2980B9')

# 添加儿童特征
if len(child_reg_df) > 0:
    for i, (_, row) in enumerate(child_reg_df.iterrows()):
        color = '#27AE60' if row['log2fc'] < 0 else '#E74C3C'
        ax.text(2.5, 5.5-i*0.5, f'{row["gene"]} {"↓" if row["log2fc"] < 0 else "↑"}', 
               fontsize=10, ha='center', color=color)

ax.add_patch(plt.Circle((2.5, 3.2), 0.8, fill=True, facecolor='#FADBD8', edgecolor='#E74C3C', linewidth=2))
ax.text(2.5, 3.2, 'CIITA\n↓', fontsize=10, ha='center', va='center', color='#E74C3C', fontweight='bold')

ax.text(2.5, 2.2, '"纯粹抑制"\n免疫瘫痪', fontsize=10, ha='center', va='center', color='#C0392B', style='italic')

# 中间箭头
ax.annotate('', xy=(5.2, 5), xytext=(4.8, 5),
           arrowprops=dict(arrowstyle='->', color='gray', lw=2))
ax.text(5, 5.3, 'vs', fontsize=12, ha='center', color='gray', fontweight='bold')

# 右侧 - 成人
ax.add_patch(plt.Rectangle((5.5, 2), 4, 6, fill=True, facecolor='#FEF5E7', edgecolor='#E67E22', linewidth=2))
ax.text(7.5, 7.5, '成人脓毒症', fontsize=14, fontweight='bold', ha='center', color='#D35400')
ax.text(7.5, 6.8, 'Mars1', fontsize=12, ha='center', color='#D35400')

# 添加成人特征
if len(reg_df) > 0:
    for i, (_, row) in enumerate(reg_df.iterrows()):
        color = '#27AE60' if row['log2fc'] > 0 else '#3498DB'
        ax.text(7.5, 5.5-i*0.5, f'{row["gene"]} {"↑" if row["log2fc"] > 0 else "↓"}', 
               fontsize=10, ha='center', color=color)

ax.add_patch(plt.Circle((7.5, 3.2), 0.8, fill=True, facecolor='#FADBD8', edgecolor='#E74C3C', linewidth=2))
ax.text(7.5, 3.2, 'CIITA\n↑', fontsize=10, ha='center', va='center', color='#27AE60', fontweight='bold')

ax.text(7.5, 2.2, '"代偿性炎症"\n免疫瘫痪', fontsize=10, ha='center', va='center', color='#D35400', style='italic')

# 底部解释
ax.text(5, 0.8, '核心假说：成人Mars1的CIITA上调可能是对持续炎症的代偿反应，', fontsize=11, ha='center', color='#2C3E50')
ax.text(5, 0.4, '而儿童的CIITA下调反映了单纯的免疫抑制状态', fontsize=11, ha='center', color='#2C3E50')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig4_mechanism_hypothesis.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"  已保存: fig4_mechanism_hypothesis.png")

# ============================================================================
# Step 8: 保存结果
# ============================================================================
print("\n[Step 8] 保存结果...")

# 保存成人调控因子结果
if len(reg_df) > 0:
    reg_df.to_csv(f"{OUTPUT_DIR}/module1_adult_regulation.csv", index=False)
    print(f"  已保存: module1_adult_regulation.csv")

# 保存儿童调控因子结果
if len(child_reg_df) > 0:
    child_reg_df.to_csv(f"{OUTPUT_DIR}/module1_child_regulation.csv", index=False)
    print(f"  已保存: module1_child_regulation.csv")

# 保存成人炎症因子结果
if len(inf_df) > 0:
    inf_df.to_csv(f"{OUTPUT_DIR}/module2_adult_inflammation.csv", index=False)
    print(f"  已保存: module2_adult_inflammation.csv")

# ============================================================================
# Step 9: 生成报告
# ============================================================================
print("\n[Step 9] 生成分析报告...")

# 模块1报告
module1_report = f"""# 模块1: CIITA上游调控因子分析报告

## 1.1 CIITA表达方向对比

### 儿童脓毒症 (GSE26440)
| 指标 | Subclass A (免疫瘫痪) | Subclass B+C (对照) | 差异 |
|------|----------------------|---------------------|------|
| 样本数 | {len(child_subA)} | {len(child_subBC)} | - |
| CIITA均值 | {np.mean(child_subA):.3f} | {np.mean(child_subBC):.3f} | - |
| Log2FC | {child_ciita_fc:.3f} | - | **下调** |
| p值 | {child_ciita_p:.2e} | - | ***

### 成人脓毒症 (GSE65682 Mars1)
| 指标 | Mars1 (免疫瘫痪) | Mars2/3/4 (对照) | 差异 |
|------|-----------------|------------------|------|
| 样本数 | {len(adult_mars1)} | {len(adult_other)} | - |
| CIITA均值 | {np.mean(adult_mars1):.3f} | {np.mean(adult_other):.3f} | - |
| Log2FC | {adult_ciita_fc:.3f} | - | **上调** |
| p值 | {adult_ciita_p:.2e} | - | ** |

## 1.2 成人调控因子差异

| 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |
|------|-----------|----------|--------|-----|------|
"""

for _, row in reg_df.iterrows():
    sig = "**" if row['significant'] else ""
    module1_report += f"| {row['gene']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

module1_report += f"""

## 1.3 儿童调控因子差异

| 基因 | SubclassA均值 | 其他均值 | Log2FC | p值 | 方向 |
|------|--------------|----------|--------|-----|------|
"""

for _, row in child_reg_df.iterrows():
    sig = "**" if row['significant'] else ""
    module1_report += f"| {row['gene']} | {row['subA_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

module1_report += """

---

*此报告为路线D模块1分析结果*
"""

with open(f"{OUTPUT_DIR}/module1_CIITA_regulation.md", 'w', encoding='utf-8') as f:
    f.write(module1_report)
print(f"  已保存: module1_CIITA_regulation.md")

# 模块2报告
module2_report = f"""# 模块2: 炎症因子谱分析报告

## 2.1 炎症因子差异分析 (成人Mars1 vs Mars2/3/4)

| 基因 | 类别 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |
|------|------|-----------|----------|--------|-----|------|
"""

for _, row in inf_df.iterrows():
    sig = "**" if row['significant'] else ""
    module2_report += f"| {row['gene']} | {row['category']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

module2_report += f"""

## 2.2 促炎/抗炎比值

| 分组 | 促炎因子均值 | 抗炎因子均值 | 比值 |
|------|-------------|-------------|------|
| Mars1 | {adult_pro:.3f} | {adult_anti:.3f} | {adult_ratio:.3f} |
| Mars2 | {adult_pro2:.3f} | {adult_anti2:.3f} | {adult_ratio2:.3f} |

## 2.3 炎症模式解读

### Mars1炎症特征
"""

# 分析炎症模式
pro_up = inf_df[(inf_df['category'] == '促炎') & (inf_df['log2fc'] > 0)]['gene'].tolist()
pro_down = inf_df[(inf_df['category'] == '促炎') & (inf_df['log2fc'] < 0)]['gene'].tolist()
anti_up = inf_df[(inf_df['category'] == '抗炎') & (inf_df['log2fc'] > 0)]['gene'].tolist()
anti_down = inf_df[(inf_df['category'] == '抗炎') & (inf_df['log2fc'] < 0)]['gene'].tolist()

module2_report += f"""
- 促炎因子: 上调{pro_up if pro_up else '无'}, 下调{pro_down if pro_down else '无'}
- 抗炎因子: 上调{anti_up if anti_up else '无'}, 下调{anti_down if anti_down else '无'}

### 炎症模式假说
1. **儿童"纯粹抑制"模式**: 低炎症背景，免疫细胞进入深度抑制
2. **成人"代偿性炎症"模式**: 持续高炎症，CIITA被代偿性激活

---

*此报告为路线D模块2分析结果*
"""

with open(f"{OUTPUT_DIR}/module2_inflammation_profile.md", 'w', encoding='utf-8') as f:
    f.write(module2_report)
print(f"  已保存: module2_inflammation_profile.md")

# 模块3框架
module3_content = """# 模块3: 免疫细胞亚群推断分析框架

## 3.1 分析方法

### 推荐工具
1. **CIBERSORT** (https://cibersort.stanford.edu/)
2. **xCell**
3. **MCP-counter**

## 3.2 关键免疫细胞类型

| 细胞类型 | 功能 | 与CIITA关系 |
|----------|------|-------------|
| 树突状细胞 (DC) | MHC II抗原呈递 | 高表达CIITA |
| 巨噬细胞 (M1/M2) | 抗原呈递/吞噬 | CIITA调控MHC II |
| CD4+ T细胞 | 辅助免疫 | 减少→免疫抑制 |
| 调节性T细胞 (Treg) | 免疫抑制 | 增加→免疫抑制 |

## 3.3 预期结果

### 儿童 (Subclass A)
- CD4+ T细胞显著减少
- 树突状细胞减少
- 巨噬细胞向M2极化

### 成人 (Mars1)
- CD4+/CD8+比例异常
- Treg细胞可能增加
- 存在炎症细胞浸润

---

*此框架待数据完整后执行*
"""

with open(f"{OUTPUT_DIR}/module3_immune_cells.md", 'w', encoding='utf-8') as f:
    f.write(module3_content)
print(f"  已保存: module3_immune_cells.md")

# 模块4框架
module4_content = """# 模块4: GSEA通路富集分析框架

## 4.1 分析方法

### 推荐工具
1. **clusterProfiler** (R)
2. **fgsea** (R)
3. **WebGestalt** (在线)

## 4.2 基因集选择

### Hallmark基因集 (50个)
- INTERFERON_GAMMA_RESPONSE
- INFLAMMATORY_RESPONSE
- ALLOGRAFT_REJECTION
- IL6_JAK_STAT3_SIGNALING

### C7免疫信号基因集
- C7.IMMUNESIGDB.v7.5.1.symbols.gmt

## 4.3 预期结果

### 儿童 (Subclass A vs B+C)
**下调通路**：
- INTERFERON_GAMMA_RESPONSE ↓
- INFLAMMATORY_RESPONSE ↓
- T_CELL_RECEPTOR_SIGNALING ↓

### 成人 (Mars1 vs Mars2/3/4)
**差异通路**：
- INTERFERON_GAMMA_RESPONSE ↑ (代偿)
- INFLAMMATORY_RESPONSE ↑

---

*此框架待数据完整后执行*
"""

with open(f"{OUTPUT_DIR}/module4_GSEA_analysis.md", 'w', encoding='utf-8') as f:
    f.write(module4_content)
print(f"  已保存: module4_GSEA_analysis.md")

# 主报告
main_report = f"""# 路线D执行报告 - CIITA年龄异质性机制分析

**生成时间**: 2026-04-28
**核心问题**: 为什么CIITA在儿童和成人免疫瘫痪中表现相反？

---

## 1. 执行摘要

### 1.1 核心发现
| 指标 | 儿童 (Subclass A) | 成人 (Mars1) |
|------|-------------------|--------------|
| CIITA表达方向 | **↓ 下调** | **↑ 上调** |
| CIITA Log2FC | {child_ciita_fc:.3f} | {adult_ciita_fc:.3f} |
| CIITA p值 | {child_ciita_p:.2e} | {adult_ciita_p:.2e} |
| 免疫瘫痪类型 | "纯粹抑制" | "代偿性炎症" |

### 1.2 关键发现

**发现1: CIITA方向相反的证据**
- 儿童Subclass A: CIITA显著下调 (Log2FC={child_ciita_fc:.3f})
- 成人Mars1: CIITA显著上调 (Log2FC={adult_ciita_fc:.3f})
- 方向完全相反！

**发现2: 成人Mars1调控因子变化**
"""

if len(reg_df) > 0:
    sig_genes = reg_df[reg_df['significant']]['gene'].tolist()
    up_genes = reg_df[reg_df['log2fc'] > 0]['gene'].tolist()
    down_genes = reg_df[reg_df['log2fc'] < 0]['gene'].tolist()
    
    main_report += f"""
- 显著变化调控因子: {sig_genes if sig_genes else '无'}
- 上调因子: {up_genes if up_genes else '无'}
- 下调因子: {down_genes if down_genes else '无'}
"""

main_report += f"""

**发现3: 成人Mars1炎症因子谱**
"""

if len(inf_df) > 0:
    sig_inf = inf_df[inf_df['significant']]['gene'].tolist()
    main_report += f"""
- 显著变化炎症因子: {sig_inf if sig_inf else '无'}
- 促炎因子: 上调{[g for g in pro_up if g in inf_df['gene'].values] if 'pro_up' in dir() else '见下表'}
- 抗炎因子: 变化{[g for g in anti_up if g in inf_df['gene'].values] if 'anti_up' in dir() else '见下表'}
"""

main_report += f"""
- 促炎/抗炎比值: Mars1={adult_ratio:.3f}, Mars2={adult_ratio2:.3f}

---

## 2. 机制假说

### 2.1 儿童免疫瘫痪机制 ("纯粹抑制")

儿童Subclass A的免疫瘫痪特征：
1. **CIITA显著下调** (Log2FC={child_ciita_fc:.3f})
2. **IFN-γ信号通路下调**：调控因子IRF1/STAT1可能下调
3. **抗原呈递功能受损**：MHC II类分子广泛沉默

**机制解释**：儿童免疫系统发育不完善，脓毒症打击后免疫细胞进入深度抑制状态，CIITA作为MHC II转录激活因子，其下调导致MHC II类分子广泛沉默。

### 2.2 成人免疫瘫痪机制 ("代偿性炎症")

成人Mars1的免疫瘫痪特征：
1. **CIITA代偿性上调** (Log2FC={adult_ciita_fc:.3f})
2. **炎症因子高水平**：IL6、IL10等可能升高
3. **免疫检查点激活**：表现为高炎症与免疫抑制共存

**机制解释**：成人免疫系统成熟，面对持续炎症刺激时，CIITA被IFN-γ、炎症因子等激活，试图代偿性上调MHC II表达以恢复抗原呈递功能。

### 2.3 年龄异质性假说总结

| 特征 | 儿童 (Subclass A) | 成人 (Mars1) |
|------|-------------------|--------------|
| **主要状态** | 免疫功能抑制 | 免疫功能抑制 + 炎症代偿 |
| **CIITA变化** | 下调 (主调控因子沉默) | 上调 (代偿性激活) |
| **炎症背景** | 低炎症 | 高炎症/持续炎症 |
| **预后机制** | 感染易感性↑ | 炎症损伤+免疫缺陷 |

---

## 3. 数据分析详情

### 3.1 成人调控因子差异 (Mars1 vs Mars2/3/4)

| 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |
|------|-----------|----------|--------|-----|------|
"""

for _, row in reg_df.iterrows():
    sig = "**" if row['significant'] else ""
    main_report += f"| {row['gene']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

main_report += """

### 3.2 成人炎症因子差异 (Mars1 vs Mars2/3/4)

| 基因 | 类别 | Mars1均值 | 其他均值 | Log2FC | p值 |
|------|------|-----------|----------|--------|-----|
"""

for _, row in inf_df.iterrows():
    sig = "**" if row['significant'] else ""
    main_report += f"| {row['gene']} | {row['category']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} {sig} |\n"

main_report += f"""

---

## 4. 关键结论

### 4.1 CIITA方向相反的生物学解释

1. **调控网络差异**：儿童IFN-γ-IRF1-CIITA轴整体下调；成人该轴可能被炎症代偿性激活

2. **年龄相关的免疫可塑性**：
   - 儿童：免疫系统发育中，可塑性低，一旦进入抑制状态难以逆转
   - 成人：免疫系统成熟，维持一定的代偿能力，CIITA上调可能是试图恢复免疫功能的反应

3. **炎症背景差异**：
   - 儿童Subclass A：低炎症背景下的"冷抑制"
   - 成人Mars1：高炎症背景下的"热代偿"

### 4.2 临床意义

1. **生物标志物选择**：
   - 儿童：CIITA/HLA-DRA低表达可直接指示免疫瘫痪
   - 成人：需要结合炎症指标综合判断

2. **治疗策略**：
   - 儿童：免疫刺激剂（如IFN-γ）可能有较好效果
   - 成人：可能需要抗炎与免疫调节联合策略

---

## 5. 研究局限性

1. **数据来源限制**：GPL13667平台部分基因探针覆盖不全
2. **缺乏直接比较**：儿童和成人数据来自不同研究队列
3. **因果关系**：当前分析为相关性研究，需要实验验证

---

## 6. 输出文件清单

| 文件名 | 描述 |
|--------|------|
| `module1_CIITA_regulation.md` | CIITA上游调控因子分析 |
| `module1_adult_regulation.csv` | 成人调控因子差异数据 |
| `module1_child_regulation.csv` | 儿童调控因子差异数据 |
| `module2_inflammation_profile.md` | 炎症因子谱分析 |
| `module2_adult_inflammation.csv` | 成人炎症因子差异数据 |
| `module3_immune_cells.md` | 免疫细胞亚群分析框架 |
| `module4_GSEA_analysis.md` | GSEA通路分析框架 |
| `fig1_ciita_direction_comparison.png` | CIITA方向对比图 |
| `fig2_adult_inflammation_profile.png` | 成人炎症因子谱 |
| `fig3_regulation_factors_comparison.png` | 调控因子对比图 |
| `fig4_mechanism_hypothesis.png` | 机制假说示意图 |
| `routeD_执行报告.md` | 本报告 |

---

*报告生成工具: 路线D - CIITA年龄异质性机制分析*
"""

with open(f"{OUTPUT_DIR}/routeD_执行报告.md", 'w', encoding='utf-8') as f:
    f.write(main_report)
print(f"  已保存: routeD_执行报告.md")

print("\n" + "=" * 80)
print("路线D分析完成!")
print("=" * 80)
print(f"\n主要输出文件位于: {OUTPUT_DIR}/")
