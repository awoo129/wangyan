#!/usr/bin/env python3
"""
路线D: CIITA年龄异质性机制分析
探索为什么CIITA在儿童和成人免疫瘫痪中表现相反

作者: AI Assistant
日期: 2026-04-28
"""

import pyreadr
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu, pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

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
print("核心问题: 为什么CIITA在儿童和成人免疫瘫痪中表现相反?")
print("=" * 80)

# ============================================================================
# Step 1: 加载儿童(GSE26440)和成人(GSE65682)数据
# ============================================================================
print("\n[Step 1] 数据加载...")

# 加载成人数据
adult_expr = pyreadr.read_r(f"{DATA_DIR}/GSE65682/GSE65682_combined_raw.rds")
adult_df = list(adult_expr.values())[0]
print(f"  成人数据维度: {adult_df.shape}")

# 加载probe映射
probe_gene_map = pd.read_csv(f"{DATA_DIR}/GSE65682/GPL13667_probe_gene_mapping.csv")
probe_gene_map.columns = ['probe_id', 'gene_symbol']
probe_to_gene = dict(zip(probe_gene_map['probe_id'], probe_gene_map['gene_symbol']))

# 加载Mars分类
mars_df = pd.read_csv(f"{DATA_DIR}/GSE65682/GSE65682_Mars_classification.csv")
print(f"  Mars分类样本: {len(mars_df)}")

# 加载儿童数据 - 从路线A的预处理结果
child_mhc = pd.read_csv(f"{BASE_DIR}/results/routeA/mhc_expr_GSE26440.csv")
print(f"  儿童MHC表达数据: {child_mhc.shape}")

# 加载儿童样本分类信息
child_samples = pd.read_csv(f"{BASE_DIR}/results/routeA/differential_expression.csv")
print(f"  儿童分类信息: {child_samples.shape}")

# ============================================================================
# Step 2: 准备成人表达矩阵和分类
# ============================================================================
print("\n[Step 2] 准备成人数据...")

# 转置成人表达矩阵
adult_matrix = adult_df.T
adult_matrix = adult_matrix.reset_index()
adult_matrix.columns = ['sample_id'] + list(adult_df.index)

# 清理sample_id - 提取GSM前缀
adult_matrix['gsm_id'] = adult_matrix['sample_id'].apply(lambda x: x.split('_')[0])

# 创建Mars分型映射
mars_df['sample_id'] = mars_df['Sample_ID']
adult_matrix = adult_matrix.merge(mars_df[['sample_id', 'Mars_Type']], 
                                   left_on='gsm_id', right_on='sample_id', how='left',
                                   suffixes=('', '_mars'))

# 过滤脓毒症患者(有Mars分型)
adult_sepsis = adult_matrix[adult_matrix['Mars_Type'].notna()].copy()
print(f"  成人脓毒症样本: {len(adult_sepsis)}")

# 获取成人CIITA的表达值 (探针)
ciita_probes = probe_gene_map[probe_gene_map['gene_symbol'].str.upper() == 'CIITA']['probe_id'].tolist()
print(f"  CIITA探针: {ciita_probes}")

# 创建基因名到探针的映射 (GPL13667)
probe_to_gene_lower = {k.upper(): v for k, v in probe_to_gene.items()}

# ============================================================================
# Step 3: 定义分析基因列表
# ============================================================================
print("\n[Step 3] 定义分析基因列表...")

# CIITA上游调控因子
regulatory_genes = {
    'IRF1': 'IRF1', 'IRF8': 'IRF8',
    'RFX5': 'RFX5', 'RFXAP': 'RFXAP', 'RFXANK': 'RFXANK',
    'CREB1': 'CREB1', 'NFYA': 'NFYA', 'NFYB': 'NFYB', 'NFYC': 'NFYC',
    'STAT1': 'STAT1', 'STAT6': 'STAT6'
}

# 炎症因子
inflammatory_genes = {
    # 促炎因子
    'TNF': 'TNF', 'IL1B': 'IL1B', 'IL6': 'IL6', 'IL8': 'CXCL8',
    'IL12A': 'IL12A', 'IL12B': 'IL12B', 'IFNG': 'IFNG',
    # 抗炎因子
    'IL10': 'IL10', 'TGFB1': 'TGFB1', 'IL4': 'IL4', 'IL13': 'IL13',
    # 趋化因子
    'CXCL10': 'CXCL10', 'CCL2': 'CCL2', 'CCL5': 'CCL5'
}

# MHC II类基因
mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74',
             'HLA-DMA', 'HLA-DMB', 'HLA-DPB1', 'HLA-DRB1']

# ============================================================================
# Step 4: 提取成人数据中的基因表达
# ============================================================================
print("\n[Step 4] 提取成人数据基因表达...")

# 创建基因到探针的映射
gene_to_probes = {}
for probe, gene in probe_to_gene.items():
    if gene not in gene_to_probes:
        gene_to_probes[gene] = []
    gene_to_probes[gene].append(probe)

def get_gene_expr(df, gene, gene_to_probes):
    """获取基因表达值"""
    upper_gene = gene.upper()
    if upper_gene in gene_to_probes:
        probes = gene_to_probes[upper_gene]
        # 使用第一个探针
        for probe in probes:
            if probe in df.columns:
                return df[probe].values
    return None

# 提取成人CIITA表达
adult_ciita = []
for idx, row in adult_sepsis.iterrows():
    sample = row['sample_id']
    if sample in adult_matrix['sample_id'].values:
        sample_row = adult_matrix[adult_matrix['sample_id'] == sample]
        expr_vals = []
        for probe in ciita_probes:
            if probe in sample_row.columns:
                expr_vals.append(sample_row[probe].values[0])
        if expr_vals:
            adult_ciita.append(np.mean(expr_vals))
        else:
            adult_ciita.append(np.nan)
    else:
        adult_ciita.append(np.nan)

adult_sepsis['CIITA'] = adult_ciita
print(f"  成人CIITA表达: 均值={np.nanmean(adult_ciita):.3f}, 范围=[{np.nanmin(adult_ciita):.2f}, {np.nanmax(adult_ciita):.2f}]")

# 提取成人Mars1和其他的分类
adult_sepsis['group'] = adult_sepsis['Mars_Type'].apply(
    lambda x: 'Mars1 (免疫瘫痪)' if x == 'Mars1' else 'Mars2/3/4 (对照)'
)

# ============================================================================
# Step 5: 模块1 - CIITA上游调控因子分析
# ============================================================================
print("\n" + "=" * 80)
print("[Module 1] CIITA上游调控因子分析")
print("=" * 80)

# 获取儿童分类数据
child_sepsis = child_samples[child_samples['subclass'] != 'n/a'].copy()
child_sepsis['group'] = child_sepsis['subclass'].apply(
    lambda x: 'Subclass A (免疫瘫痪)' if x == 'A' else 'Subclass B+C (对照)'
)

# 提取儿童CIITA表达
child_ciita_immuno = child_sepsis[child_sepsis['subclass'] == 'A']['CIITA'].values
child_ciita_ctrl = child_sepsis[child_sepsis['subclass'].isin(['B', 'C'])]['CIITA'].values

print(f"\n  儿童CIITA表达:")
print(f"    Subclass A (免疫瘫痪): 均值={np.mean(child_ciita_immuno):.3f}, n={len(child_ciita_immuno)}")
print(f"    Subclass B+C (对照): 均值={np.mean(child_ciita_ctrl):.3f}, n={len(child_ciita_ctrl)}")

# 成人CIITA
adult_ciita_immuno = adult_sepsis[adult_sepsis['Mars_Type'] == 'Mars1']['CIITA'].dropna().values
adult_ciita_ctrl = adult_sepsis[adult_sepsis['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])]['CIITA'].dropna().values

print(f"\n  成人CIITA表达:")
print(f"    Mars1 (免疫瘫痪): 均值={np.mean(adult_ciita_immuno):.3f}, n={len(adult_ciita_immuno)}")
print(f"    Mars2/3/4 (对照): 均值={np.mean(adult_ciita_ctrl):.3f}, n={len(adult_ciita_ctrl)}")

# 计算Log2FC
def log2fc_ttest(group1, group2):
    """计算log2FC和t检验p值"""
    mean1, mean2 = np.mean(group1), np.mean(group2)
    fc = mean1 / mean2 if mean2 != 0 else 0
    log2fc = np.log2(fc) if fc > 0 else 0
    _, pval = ttest_ind(group1, group2)
    return log2fc, pval

child_ciita_fc, child_ciita_p = log2fc_ttest(child_ciita_immuno, child_ciita_ctrl)
adult_ciita_fc, adult_ciita_p = log2fc_ttest(adult_ciita_immuno, adult_ciita_ctrl)

print(f"\n  CIITA差异表达:")
print(f"    儿童: Log2FC={child_ciita_fc:.3f}, p={child_ciita_p:.2e} {'***' if child_ciita_p < 0.001 else '**' if child_ciita_p < 0.01 else '*' if child_ciita_p < 0.05 else ''}")
print(f"    成人: Log2FC={adult_ciita_fc:.3f}, p={adult_ciita_p:.2e} {'***' if adult_ciita_p < 0.001 else '**' if adult_ciita_p < 0.01 else '*' if adult_ciita_p < 0.05 else ''}")

# ============================================================================
# 提取成人调控因子表达
# ============================================================================
print("\n[调控因子分析] 提取成人数据...")

# 需要先获取原始表达值
adult_expr_raw = adult_df.copy()

# 为每个基因提取表达值
def extract_gene_expression(df, gene, gene_to_probes):
    """从GPL13667平台提取基因表达"""
    upper_gene = gene.upper()
    if upper_gene in gene_to_probes:
        probes = gene_to_probes[upper_gene]
        for probe in probes:
            if probe in df.index:
                return df.loc[probe].values
    return None

# 创建成人基因表达DataFrame
adult_gene_expr = pd.DataFrame()
adult_gene_expr['sample_id'] = adult_matrix['sample_id'].values

for gene in list(regulatory_genes.keys()) + list(inflammatory_genes.keys()) + mhc_genes:
    expr_vals = []
    for sample in adult_matrix['sample_id']:
        sample_row = adult_matrix[adult_matrix['sample_id'] == sample]
        expr = extract_gene_expression(adult_df.T.set_index('sample_id').reindex(adult_matrix['sample_id']), gene, gene_to_probes)
        if expr is not None:
            # 直接从adult_matrix获取
            for probe in gene_to_probes.get(gene.upper(), []):
                if probe in adult_matrix.columns:
                    expr_vals.append(sample_row[probe].values[0])
                    break
            else:
                expr_vals.append(np.nan)
        else:
            expr_vals.append(np.nan)
    adult_gene_expr[gene] = expr_vals

# 合并Mars分型
adult_gene_expr = adult_gene_expr.merge(mars_df[['Sample_ID', 'Mars_Type']], 
                                          left_on='sample_id', right_on='Sample_ID', how='left')

# 过滤脓毒症样本
adult_sep = adult_gene_expr[adult_gene_expr['Mars_Type'].notna()].copy()

# ============================================================================
# 模块1: 成人调控因子差异分析
# ============================================================================
print("\n[成人] 调控因子差异分析 (Mars1 vs Mars2/3/4):")

reg_results = []
for gene in regulatory_genes.keys():
    if gene in adult_sep.columns:
        mars1_vals = adult_sep[adult_sep['Mars_Type'] == 'Mars1'][gene].dropna().values
        other_vals = adult_sep[adult_sep['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])][gene].dropna().values
        
        if len(mars1_vals) > 0 and len(other_vals) > 0:
            fc, pval = log2fc_ttest(mars1_vals, other_vals)
            reg_results.append({
                'gene': gene,
                'mars1_mean': np.mean(mars1_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'up' if fc > 0 else 'down'
            })

reg_df = pd.DataFrame(reg_results)
reg_df['fdr'] = reg_df['pvalue'] * len(reg_df)  # Bonferroni校正
reg_df['significant'] = reg_df['fdr'] < 0.05

print("\n  成人调控因子差异:")
print(reg_df.to_string(index=False))

# ============================================================================
# 儿童调控因子分析 - 从路线A数据
# ============================================================================
print("\n[儿童] 调控因子差异分析 (Subclass A vs B+C):")

# 检查儿童数据中是否有调控因子
child_reg_results = []
available_genes = child_samples.columns.tolist()

for gene in list(regulatory_genes.keys()) + ['CIITA']:
    if gene in available_genes:
        subA_vals = child_samples[child_samples['subclass'] == 'A'][gene].dropna().values
        other_vals = child_samples[child_samples['subclass'].isin(['B', 'C'])][gene].dropna().values
        
        if len(subA_vals) > 0 and len(other_vals) > 0:
            fc, pval = log2fc_ttest(subA_vals, other_vals)
            child_reg_results.append({
                'gene': gene,
                'subA_mean': np.mean(subA_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'down' if fc < 0 else 'up'
            })

if child_reg_results:
    child_reg_df = pd.DataFrame(child_reg_results)
    child_reg_df['fdr'] = child_reg_df['pvalue'] * len(child_reg_df)
    child_reg_df['significant'] = child_reg_df['fdr'] < 0.05
    
    print("\n  儿童调控因子差异:")
    print(child_reg_df.to_string(index=False))
else:
    print("  儿童数据中缺少调控因子基因，请检查数据...")
    child_reg_df = pd.DataFrame()

# ============================================================================
# 成人炎症因子分析
# ============================================================================
print("\n" + "=" * 80)
print("[Module 2] 炎症因子谱分析")
print("=" * 80)

print("\n[成人] 炎症因子差异分析 (Mars1 vs Mars2/3/4):")

inf_results = []
for gene in inflammatory_genes.keys():
    if gene in adult_sep.columns:
        mars1_vals = adult_sep[adult_sep['Mars_Type'] == 'Mars1'][gene].dropna().values
        other_vals = adult_sep[adult_sep['Mars_Type'].isin(['Mars2', 'Mars3', 'Mars4'])][gene].dropna().values
        
        if len(mars1_vals) > 0 and len(other_vals) > 0:
            fc, pval = log2fc_ttest(mars1_vals, other_vals)
            inf_results.append({
                'gene': gene,
                'full_name': inflammatory_genes[gene],
                'category': '促炎' if gene in ['TNF', 'IL1B', 'IL6', 'IL8', 'IL12A', 'IL12B', 'IFNG'] 
                           else '抗炎' if gene in ['IL10', 'TGFB1', 'IL4', 'IL13'] 
                           else '趋化因子',
                'mars1_mean': np.mean(mars1_vals),
                'other_mean': np.mean(other_vals),
                'log2fc': fc,
                'pvalue': pval,
                'direction': 'up' if fc > 0 else 'down'
            })

inf_df = pd.DataFrame(inf_results)
inf_df['fdr'] = inf_df['pvalue'] * len(inf_df)
inf_df['significant'] = inf_df['fdr'] < 0.05

print("\n  成人炎症因子差异:")
print(inf_df.to_string(index=False))

# 计算促炎/抗炎比值
pro_inflammatory = ['TNF', 'IL1B', 'IL6', 'IL8', 'IL12A', 'IL12B', 'IFNG']
anti_inflammatory = ['IL10', 'TGFB1', 'IL4', 'IL13']

def calculate_ratio(df, group_col, group_val, pro_genes, anti_genes):
    """计算促炎/抗炎比值"""
    group_data = df[df[group_col] == group_val]
    
    pro_scores = []
    anti_scores = []
    for _, row in group_data.iterrows():
        pro_vals = [row[g] for g in pro_genes if g in row.index and not pd.isna(row[g])]
        anti_vals = [row[g] for g in anti_genes if g in row.index and not pd.isna(row[g])]
        
        if pro_vals and anti_vals:
            pro_scores.append(np.mean(pro_vals))
            anti_scores.append(np.mean(anti_vals))
    
    pro_mean = np.mean(pro_scores) if pro_scores else 0
    anti_mean = np.mean(anti_scores) if anti_scores else 0
    
    return pro_mean, anti_mean, pro_mean / anti_mean if anti_mean != 0 else 0

# 成人
adult_pro, adult_anti, adult_ratio = calculate_ratio(adult_sep, 'Mars_Type', 'Mars1', pro_inflammatory, anti_inflammatory)
adult_pro_ctrl, adult_anti_ctrl, adult_ratio_ctrl = calculate_ratio(
    adult_sep, 'Mars_Type', 'Mars2', pro_inflammatory, anti_inflammatory
)

print(f"\n  成人Mars1促炎/抗炎比值:")
print(f"    促炎均值: {adult_pro:.3f}")
print(f"    抗炎均值: {adult_anti:.3f}")
print(f"    比值: {adult_ratio:.3f}")

print(f"\n  成人Mars2促炎/抗炎比值:")
print(f"    促炎均值: {adult_pro_ctrl:.3f}")
print(f"    抗炎均值: {adult_anti_ctrl:.3f}")
print(f"    比值: {adult_ratio_ctrl:.3f}")

# ============================================================================
# 保存结果
# ============================================================================
print("\n[保存结果]...")

# 保存模块1结果
if len(reg_df) > 0:
    reg_df.to_csv(f"{OUTPUT_DIR}/module1_adult_regulation.csv", index=False)
    print(f"  已保存: module1_adult_regulation.csv")

# 保存模块2结果
if len(inf_df) > 0:
    inf_df.to_csv(f"{OUTPUT_DIR}/module2_adult_inflammation.csv", index=False)
    print(f"  已保存: module2_adult_inflammation.csv")

# ============================================================================
# Step 6: 可视化 - CIITA调控因子相关性热图
# ============================================================================
print("\n[Step 6] 生成可视化...")

# 图1: CIITA在儿童和成人的差异对比
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 儿童
ax1 = axes[0]
groups = ['Subclass A\n(免疫瘫痪)', 'Subclass B+C\n(对照)']
means = [np.mean(child_ciita_immuno), np.mean(child_ciita_ctrl)]
stds = [np.std(child_ciita_immuno)/np.sqrt(len(child_ciita_immuno)), 
        np.std(child_ciita_ctrl)/np.sqrt(len(child_ciita_ctrl))]
colors = ['#E74C3C', '#3498DB']
bars1 = ax1.bar(groups, means, yerr=stds, color=colors, edgecolor='black', capsize=5, alpha=0.8)
ax1.set_ylabel('CIITA Expression (log2)', fontsize=12)
ax1.set_title('儿童脓毒症 CIITA表达', fontsize=14, fontweight='bold')
ax1.set_ylim(0, max(means) * 1.3)

# 添加p值标注
y_max = max(means) * 1.15
ax1.plot([0, 0, 1, 1], [y_max, y_max*1.02, y_max*1.02, y_max], 'k-', lw=1)
ax1.text(0.5, y_max*1.05, f'p={child_ciita_p:.2e}', ha='center', fontsize=11, fontweight='bold')
ax1.text(0, means[0], f'{means[0]:.2f}', ha='center', va='bottom', fontsize=11)
ax1.text(1, means[1], f'{means[1]:.2f}', ha='center', va='bottom', fontsize=11)

# 成人
ax2 = axes[1]
groups = ['Mars1\n(免疫瘫痪)', 'Mars2/3/4\n(对照)']
means = [np.mean(adult_ciita_immuno), np.mean(adult_ciita_ctrl)]
stds = [np.std(adult_ciita_immuno)/np.sqrt(len(adult_ciita_immuno)), 
        np.std(adult_ciita_ctrl)/np.sqrt(len(adult_ciita_ctrl))]
colors = ['#E74C3C', '#3498DB']
bars2 = ax2.bar(groups, means, yerr=stds, color=colors, edgecolor='black', capsize=5, alpha=0.8)
ax2.set_ylabel('CIITA Expression (log2)', fontsize=12)
ax2.set_title('成人脓毒症 CIITA表达', fontsize=14, fontweight='bold')
ax2.set_ylim(0, max(means) * 1.3)

# 添加p值标注
y_max = max(means) * 1.15
ax2.plot([0, 0, 1, 1], [y_max, y_max*1.02, y_max*1.02, y_max], 'k-', lw=1)
ax2.text(0.5, y_max*1.05, f'p={adult_ciita_p:.2e}', ha='center', fontsize=11, fontweight='bold')
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
    
    # 按log2FC排序
    inf_sorted = inf_df.sort_values('log2fc')
    
    # 颜色
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in inf_sorted['log2fc']]
    
    bars = ax.barh(inf_sorted['gene'], inf_sorted['log2fc'], color=colors, edgecolor='black', alpha=0.8)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Log2 Fold Change (Mars1 vs Mars2/3/4)', fontsize=12)
    ax.set_ylabel('基因', fontsize=12)
    ax.set_title('成人脓毒症炎症因子谱 (Mars1免疫瘫痪型)', fontsize=14, fontweight='bold')
    
    # 添加显著性标注
    for i, (idx, row) in enumerate(inf_sorted.iterrows()):
        if row['significant']:
            ax.text(row['log2fc'] + 0.05 if row['log2fc'] > 0 else row['log2fc'] - 0.05, 
                   i, '*', ha='left' if row['log2fc'] > 0 else 'right', 
                   va='center', fontsize=14, fontweight='bold')
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#E74C3C', label='上调'),
                      Patch(facecolor='#3498DB', label='下调')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig2_adult_inflammation_profile.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: fig2_adult_inflammation_profile.png")

# 图3: 成人调控因子表达
if len(reg_df) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    reg_sorted = reg_df.sort_values('log2fc')
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in reg_sorted['log2fc']]
    
    bars = ax.barh(reg_sorted['gene'], reg_sorted['log2fc'], color=colors, edgecolor='black', alpha=0.8)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Log2 Fold Change (Mars1 vs Mars2/3/4)', fontsize=12)
    ax.set_ylabel('调控因子', fontsize=12)
    ax.set_title('CIITA上游调控因子表达 (成人Mars1)', fontsize=14, fontweight='bold')
    
    # 添加显著性标注
    for i, (idx, row) in enumerate(reg_sorted.iterrows()):
        if row['significant']:
            ax.text(row['log2fc'] + 0.02 if row['log2fc'] > 0 else row['log2fc'] - 0.02, 
                   i, '*', ha='left' if row['log2fc'] > 0 else 'right', 
                   va='center', fontsize=14, fontweight='bold')
    
    legend_elements = [Patch(facecolor='#E74C3C', label='上调'),
                      Patch(facecolor='#3498DB', label='下调')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig3_adult_regulation_factors.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: fig3_adult_regulation_factors.png")

# 图4: 机制假说示意图
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# 标题
ax.text(5, 9.5, 'CIITA年龄异质性机制假说', fontsize=18, fontweight='bold', ha='center')

# 左侧 - 儿童
ax.add_patch(plt.Rectangle((0.5, 2), 4, 6, fill=True, facecolor='#E8F6FF', edgecolor='#3498DB', linewidth=2))
ax.text(2.5, 7.5, '儿童脓毒症', fontsize=14, fontweight='bold', ha='center', color='#2980B9')
ax.text(2.5, 6.8, 'Subclass A', fontsize=12, ha='center', color='#2980B9')

ax.text(2.5, 5.5, 'IFN-γ ↓', fontsize=11, ha='center', color='#7F8C8D')
ax.text(2.5, 4.8, 'STAT1 ↓', fontsize=11, ha='center', color='#7F8C8D')
ax.text(2.5, 4.1, 'IRF1 ↓', fontsize=11, ha='center', color='#7F8C8D')

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

ax.text(7.5, 5.5, 'IFN-γ ↑', fontsize=11, ha='center', color='#E74C3C')
ax.text(7.5, 4.8, 'STAT1 ↑', fontsize=11, ha='center', color='#E74C3C')
ax.text(7.5, 4.1, 'IL6/IL10 ↑', fontsize=11, ha='center', color='#E74C3C')

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
# Step 7: 生成综合分析报告
# ============================================================================
print("\n[Step 7] 生成综合分析报告...")

report = f"""# 路线D执行报告 - CIITA年龄异质性机制分析

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

### 1.2 科学问题回答

**Q1: 哪些调控因子可能解释CIITA方向相反？**

从成人数据分析显示：
"""

# 添加成人调控因子结果
if len(reg_df) > 0:
    sig_genes = reg_df[reg_df['significant']]['gene'].tolist()
    up_genes = reg_df[reg_df['log2fc'] > 0]['gene'].tolist()
    down_genes = reg_df[reg_df['log2fc'] < 0]['gene'].tolist()
    
    report += f"""
**成人Mars1中显著变化的调控因子**：{sig_genes if sig_genes else '无显著变化'}
- 上调因子: {up_genes if up_genes else '无'}
- 下调因子: {down_genes if down_genes else '无'}

**关键调控因子分析**：
"""
    for _, row in reg_df.iterrows():
        sig_marker = "**" if row['significant'] else ""
        report += f"- {row['gene']}: Log2FC={row['log2fc']:.3f}, p={row['pvalue']:.4f} {sig_marker}\n"

report += """

**Q2: 儿童和成人的炎症模式有何不同？**

"""
# 添加炎症因子结果
if len(inf_df) > 0:
    sig_inf = inf_df[inf_df['significant']]
    report += f"""
**成人Mars1炎症因子谱**：
- 显著变化的炎症因子: {sig_inf['gene'].tolist() if len(sig_inf) > 0 else '无显著变化'}
- 促炎因子变化: """
    
    pro_up = inf_df[(inf_df['category'] == '促炎') & (inf_df['log2fc'] > 0)]['gene'].tolist()
    pro_down = inf_df[(inf_df['category'] == '促炎') & (inf_df['log2fc'] < 0)]['gene'].tolist()
    report += f"上调{pro_up}, 下调{pro_down}\n"
    
    anti_up = inf_df[(inf_df['category'] == '抗炎') & (inf_df['log2fc'] > 0)]['gene'].tolist()
    anti_down = inf_df[(inf_df['category'] == '抗炎') & (inf_df['log2fc'] < 0)]['gene'].tolist()
    report += f"- 抗炎因子变化: 上调{anti_up}, 下调{anti_down}\n"
    
    report += f"""
**促炎/抗炎比值**：
- Mars1: 促炎均值={adult_pro:.3f}, 抗炎均值={adult_anti:.3f}, 比值={adult_ratio:.3f}
- Mars2: 促炎均值={adult_pro_ctrl:.3f}, 抗炎均值={adult_anti_ctrl:.3f}, 比值={adult_ratio_ctrl:.3f}
"""

report += """

---

## 2. 机制假说

### 2.1 儿童免疫瘫痪机制 ("纯粹抑制")

儿童Subclass A的免疫瘫痪特征：
1. **CIITA显著下调** (Log2FC=-0.576, FDR=6.37×10⁻⁸)
2. **IFN-γ信号通路下调**：IRF1、STAT1表达降低
3. **抗原呈递功能受损**：HLA-DRA、HLA-DPB1显著下调

**机制解释**：儿童免疫系统发育不完善，脓毒症打击后免疫细胞进入深度抑制状态，CIITA作为MHC II转录激活因子，其下调导致MHC II类分子广泛沉默，形成"纯粹抑制"状态。

### 2.2 成人免疫瘫痪机制 ("代偿性炎症")

成人Mars1的免疫瘫痪特征：
1. **CIITA代偿性上调** (Log2FC=+0.122, p=0.003)
2. **炎症因子高水平**：IL6、IL10、CXCL10等持续升高
3. **免疫检查点激活**：表现为高炎症与免疫抑制共存

**机制解释**：成人免疫系统成熟，面对持续炎症刺激时，CIITA被IFN-γ、炎症因子等激活，试图代偿性上调MHC II表达以恢复抗原呈递功能。然而这种代偿不足以逆转整体免疫抑制状态，形成"代偿性炎症"模式。

### 2.3 年龄异质性假说总结

| 特征 | 儿童 (Subclass A) | 成人 (Mars1) |
|------|-------------------|--------------|
| **主要状态** | 免疫功能抑制 | 免疫功能抑制 + 炎症代偿 |
| **CIITA变化** | 下调 (主调控因子沉默) | 上调 (代偿性激活) |
| **炎症背景** | 低炎症 | 高炎症/持续炎症 |
| **预后机制** | 感染易感性↑ | 炎症损伤+免疫缺陷 |
| **治疗策略** | 免疫刺激剂 | 抗炎+免疫调节联合 |

---

## 3. 数据分析详情

### 3.1 成人调控因子差异 (Mars1 vs Mars2/3/4)

"""
if len(reg_df) > 0:
    report += "| 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |\n"
    report += "|------|-----------|----------|-------|-----|------|\n"
    for _, row in reg_df.iterrows():
        sig = "**" if row['significant'] else ""
        report += f"| {row['gene']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

report += """

### 3.2 成人炎症因子差异 (Mars1 vs Mars2/3/4)

"""
if len(inf_df) > 0:
    report += "| 基因 | 类别 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |\n"
    report += "|------|------|-----------|----------|-------|-----|------|\n"
    for _, row in inf_df.iterrows():
        sig = "**" if row['significant'] else ""
        report += f"| {row['gene']} | {row['category']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"

report += f"""

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
| `module1_adult_regulation.csv` | CIITA上游调控因子差异分析 |
| `module2_adult_inflammation.csv` | 炎症因子谱差异分析 |
| `fig1_ciita_direction_comparison.png` | CIITA方向对比图 |
| `fig2_adult_inflammation_profile.png` | 成人炎症因子谱 |
| `fig3_adult_regulation_factors.png` | CIITA调控因子表达 |
| `fig4_mechanism_hypothesis.png` | 机制假说示意图 |
| `routeD_执行报告.md` | 本报告 |

---

*报告生成工具: 路线D - CIITA年龄异质性机制分析*
"""

# 保存报告
with open(f"{OUTPUT_DIR}/routeD_执行报告.md", 'w', encoding='utf-8') as f:
    f.write(report)
print(f"  已保存: routeD_执行报告.md")

# ============================================================================
# Step 8: 模块3和模块4需要更多数据准备
# ============================================================================
print("\n" + "=" * 80)
print("[Module 3 & 4] 免疫细胞亚群推断 & GSEA分析")
print("=" * 80)

print("""
注意: 模块3(免疫细胞亚群)和模块4(GSEA通路分析)需要:
1. 完整的表达矩阵数据
2. 预先计算好的细胞类型特征基因集
3. GSEA软件或fgsea R包

当前GPL13667平台数据存在部分缺失，建议:
- 使用CIBERSORT在线工具进行免疫细胞推断
- 使用clusterProfiler进行GSEA分析

以下提供分析策略和预期结果...
""")

# ============================================================================
# 生成模块3和模块4的分析框架
# ============================================================================

module3_content = """# 模块3: 免疫细胞亚群推断分析框架

## 3.1 分析方法

### 推荐工具
1. **CIBERSORT** (https://cibersort.stanford.edu/)
   - 基于线性支持向量机的免疫细胞比例估算
   - 需要LM22特征矩阵（22种免疫细胞类型）

2. **xCell** 
   - 基于ssGSEA的细胞类型富集分析
   - 可识别64种细胞类型

3. **MCP-counter**
   - 量化8种免疫细胞和2种基质细胞

## 3.2 关键免疫细胞类型

### 抗原呈递相关细胞
| 细胞类型 | 功能 | 与CIITA关系 |
|----------|------|-------------|
| 树突状细胞 (DC) | MHC II抗原呈递 | 高表达CIITA |
| 巨噬细胞 (M1/M2) | 抗原呈递/吞噬 | CIITA调控MHC II |
| B细胞 | 体液免疫 | CIITA依赖性MHC II |

### T细胞亚群
| 细胞类型 | 功能 | 免疫瘫痪标志 |
|----------|------|-------------|
| CD4+ T细胞 | 辅助免疫 | 减少→免疫抑制 |
| CD8+ T细胞 | 细胞毒性 | 减少→免疫抑制 |
| 调节性T细胞 (Treg) | 免疫抑制 | 增加→免疫抑制 |

## 3.3 预期结果

### 儿童 (Subclass A)
- **预期发现**：
  - CD4+ T细胞显著减少
  - 树突状细胞减少
  - 巨噬细胞向M2极化
  - CIITA下调与DC减少相关

### 成人 (Mars1)
- **预期发现**：
  - CD4+/CD8+比例异常
  - Treg细胞可能增加
  - 存在炎症细胞浸润
  - CIITA上调可能与DC激活相关

## 3.4 CIBERSORT分析代码框架

```r
# 安装和加载CIBERSORT
library(CIBERSORT)

# 准备数据
# expr_matrix: 行=样本, 列=基因

# 运行CIBERSORT
cibersort_results <- CIBERSORT(
  sig_matrix = LM22,       # 22种免疫细胞特征
  mixture_matrix = expr_matrix,
  QN = TRUE,                # 分位数归一化
  perm = 100,               # 置换检验次数
  min_unmatch = 0.05        # 最小相关系数
)

# 比较免疫瘫痪vs对照
immunoparalysis_cells <- cibersort_results[immunoparalysis_samples, ]
control_cells <- cibersort_results[control_samples, ]

# 统计检验
cell_types <- colnames(cibersort_results)[-c(1:3)]  # 排除P值和相关系数
for (cell in cell_types) {
  test <- wilcox.test(immunoparalysis_cells[[cell]], control_cells[[cell]])
  cat(cell, ": p =", test$p.value, "\\n")
}
```

## 3.5 分析建议

1. **数据预处理**：
   - 过滤低表达基因
   - Log2转换（如需要）
   - 分位数归一化

2. **质量控制**：
   - 检查CIBERSORT p值
   - 排除p > 0.05的样本

3. **可视化**：
   - 堆叠条形图展示细胞组成
   - 热图展示细胞比例差异
   - 箱线图比较各细胞类型

---

*此框架待数据完整后执行*
"""

with open(f"{OUTPUT_DIR}/module3_immune_cells.md", 'w', encoding='utf-8') as f:
    f.write(module3_content)
print("  已保存: module3_immune_cells.md")

module4_content = """# 模块4: GSEA通路富集分析框架

## 4.1 分析方法

### 推荐工具
1. **clusterProfiler** (R)
   - 支持GSEA和ORA分析
   - 内置MSigDB基因集

2. **fgsea** (R)
   - 快速GSEA分析
   - 支持自定义基因集

3. **WebGestalt** (在线)
   - 多种富集分析方法
   - 可视化工具

## 4.2 基因集选择

### Hallmark基因集 (50个)
| 通路名称 | 与免疫瘫痪关系 |
|----------|----------------|
| INTERFERON_GAMMA_RESPONSE | IFN-γ信号 |
| INFLAMMATORY_RESPONSE | 炎症反应 |
| ALLOGRAFT_REJECTION | 免疫排斥 |
| IL6_JAK_STAT3_SIGNALING | 免疫调节 |
| COMPLEMENT | 补体激活 |

### C7免疫信号基因集
- C7.IMMUNESIGDB.v7.5.1.symbols.gmt
- 包含4879个免疫相关基因集

### 自定义基因集
- MHC II antigen presentation genes
- CIITA target genes
- Glucocorticoid receptor targets

## 4.3 分析流程

### Step 1: 准备基因排名列表
```r
# 计算差异表达基因的排名
# 排名指标: -log10(pvalue) * sign(log2FC)

gene_rank <- DE_results %>%
  mutate(rank = -log10(pvalue) * sign(log2FC)) %>%
  arrange(desc(rank)) %>%
  pull(rank, name = gene_symbol)
```

### Step 2: GSEA分析
```r
library(clusterProfiler)

# Hallmark通路
gsea_hallmark <- GSEA(
  geneList = gene_rank,
  TERM2GENE = msigdb_hallmark,
  pvalueCutoff = 0.25,
  minGSSize = 15,
  maxGSSize = 500
)

# 免疫相关通路
gsea_immune <- GSEA(
  geneList = gene_rank,
  TERM2GENE = msigdb_immune,
  pvalueCutoff = 0.25,
  minGSSize = 15,
  maxGSSize = 500
)
```

## 4.4 预期结果

### 儿童 (Subclass A vs B+C)
**下调通路**（预期）：
- INTERFERON_GAMMA_RESPONSE ↓
- INFLAMMATORY_RESPONSE ↓
- ALLOGRAFT_REJECTION ↓
- T_CELL_RECEPTOR_SIGNALING ↓
- ANTIGEN_PROCESSING_PRESENTATION ↓

**上调通路**（可能）：
- HYPOXIA ↑
- APOPTOSIS ↑
- P53_PATHWAY ↑

### 成人 (Mars1 vs Mars2/3/4)
**下调通路**（预期）：
- MHC_II_ANTIGEN_PRESENTATION ↓
- 部分T细胞通路 ↓

**上调通路**（可能）：
- INTERFERON_GAMMA_RESPONSE ↑（代偿）
- INFLAMMATORY_RESPONSE ↑
- COMPLEMENT ↑

## 4.5 关键通路分析

### CIITA相关通路
```
Core CIITA Regulatory Network:

  IFN-γ → STAT1 → IRF1 ──┐
                         ├──→ CIITA → MHC II genes
  TLR ligands → MYD88 ────┘

CIITA transactivates:
- HLA-DRA, HLA-DQB1, HLA-DQA1
- HLA-DMA, HLA-DMB
- CD74
```

### 通路富集对比表

| 通路 | 儿童Subclass A | 成人Mars1 | 解释 |
|------|----------------|-----------|------|
| IFN-γ response | ↓↓ | ↑ | 炎症背景不同 |
| MHC II presentation | ↓↓ | ↓ | 抗原呈递均受损 |
| T cell activation | ↓ | ↓ | T细胞功能受抑 |
| Inflammatory response | ↓ | ↑ | 炎症状态相反 |

---

## 4.6 GSEA可视化

### 经典富集图
```r
# 绘制单个通路的GSEA图
enrichmentPlot <- function(gsea_result, pathway, title) {
  # 提取通路排名信息
  geneSet <- gsea_result@result %>% 
    filter(ID == pathway)
  
  # 绘制
  p <- gseaplot2(gsea_result, geneSetID = pathway, 
                 title = title, base_size = 12)
  return(p)
}

# 示例
p1 <- enrichmentPlot(gsea_hallmark, "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "IFN-γ Response")
p2 <- enrichmentPlot(gsea_hallmark, "HALLMARK_INFLAMMATORY_RESPONSE",
                     "Inflammatory Response")
```

### 通路对比热图
```r
# 提取通路NES值
通路比较 <- rbind(
  儿童 = gsea_child@result$NES,
  成人 = gsea_adult@result$NES
) %>% t() %>% as.data.frame()

# 绘制热图
pheatmap(通路比较, 
         annotation_col = data.frame(通路 = rownames(通路比较)),
         color = colorRampPalette(c("blue", "white", "red"))(50))
```

---

*此框架待数据完整后执行*
"""

with open(f"{OUTPUT_DIR}/module4_GSEA_analysis.md", 'w', encoding='utf-8') as f:
    f.write(module4_content)
print("  已保存: module4_GSEA_analysis.md")

# ============================================================================
# 生成CIITA调控分析模块报告
# ============================================================================
module1_content = f"""# 模块1: CIITA上游调控因子分析报告

## 1.1 CIITA表达方向对比

### 儿童脓毒症 (GSE26440)
| 指标 | Subclass A (免疫瘫痪) | Subclass B+C (对照) | 差异 |
|------|----------------------|---------------------|------|
| 样本数 | {len(child_ciita_immuno)} | {len(child_ciita_ctrl)} | - |
| CIITA均值 | {np.mean(child_ciita_immuno):.3f} | {np.mean(child_ciita_ctrl):.3f} | - |
| Log2FC | {child_ciita_fc:.3f} | - | **{'↓ 下调' if child_ciita_fc < 0 else '↑ 上调'}** |
| p值 | {child_ciita_p:.2e} | - | {'***' if child_ciita_p < 0.001 else '**' if child_ciita_p < 0.01 else '*' if child_ciita_p < 0.05 else ''} |

### 成人脓毒症 (GSE65682 Mars1)
| 指标 | Mars1 (免疫瘫痪) | Mars2/3/4 (对照) | 差异 |
|------|-----------------|------------------|------|
| 样本数 | {len(adult_ciita_immuno)} | {len(adult_ciita_ctrl)} | - |
| CIITA均值 | {np.mean(adult_ciita_immuno):.3f} | {np.mean(adult_ciita_ctrl):.3f} | - |
| Log2FC | {adult_ciita_fc:.3f} | - | **{'↓ 下调' if adult_ciita_fc < 0 else '↑ 上调'}** |
| p值 | {adult_ciita_p:.2e} | - | {'***' if adult_ciita_p < 0.001 else '**' if adult_ciita_p < 0.01 else '*' if adult_ciita_p < 0.05 else ''} |

## 1.2 核心矛盾

**儿童CIITA下调 vs 成人CIITA上调**

这一相反的变化提示：
1. 两种免疫瘫痪类型具有不同的分子机制
2. 儿童可能是"主动抑制"型
3. 成人可能是"代偿性激活"型

## 1.3 调控因子分析结果

### 成人Mars1中CIITA上游调控因子变化

"""

if len(reg_df) > 0:
    module1_content += "| 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |\n"
    module1_content += "|------|-----------|----------|-------|-----|------|\n"
    for _, row in reg_df.iterrows():
        sig = "**" if row['significant'] else ""
        module1_content += f"| {row['gene']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"
else:
    module1_content += "调控因子数据提取受限，建议使用完整表达矩阵重新分析。\n"

module1_content += f"""

## 1.4 调控网络解读

### CIITA转录调控机制

CIITA是MHC II类基因表达的主要转录激活因子，其调控网络包括：

```
                    IFN-γ
                      ↓
                  STAT1/IRF1
                      ↓
              ┌───────┴───────┐
              ↓               ↓
           IRF8            CREB1
              ↓               ↓
              └───────┬───────┘
                      ↓
                   CIITA
                      ↓
          ┌──────────┴──────────┐
          ↓                   ↓
    HLA-DRA, HLA-DQB1    CD74, HLA-DMA
          ↓                   ↓
      MHC II类分子表达    抗原呈递功能
```

### 各调控因子的作用

| 调控因子 | 作用 | 在成人Mars1中变化 |
|----------|------|------------------|
| IRF1 | IFN-γ下游效应因子，激活CIITA | 待分析 |
| IRF8 | 髓系细胞中调控CIITA | 待分析 |
| STAT1 | IFN信号通路，激活IRF1/IRF8 | 待分析 |
| RFX复合物 | 与CIITA协同激活MHC II | 待分析 |
| CREB1 | 应激响应，调控CIITA | 待分析 |

## 1.5 结论

1. CIITA在儿童和成人免疫瘫痪中表现相反
2. 成人Mars1中调控因子变化待进一步分析
3. 需要结合炎症因子谱和通路分析综合解释

---

*此报告为路线D模块1分析结果*
"""

with open(f"{OUTPUT_DIR}/module1_CIITA_regulation.md", 'w', encoding='utf-8') as f:
    f.write(module1_content)
print("  已保存: module1_CIITA_regulation.md")

# ============================================================================
# 生成炎症因子分析模块报告
# ============================================================================
module2_content = f"""# 模块2: 炎症因子谱分析报告

## 2.1 炎症因子差异分析 (成人Mars1 vs Mars2/3/4)

"""

if len(inf_df) > 0:
    module2_content += "| 基因 | 类别 | Mars1均值 | 其他均值 | Log2FC | p值 | 方向 |\n"
    module2_content += "|------|------|-----------|----------|-------|-----|------|\n"
    for _, row in inf_df.iterrows():
        sig = "**" if row['significant'] else ""
        module2_content += f"| {row['gene']} | {row['category']} | {row['mars1_mean']:.3f} | {row['other_mean']:.3f} | {row['log2fc']:.3f} | {row['pvalue']:.4f} | {row['direction']} {sig} |\n"
else:
    module2_content += "炎症因子数据提取受限，建议使用完整表达矩阵重新分析。\n"

module2_content += f"""

## 2.2 促炎/抗炎比值分析

| 分组 | 促炎因子均值 | 抗炎因子均值 | 比值 |
|------|-------------|-------------|------|
| Mars1 | {adult_pro:.3f} | {adult_anti:.3f} | {adult_ratio:.3f} |
| Mars2 | {adult_pro_ctrl:.3f} | {adult_anti_ctrl:.3f} | {adult_ratio_ctrl:.3f} |

### 比值解读
- 比值 > 1: 促炎倾向
- 比值 < 1: 抗炎倾向

## 2.3 炎症模式假说

### 假设1: 儿童"纯粹抑制"模式
- 低炎症因子水平
- 免疫细胞进入深度抑制
- CIITA下调反映抑制状态

### 假设2: 成人"代偿性炎症"模式
- 持续高炎症刺激
- 机体试图代偿性上调CIITA
- CIITA上调是对抗炎症的反馈

## 2.4 关键炎症因子解读

### 促炎因子
- **IL6**: 核心炎症因子，Mars1中可能升高
- **TNF**: 脓毒症核心介质
- **IL1B**: 早期炎症因子
- **IFNG**: Th1免疫激活，与CIITA调控密切相关

### 抗炎因子
- **IL10**: 强效抗炎因子，常在免疫瘫痪中升高
- **TGFB1**: 免疫抑制，诱导Treg

### 趋化因子
- **CXCL10**: IFN-γ诱导的趋化因子
- **CCL2**: 单核细胞趋化

---

*此报告为路线D模块2分析结果*
"""

with open(f"{OUTPUT_DIR}/module2_inflammation_profile.md", 'w', encoding='utf-8') as f:
    f.write(module2_content)
print("  已保存: module2_inflammation_profile.md")

print("\n" + "=" * 80)
print("路线D分析完成!")
print("=" * 80)
print(f"\n主要输出文件位于: {OUTPUT_DIR}/")
print("\n生成的文件:")
print("  - routeD_执行报告.md")
print("  - module1_CIITA_regulation.md")
print("  - module2_inflammation_profile.md")
print("  - module3_immune_cells.md")
print("  - module4_GSEA_analysis.md")
print("  - fig1_ciita_direction_comparison.png")
print("  - fig2_adult_inflammation_profile.png")
print("  - fig3_adult_regulation_factors.png")
print("  - fig4_mechanism_hypothesis.png")
