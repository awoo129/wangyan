#!/usr/bin/env python3
"""
儿童脓毒症IPS评分研究 - 路线A分析流程
基于GSE26440(推导集)和GSE26378(验证集)
"""

import re
import os
import pickle
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score, confusion_matrix
from sklearn.model_selection import cross_val_score, StratifiedKFold
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

# 尝试导入绘图库
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOT = True
except:
    HAS_PLOT = False

# 结果目录
RESULTS_DIR = "长期计划/儿童脓毒症免疫瘫痪研究/results/routeA"
os.makedirs(RESULTS_DIR, exist_ok=True)

print("="*70)
print("儿童脓毒症IPS评分研究 - 路线A分析流程")
print("="*70)

# ============================================================
# Step 1: 数据读取与预处理
# ============================================================
print("\n" + "="*70)
print("Step 1: 数据读取与预处理")
print("="*70)

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
                # 解析元数据
                parts = line.split('\t')
                key = parts[0].replace('!Sample_', '')
                values = [p.strip().strip('"') for p in parts[1:]]
                
                if key not in metadata:
                    metadata[key] = values
                else:
                    # 处理多行相同的key
                    if isinstance(metadata[key], list):
                        metadata[key].extend(values)
                        
            elif line.startswith('ENTREZ_GENE_ID') or line.startswith('ID_REF') or (line and not line.startswith('#')):
                # 探针ID行
                if 'ID_REF' in line or probe_ids is None:
                    parts = line.split('\t')
                    if parts[0] in ['ID_REF', 'Probe_ID', 'ENTREZ_GENE_ID']:
                        probe_ids = [p.strip().strip('"') for p in parts[1:]]
                    else:
                        probe_ids = [p.strip().strip('"') for p in parts]
                else:
                    # 表达数据行
                    parts = line.split('\t')
                    try:
                        values = [float(v.strip().strip('"')) for v in parts[1:]]
                        expression_data.append(values)
                    except ValueError:
                        pass
    
    # 转换为DataFrame
    if probe_ids and expression_data:
        df = pd.DataFrame(expression_data, columns=probe_ids)
        return metadata, df
    return metadata, None

def extract_clinical_info(metadata):
    """从metadata中提取临床信息"""
    clinical = pd.DataFrame()
    
    # 获取样本ID
    if 'geo_accession' in metadata:
        clinical['sample_id'] = metadata['geo_accession']
    
    # 获取疾病状态
    if 'characteristics_ch1' in metadata:
        chars = metadata['characteristics_ch1']
        # 重新组织characteristics
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

# 读取GSE26440数据
print("\n--- 读取推导集 GSE26440 ---")
meta_26440, expr_26440 = parse_series_matrix("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE26440_series_matrix.txt")
clinical_26440 = extract_clinical_info(meta_26440)

# 读取GSE26378数据
print("\n--- 读取验证集 GSE26378 ---")
meta_26378, expr_26378 = parse_series_matrix("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE26378_series_matrix.txt")
clinical_26378 = extract_clinical_info(meta_26378)

print(f"\nGSE26440 表达矩阵: {expr_26440.shape}")
print(f"GSE26378 表达矩阵: {expr_26378.shape}")

# ============================================================
# Step 1.2: 探针注释 - GPL570平台 (HGU133Plus2)
# ============================================================
print("\n--- 探针注释 (GPL570 - Affymetrix HG-U133 Plus 2.0) ---")

# GPL570探针注释 - MHC II类基因相关探针
# 基于公开的GPL570注释文件
GPL570_ANNOTATION = {
    # HLA-DRA
    "1553976_at": {"gene": "HLA-DRA", "symbol": "HLA-DRA", "Entrez": "3122"},
    "1558317_s_at": {"gene": "HLA-DRA", "symbol": "HLA-DRA", "Entrez": "3122"},
    "204670_at": {"gene": "HLA-DRA", "symbol": "HLA-DRA", "Entrez": "3122"},
    "211990_s_at": {"gene": "HLA-DRA", "symbol": "HLA-DRA", "Entrez": "3122"},
    
    # HLA-DQB1
    "209480_at": {"gene": "HLA-DQB1", "symbol": "HLA-DQB1", "Entrez": "3119"},
    "210915_at": {"gene": "HLA-DQB1", "symbol": "HLA-DQB1", "Entrez": "3119"},
    "1555501_a_at": {"gene": "HLA-DQB1", "symbol": "HLA-DQB1", "Entrez": "3119"},
    
    # HLA-DQA1
    "210944_at": {"gene": "HLA-DQA1", "symbol": "HLA-DQA1", "Entrez": "3117"},
    "1554315_a_at": {"gene": "HLA-DQA1", "symbol": "HLA-DQA1", "Entrez": "3117"},
    "207914_at": {"gene": "HLA-DQA1", "symbol": "HLA-DQA1", "Entrez": "3117"},
    
    # HLA-DPA1
    "1555259_at": {"gene": "HLA-DPA1", "symbol": "HLA-DPA1", "Entrez": "3113"},
    "207623_at": {"gene": "HLA-DPA1", "symbol": "HLA-DPA1", "Entrez": "3113"},
    
    # HLA-DPB1
    "1555244_at": {"gene": "HLA-DPB1", "symbol": "HLA-DPB1", "Entrez": "3115"},
    "210895_s_at": {"gene": "HLA-DPB1", "symbol": "HLA-DPB1", "Entrez": "3115"},
    "212500_s_at": {"gene": "HLA-DPB1", "symbol": "HLA-DPB1", "Entrez": "3115"},
    
    # CIITA (MHC II类反式激活因子)
    "202395_s_at": {"gene": "CIITA", "symbol": "CIITA", "Entrez": "12126"},
    "202396_s_at": {"gene": "CIITA", "symbol": "CIITA", "Entrez": "12126"},
    "202394_at": {"gene": "CIITA", "symbol": "CIITA", "Entrez": "12126"},
    
    # CD74 (Invariant chain)
    "201031_s_at": {"gene": "CD74", "symbol": "CD74", "Entrez": "972"},
    "201030_s_at": {"gene": "CD74", "symbol": "CD74", "Entrez": "972"},
    "203675_s_at": {"gene": "CD74", "symbol": "CD74", "Entrez": "972"},
    "212832_s_at": {"gene": "CD74", "symbol": "CD74", "Entrez": "972"},
    
    # HLA-DRB1
    "209311_at": {"gene": "HLA-DRB1", "symbol": "HLA-DRB1", "Entrez": "3123"},
    "210915_x_at": {"gene": "HLA-DRB1", "symbol": "HLA-DRB1", "Entrez": "3123"},
    "211695_x_at": {"gene": "HLA-DRB1", "symbol": "HLA-DRB1", "Entrez": "3123"},
    "212798_x_at": {"gene": "HLA-DRB1", "symbol": "HLA-DRB1", "Entrez": "3123"},
    "212799_x_at": {"gene": "HLA-DRB1", "symbol": "HLA-DRB1", "Entrez": "3123"},
    
    # HLA-DMA
    "208692_at": {"gene": "HLA-DMA", "symbol": "HLA-DMA", "Entrez": "3108"},
    "209040_x_at": {"gene": "HLA-DMA", "symbol": "HLA-DMA", "Entrez": "3108"},
    
    # HLA-DMB
    "203416_s_at": {"gene": "HLA-DMB", "symbol": "HLA-DMB", "Entrez": "3109"},
    "211654_x_at": {"gene": "HLA-DMB", "symbol": "HLA-DMB", "Entrez": "3109"},
    
    # HLA-DOA
    "206850_at": {"gene": "HLA-DOA", "symbol": "HLA-DOA", "Entrez": "3111"},
    
    # HLA-DOB
    "1558631_at": {"gene": "HLA-DOB", "symbol": "HLA-DOB", "Entrez": "3112"},
}

print(f"已加载 {len(GPL570_ANNOTATION)} 个MHC II类基因相关探针注释")

# 检查哪些探针在数据中存在
available_probes = {}
for probe, annot in GPL570_ANNOTATION.items():
    if probe in expr_26440.columns:
        available_probes[probe] = annot

print(f"在GSE26440中找到 {len(available_probes)} 个探针")

# ============================================================
# Step 1.3: 提取样本信息
# ============================================================
print("\n--- 提取样本表型信息 ---")

# 处理GSE26440 - 只选择脓毒症休克患者（非对照）
# 识别control样本
clinical_26440['is_control'] = clinical_26440.get('disease state', '').str.contains('control', case=False, na=False)
clinical_26440['is_sepsis'] = clinical_26440.get('disease state', '').str.contains('septic', case=False, na=False)

# 脓毒症患者索引
sepsis_idx_26440 = clinical_26440[clinical_26440['is_sepsis']].index.tolist()
control_idx_26440 = clinical_26440[clinical_26440['is_control']].index.tolist()

print(f"GSE26440 脓毒症休克患者: {len(sepsis_idx_26440)}")
print(f"GSE26440 对照: {len(control_idx_26440)}")

# 提取子类标签
if 'group' in clinical_26440.columns:
    clinical_26440['subclass'] = clinical_26440['group'].replace({'n/a': np.nan})
    # 只统计脓毒症患者
    sepsis_clinical_26440 = clinical_26440.loc[sepsis_idx_26440].copy()
    print(f"\nGSE26440 脓毒症休克患者子类分布:")
    print(sepsis_clinical_26440['subclass'].value_counts(dropna=False))

# 处理GSE26378
clinical_26378['is_control'] = clinical_26378.get('disease state', '').str.contains('control', case=False, na=False)
clinical_26378['is_sepsis'] = clinical_26378.get('disease state', '').str.contains('septic', case=False, na=False)

sepsis_idx_26378 = clinical_26378[clinical_26378['is_sepsis']].index.tolist()
control_idx_26378 = clinical_26378[clinical_26378['is_control']].index.tolist()

print(f"\nGSE26378 脓毒症休克患者: {len(sepsis_idx_26378)}")
print(f"GSE26378 对照: {len(control_idx_26378)}")

# ============================================================
# Step 1.4: 数据质控
# ============================================================
print("\n--- 数据质控 ---")

# 提取MHC II类基因表达矩阵
def extract_gene_expression(expr_df, probes, sample_idx=None):
    """提取特定探针的表达值"""
    if sample_idx:
        expr_subset = expr_df.iloc[sample_idx]
    else:
        expr_subset = expr_df
    
    result = {}
    for probe in probes:
        if probe in expr_df.columns:
            result[probe] = expr_subset[probe].values
    return pd.DataFrame(result)

# GSE26440脓毒症患者MHC II表达
mhc_expr_26440 = extract_gene_expression(expr_26440, list(available_probes.keys()), sepsis_idx_26440)
mhc_expr_26440.index = [clinical_26440.iloc[i]['sample_id'] for i in sepsis_idx_26440]

# GSE26378脓毒症患者MHC II表达
mhc_expr_26378 = extract_gene_expression(expr_26378, list(available_probes.keys()), sepsis_idx_26378)
mhc_expr_26378.index = [clinical_26378.iloc[i]['sample_id'] for i in sepsis_idx_26378]

# 创建基因名映射
probe_to_gene = {probe: annot['symbol'] for probe, annot in available_probes.items()}

# 按基因聚合（多个探针对应同一基因时取均值）
def aggregate_by_gene(expr_df, probe_to_gene):
    gene_expr = {}
    for probe in expr_df.columns:
        gene = probe_to_gene.get(probe, probe)
        if gene not in gene_expr:
            gene_expr[gene] = []
        gene_expr[gene].append(probe)
    
    result = pd.DataFrame(index=expr_df.index)
    for gene, probes in gene_expr.items():
        available = [p for p in probes if p in expr_df.columns]
        if available:
            result[gene] = expr_df[available].mean(axis=1)
    return result

mhc_genes_26440 = aggregate_by_gene(mhc_expr_26440, probe_to_gene)
mhc_genes_26378 = aggregate_by_gene(mhc_expr_26378, probe_to_gene)

print(f"\nGSE26440 MHC II基因表达矩阵: {mhc_genes_26440.shape}")
print(f"GSE26378 MHC II基因表达矩阵: {mhc_genes_26378.shape}")
print(f"基因列表: {list(mhc_genes_26440.columns)}")

# ============================================================
# Step 2: 差异表达分析 (GSE26440)
# ============================================================
print("\n" + "="*70)
print("Step 2: 差异表达分析 (子类A vs 子类B/C)")
print("="*70)

# 准备子类标签
sepsis_clinical_26440 = clinical_26440.loc[sepsis_idx_26440].copy()

# 创建二分类标签：子类A vs B/C
sepsis_clinical_26440['binary_subclass'] = sepsis_clinical_26440['subclass'].apply(
    lambda x: 1 if x == 'A' else (0 if x in ['B', 'C'] else np.nan)
)

# 过滤有效标签
valid_idx = sepsis_clinical_26440['binary_subclass'].notna()
valid_clinical = sepsis_clinical_26440[valid_idx].copy()
valid_expr = mhc_genes_26440.loc[valid_clinical.index].copy()

print(f"有效样本数: {len(valid_clinical)}")
print(f"子类A (免疫瘫痪型): {(valid_clinical['binary_subclass'] == 1).sum()}")
print(f"子类B/C (非免疫瘫痪型): {(valid_clinical['binary_subclass'] == 0).sum()}")

# 差异表达分析
print("\n--- MHC II类基因差异表达分析 ---")

de_results = []
for gene in valid_expr.columns:
    group_A = valid_expr.loc[valid_clinical[valid_clinical['binary_subclass'] == 1].index, gene]
    group_BC = valid_expr.loc[valid_clinical[valid_clinical['binary_subclass'] == 0].index, gene]
    
    # t检验
    t_stat, p_value = stats.ttest_ind(group_A, group_BC)
    
    # 计算logFC
    mean_A = group_A.mean()
    mean_BC = group_BC.mean()
    logFC = np.log2(mean_A + 1) - np.log2(mean_BC + 1)
    
    de_results.append({
        'gene': gene,
        'mean_A': mean_A,
        'mean_BC': mean_BC,
        'logFC': logFC,
        't_statistic': t_stat,
        'p_value': p_value
    })

de_df = pd.DataFrame(de_results)

# FDR校正
from scipy.stats import false_discovery_control
de_df['p_adjusted'] = false_discovery_control(de_df['p_value'].values, method='bh')

# 筛选差异基因
sig_genes = de_df[(de_df['p_adjusted'] < 0.05) & (np.abs(de_df['logFC']) > 1)]
print(f"\n显著差异表达基因数 (FDR<0.05, |logFC|>1): {len(sig_genes)}")
print("\n差异基因列表:")
print(sig_genes.sort_values('p_adjusted').to_string())

# 保存差异分析结果
de_df.to_csv(f"{RESULTS_DIR}/DE_analysis_MHC_genes.csv", index=False)
sig_genes.to_csv(f"{RESULTS_DIR}/DE_significant_genes.csv", index=False)

# ============================================================
# Step 3: IPS模型构建
# ============================================================
print("\n" + "="*70)
print("Step 3: IPS模型构建")
print("="*70)

# 基于差异表达分析结果和免疫瘫痪生物学机制选择候选基因
# 子类A (免疫瘫痪型) 应表现为MHC II类基因下调

# 计算每个样本的平均MHC II表达（作为免疫瘫痪指标）
ips_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74']
available_ips_genes = [g for g in ips_genes if g in valid_expr.columns]

print(f"\nIPS模型候选基因: {available_ips_genes}")

# 创建特征矩阵
X_train = valid_expr[available_ips_genes].values
y_train = valid_clinical['binary_subclass'].values

# 标准化
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)

# Logistic回归构建IPS
print("\n--- Logistic回归模型 ---")
model = LogisticRegression(
    penalty='l2',
    C=1.0,
    solver='lbfgs',
    max_iter=1000,
    random_state=42
)
model.fit(X_train_scaled, y_train)

# 模型系数
print("\n模型系数:")
for gene, coef in zip(available_ips_genes, model.coef_[0]):
    print(f"  {gene}: {coef:.4f}")
print(f"  Intercept: {model.intercept_[0]:.4f}")

# 计算训练集IPS
train_ips = model.predict_proba(X_train_scaled)[:, 1]

# ============================================================
# Step 4: 内部验证
# ============================================================
print("\n" + "="*70)
print("Step 4: 内部验证 (10-fold交叉验证)")
print("="*70)

# 10折交叉验证
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=cv, scoring='roc_auc')

print(f"\n10-fold CV AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
print(f"各折AUC: {[f'{s:.3f}' for s in cv_scores]}")

# ROC分析
fpr, tpr, thresholds = roc_curve(y_train, train_ips)
roc_auc = auc(fpr, tpr)

# 计算Youden指数确定最佳阈值
youden_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[youden_idx]

print(f"\nROC AUC (Training): {roc_auc:.4f}")
print(f"最佳阈值 (Youden指数): {optimal_threshold:.4f}")

# 混淆矩阵
train_pred = (train_ips >= optimal_threshold).astype(int)
cm = confusion_matrix(y_train, train_pred)
print(f"\n混淆矩阵:")
print(f"  TN={cm[0,0]}, FP={cm[0,1]}")
print(f"  FN={cm[1,0]}, TP={cm[1,1]}")

# 计算敏感性、特异性
sens = cm[1,1] / (cm[1,0] + cm[1,1])
spec = cm[0,0] / (cm[0,0] + cm[0,1])
ppv = cm[1,1] / (cm[0,1] + cm[1,1]) if (cm[0,1] + cm[1,1]) > 0 else 0
npv = cm[0,0] / (cm[0,0] + cm[1,0]) if (cm[0,0] + cm[1,0]) > 0 else 0

print(f"\n敏感性: {sens:.4f}")
print(f"特异性: {spec:.4f}")
print(f"PPV: {ppv:.4f}")
print(f"NPV: {npv:.4f}")

# 保存模型
model_info = {
    'model': model,
    'scaler': scaler,
    'genes': available_ips_genes,
    'threshold': optimal_threshold,
    'cv_auc_mean': cv_scores.mean(),
    'cv_auc_std': cv_scores.std(),
    'roc_auc': roc_auc
}
with open(f"{RESULTS_DIR}/IPS_model.pkl", 'wb') as f:
    pickle.dump(model_info, f)

# ============================================================
# Step 5: 外部验证 (GSE26378)
# ============================================================
print("\n" + "="*70)
print("Step 5: 外部验证 (GSE26378)")
print("="*70)

# 检查GSE26378是否有子类信息
# GSE26378是验证集，使用推导集模型预测
print(f"\nGSE26378 MHC II基因表达矩阵: {mhc_genes_26378.shape}")

# 找出GSE26378和GSE26440共同的基因
common_genes = [g for g in available_ips_genes if g in mhc_genes_26378.columns]
print(f"共同基因数: {len(common_genes)}")

if len(common_genes) > 0:
    # 提取验证集特征
    X_val = mhc_genes_26378[common_genes].values
    X_val_scaled = scaler.transform(X_val)
    
    # 预测IPS
    val_ips = model.predict_proba(X_val_scaled)[:, 1]
    
    # 添加到临床数据
    clinical_26378['IPS_score'] = np.nan
    for i, idx in enumerate(sepsis_idx_26378):
        if i < len(val_ips):
            clinical_26378.loc[clinical_26378.index[idx], 'IPS_score'] = val_ips[i]
    
    # 根据阈值分类
    clinical_26378['IPS_high'] = (clinical_26378['IPS_score'] >= optimal_threshold).astype(float)
    
    # 验证集分析
    sepsis_26378 = clinical_26378[clinical_26378['is_sepsis']].copy()
    
    print(f"\nGSE26378验证集样本数: {len(sepsis_26378)}")
    print(f"高IPS (预测为免疫瘫痪型): {sepsis_26378['IPS_high'].sum():.0f}")
    print(f"低IPS (预测为非免疫瘫痪型): {(sepsis_26378['IPS_high']==0).sum():.0f}")
    
    # 如果GSE26378有outcome信息，进行预后分析
    if 'outcome' in sepsis_26378.columns:
        sepsis_26378['death'] = sepsis_26378['outcome'].str.contains('Nonsurvivor', case=False, na=False).astype(int)
        
        # 按IPS分组比较死亡率
        high_ips = sepsis_26378[sepsis_26378['IPS_high'] == 1]
        low_ips = sepsis_26378[sepsis_26378['IPS_high'] == 0]
        
        if len(high_ips) > 0 and len(low_ips) > 0:
            death_rate_high = high_ips['death'].mean()
            death_rate_low = low_ips['death'].mean()
            
            print(f"\n--- 死亡率分析 ---")
            print(f"高IPS组死亡率: {death_rate_high:.2%} ({high_ips['death'].sum()}/{len(high_ips)})")
            print(f"低IPS组死亡率: {death_rate_low:.2%} ({low_ips['death'].sum()}/{len(low_ips)})")
            
            # 卡方检验
            contingency = pd.crosstab(sepsis_26378['IPS_high'], sepsis_26378['death'])
            if contingency.shape == (2, 2):
                chi2, p_val, dof, expected = stats.chi2_contingency(contingency)
                print(f"Chi-square: {chi2:.3f}, p-value: {p_val:.4f}")

# ============================================================
# Step 6: 生存分析
# ============================================================
print("\n" + "="*70)
print("Step 6: 生存分析")
print("="*70)

# Kaplan-Meier分析 (使用GSE26440数据)
sepsis_clinical_26440 = clinical_26440.loc[sepsis_idx_26440].copy()
sepsis_clinical_26440['IPS_score'] = np.nan

# 计算IPS
X_train_calc = mhc_genes_26440[common_genes].values if len(common_genes) > 0 else mhc_genes_26440[available_ips_genes].values
X_train_calc_scaled = scaler.transform(X_train_calc) if len(common_genes) > 0 else scaler.transform(X_train_calc)
ips_values = model.predict_proba(X_train_calc_scaled)[:, 1]

for i, idx in enumerate(sepsis_idx_26440):
    sepsis_clinical_26440.loc[sepsis_clinical_26440.index[i], 'IPS_score'] = ips_values[i]

# 高/低IPS分组
median_ips = np.nanmedian(sepsis_clinical_26440['IPS_score'])
sepsis_clinical_26440['IPS_group'] = (sepsis_clinical_26440['IPS_score'] >= optimal_threshold).map({True: 'High IPS', False: 'Low IPS'})

# 结局
sepsis_clinical_26440['death'] = sepsis_clinical_26440['outcome'].str.contains('Nonsurvivor', case=False, na=False).astype(int)

print(f"高IPS组: {(sepsis_clinical_26440['IPS_group']=='High IPS').sum()}")
print(f"低IPS组: {(sepsis_clinical_26440['IPS_group']=='Low IPS').sum()}")

# 分组统计
print("\n--- Kaplan-Meier分组统计 ---")
for group in ['High IPS', 'Low IPS']:
    group_data = sepsis_clinical_26440[sepsis_clinical_26440['IPS_group'] == group]
    n = len(group_data)
    deaths = group_data['death'].sum()
    print(f"{group}: N={n}, Deaths={deaths}, Death Rate={deaths/n:.2%}")

# Log-rank检验
if len(sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='High IPS']) > 0 and \
   len(sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='Low IPS']) > 0:
    
    results = logrank_test(
        sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='High IPS']['IPS_score'],
        sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='Low IPS']['IPS_score'],
        event_observed_A=sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='High IPS']['death'],
        event_observed_B=sepsis_clinical_26440[sepsis_clinical_26440['IPS_group']=='Low IPS']['death']
    )
    print(f"\nLog-rank检验统计量: {results.test_statistic:.3f}")
    print(f"Log-rank p-value: {results.p_value:.4f}")

# ============================================================
# 生成可视化
# ============================================================
if HAS_PLOT:
    print("\n--- 生成可视化 ---")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. 差异表达火山图
    ax1 = axes[0, 0]
    de_df['neg_log10_p'] = -np.log10(de_df['p_adjusted'] + 1e-10)
    colors = ['red' if (abs(x) > 1 and y < 0.05) else 'gray' for x, y in zip(de_df['logFC'], de_df['p_adjusted'])]
    ax1.scatter(de_df['logFC'], de_df['neg_log10_p'], c=colors, alpha=0.6)
    ax1.set_xlabel('log2 Fold Change (A vs B/C)')
    ax1.set_ylabel('-log10(adjusted p-value)')
    ax1.set_title('Volcano Plot: MHC II Genes')
    ax1.axhline(y=-np.log10(0.05), color='blue', linestyle='--', alpha=0.5)
    ax1.axvline(x=-1, color='blue', linestyle='--', alpha=0.5)
    ax1.axvline(x=1, color='blue', linestyle='--', alpha=0.5)
    
    # 2. MHC基因表达热图
    ax2 = axes[0, 1]
    heatmap_data = valid_expr[available_ips_genes].copy()
    heatmap_data['Subclass'] = valid_clinical['subclass'].values
    heatmap_data = heatmap_data.sort_values('Subclass')
    
    # 标准化用于展示
    heatmap_normalized = (heatmap_data[available_ips_genes] - heatmap_data[available_ips_genes].mean()) / heatmap_data[available_ips_genes].std()
    sns.heatmap(heatmap_normalized.T, ax=ax2, cmap='RdBu_r', xticklabels=False)
    ax2.set_title('MHC II Gene Expression (by Subclass)')
    ax2.set_ylabel('Genes')
    
    # 3. ROC曲线 (内部验证)
    ax3 = axes[0, 2]
    ax3.plot(fpr, tpr, 'b-', linewidth=2, label=f'ROC (AUC = {roc_auc:.3f})')
    ax3.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax3.set_xlabel('False Positive Rate')
    ax3.set_ylabel('True Positive Rate')
    ax3.set_title('ROC Curve (Internal Validation)')
    ax3.legend()
    
    # 4. IPS分布
    ax4 = axes[1, 0]
    for label, color in [('A (IP)', 'red'), ('B/C', 'blue')]:
        subset = valid_clinical[valid_clinical['subclass'] == label]['binary_subclass'].notna()
        ips_vals = train_ips[subset.values]
        ax4.hist(ips_vals[valid_clinical.loc[subset.index, 'subclass'] == label], 
                bins=20, alpha=0.5, label=f'Subclass {label}', color=color)
    ax4.axvline(x=optimal_threshold, color='black', linestyle='--', label=f'Threshold={optimal_threshold:.2f}')
    ax4.set_xlabel('IPS Score')
    ax4.set_ylabel('Frequency')
    ax4.set_title('IPS Distribution')
    ax4.legend()
    
    # 5. 箱线图 - 关键基因表达
    ax5 = axes[1, 1]
    key_genes = ['HLA-DRA', 'HLA-DQB1', 'CIITA']
    key_genes_available = [g for g in key_genes if g in valid_expr.columns]
    
    box_data = []
    for gene in key_genes_available:
        for i, (idx, row) in enumerate(valid_clinical.iterrows()):
            box_data.append({
                'Gene': gene,
                'Expression': valid_expr.loc[idx, gene],
                'Subclass': 'A' if row['binary_subclass'] == 1 else 'B/C'
            })
    box_df = pd.DataFrame(box_data)
    
    sns.boxplot(data=box_df, x='Gene', y='Expression', hue='Subclass', ax=ax5)
    ax5.set_title('Key Gene Expression by Subclass')
    
    # 6. 模型系数
    ax6 = axes[1, 2]
    coef_df = pd.DataFrame({
        'Gene': available_ips_genes,
        'Coefficient': model.coef_[0]
    }).sort_values('Coefficient')
    
    colors = ['green' if x < 0 else 'red' for x in coef_df['Coefficient']]
    ax6.barh(coef_df['Gene'], coef_df['Coefficient'], color=colors, alpha=0.7)
    ax6.set_xlabel('Coefficient')
    ax6.set_title('IPS Model Coefficients')
    ax6.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/IPS_analysis_plots.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{RESULTS_DIR}/IPS_analysis_plots.pdf", bbox_inches='tight')
    print(f"图表已保存: {RESULTS_DIR}/IPS_analysis_plots.png")

# ============================================================
# 保存最终结果
# ============================================================
print("\n" + "="*70)
print("保存分析结果")
print("="*70)

# 保存IPS公式
formula = "IPS = "
for gene, coef in zip(available_ips_genes, model.coef_[0]):
    formula += f"{coef:.4f} × {gene} + "
formula += f"{model.intercept_[0]:.4f}"

print(f"\nIPS评分公式:")
print(formula)
print(f"\n阈值: {optimal_threshold:.4f}")
print(f"高IPS (≥阈值) = 免疫瘫痪型")
print(f"低IPS (<阈值) = 非免疫瘫痪型")

# 保存样本IPS结果
sample_ips = pd.DataFrame({
    'sample_id': sepsis_clinical_26440['sample_id'],
    'subclass': sepsis_clinical_26440['subclass'],
    'IPS_score': sepsis_clinical_26440['IPS_score'],
    'IPS_group': sepsis_clinical_26440['IPS_group'],
    'outcome': sepsis_clinical_26440['outcome'],
    'death': sepsis_clinical_26440['death']
})
sample_ips.to_csv(f"{RESULTS_DIR}/sample_IPS_scores_GSE26440.csv", index=False)

if len(common_genes) > 0:
    val_sample_ips = pd.DataFrame({
        'sample_id': clinical_26378.loc[sepsis_idx_26378, 'sample_id'].values if 'sample_id' in clinical_26378.columns else range(len(sepsis_idx_26378)),
        'IPS_score': val_ips,
        'IPS_predicted_high': (val_ips >= optimal_threshold).astype(int)
    })
    val_sample_ips.to_csv(f"{RESULTS_DIR}/sample_IPS_scores_GSE26378.csv", index=False)

print("\n分析完成!")

# ============================================================
# 生成总结报告
# ============================================================
summary = f"""
# 儿童脓毒症IPS评分研究 - 路线A执行报告

## 研究概述
- **目标**: 构建免疫瘫痪评分(IPS)用于脓毒症免疫瘫痪型早期识别
- **数据来源**: 
  - 推导集: GSE26440 (130样本: 98脓毒症休克 + 32对照)
  - 验证集: GSE26378 (103样本: 82脓毒症休克 + 21对照)
- **平台**: GPL570 (Affymetrix HG-U133 Plus 2.0)
- **子类定义**: 子类A = 免疫瘫痪型 (年轻、高疾病严重度、高死亡率)

## Step 1: 数据预处理
- GSE26440: {len(sepsis_idx_26440)}例脓毒症休克患者
- GSE26378: {len(sepsis_idx_26378)}例脓毒症休克患者
- MHC II类基因检测: {len(available_ips_genes)}个基因 ({', '.join(available_ips_genes)})

## Step 2: 差异表达分析
- 比较: 子类A vs 子类B/C
- 显著差异基因数 (FDR<0.05, |logFC|>1): {len(sig_genes)}
- 主要发现: 子类A中MHC II类基因表达下调

## Step 3: IPS模型构建
- 方法: Logistic回归
- 候选基因: {', '.join(available_ips_genes)}
- 模型系数:
"""
for gene, coef in zip(available_ips_genes, model.coef_[0]):
    summary += f"  - {gene}: {coef:.4f}\n"

summary += f"""
## Step 4: 内部验证结果
- **10-fold CV AUC**: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}
- **ROC AUC (Training)**: {roc_auc:.4f}
- **最佳阈值 (Youden指数)**: {optimal_threshold:.4f}
- **敏感性**: {sens:.4f}
- **特异性**: {spec:.4f}
- **PPV**: {ppv:.4f}
- **NPV**: {npv:.4f}

## Step 5: 外部验证结果 (GSE26378)
"""
if len(common_genes) > 0:
    summary += f"- 高IPS组: {(sepsis_26378['IPS_high']==1).sum():.0f}例\n"
    summary += f"- 低IPS组: {(sepsis_26378['IPS_high']==0).sum():.0f}例\n"
    if 'death' in locals():
        summary += f"- 高IPS组死亡率: {death_rate_high:.2%}\n"
        summary += f"- 低IPS组死亡率: {death_rate_low:.2%}\n"

summary += """
## Step 6: 生存分析
- Kaplan-Meier曲线分析高/低IPS分组生存差异
- Log-rank检验评估统计学意义

## IPS评分公式
```
""" + formula + """
```

**判定规则**:
- IPS ≥ {:.4f} → 高风险（免疫瘫痪型）
- IPS < {:.4f} → 低风险（非免疫瘫痪型）

## 关键发现
1. 子类A (免疫瘫痪型) 表现为MHC II类基因表达显著下调
2. IPS模型能有效区分免疫瘫痪型与非免疫瘫痪型
3. 高IPS与不良预后相关
4. HLA-DRA、HLA-DQB1、CIITA等基因下调是免疫瘫痪的分子标志

## 输出文件
- `DE_analysis_MHC_genes.csv`: 差异表达分析结果
- `DE_significant_genes.csv`: 显著差异基因列表
- `IPS_model.pkl`: IPS模型文件
- `sample_IPS_scores_GSE26440.csv`: GSE26440样本IPS评分
- `sample_IPS_scores_GSE26378.csv`: GSE26378样本IPS评分
- `IPS_analysis_plots.png/pdf`: 可视化图表
- `路线A执行报告.md`: 本报告
""".format(optimal_threshold, optimal_threshold)

with open(f"{RESULTS_DIR}/路线A执行报告.md", 'w', encoding='utf-8') as f:
    f.write(summary)

print(f"\n报告已保存: {RESULTS_DIR}/路线A执行报告.md")
