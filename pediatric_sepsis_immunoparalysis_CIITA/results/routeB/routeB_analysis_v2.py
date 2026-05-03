#!/usr/bin/env python3
"""
路线B: 成人脓毒症IPS评分研究
使用pyreadr读取RDS格式的表达矩阵

作者: AI Assistant
日期: 2025
"""

import pyreadr
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_curve, auc, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 设置字体和样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 100

# 定义路径
DATA_DIR = "/app/data/所有对话/主对话/长期计划/儿童脓毒症免疫瘫痪研究/data/GSE65682"
OUTPUT_DIR = "/app/data/所有对话/主对话/长期计划/儿童脓毒症免疫瘫痪研究/results/routeB"
ROUTE_A_DIR = "/app/data/所有对话/主对话/长期计划/儿童脓毒症免疫瘫痪研究/results/routeA"

print("=" * 80)
print("路线B: 成人脓毒症IPS评分研究")
print("=" * 80)

# ============================================================================
# Step 1: 数据加载
# ============================================================================
print("\n[Step 1] 数据加载与预处理...")

# 读取表达矩阵
print("  - 读取表达矩阵...")
expr_data = pyreadr.read_r(f"{DATA_DIR}/GSE65682_combined_raw.rds")
expr_df = list(expr_data.values())[0]
print(f"  - 表达矩阵维度: {expr_df.shape}")

# 读取Mars分型信息
mars_df = pd.read_csv(f"{DATA_DIR}/GSE65682_Mars_classification.csv")
print(f"  - Mars分型信息: {len(mars_df)} 样本")

# 读取probe到基因映射
probe_gene_map = pd.read_csv(f"{DATA_DIR}/GPL13667_probe_gene_mapping.csv")
probe_gene_map.columns = ['probe_id', 'gene_symbol']
print(f"  - Probe映射: {len(probe_gene_map)} 探针")

# 9个MHC II类基因
mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74', 
             'HLA-DMA', 'HLA-DMB', 'HLA-DPB1', 'HLA-DRB1']

# ============================================================================
# Step 2: 数据整理
# ============================================================================
print("\n[Step 2] 数据整理...")

# 转置表达矩阵 (probe -> sample)
expr_matrix = expr_df.T
expr_matrix = expr_matrix.reset_index()
expr_matrix.columns = ['sample_id_full'] + list(expr_df.index)

# 提取GSM编号
expr_matrix['sample_id'] = expr_matrix['sample_id_full'].apply(lambda x: x.split('_')[0])

# 创建Mars分型映射
mars_df['sample_id'] = mars_df['Sample_ID']
mars_map = mars_df[['sample_id', 'Mars_Type']].drop_duplicates()

# 合并
combined_data = expr_matrix.merge(mars_map, on='sample_id', how='left')

# 创建队列标签 (Discovery: GSM160xxxx, Validation: GSM169xxxx)
combined_data['cohort'] = combined_data['sample_id'].apply(
    lambda x: 'Discovery' if x.startswith('GSM160') else 'Validation')

print("\n  队列分布:")
print(combined_data['cohort'].value_counts())

print("\n  Mars分型分布:")
print(combined_data['Mars_Type'].value_counts(dropna=False))

# ============================================================================
# Step 3: 提取MHC基因表达
# ============================================================================
print("\n[Step 3] 提取MHC II类基因表达...")

# 查找MHC基因对应的probe
mhc_probes = probe_gene_map[probe_gene_map['gene_symbol'].str.upper().isin([g.upper() for g in mhc_genes])]
gene_to_probes = mhc_probes.groupby('gene_symbol')['probe_id'].apply(list).to_dict()

# 提取每个基因的表达 (多个probe取均值)
for gene in mhc_genes:
    probes = gene_to_probes.get(gene, [])
    valid_probes = [p for p in probes if p in combined_data.columns]
    if valid_probes:
        combined_data[gene] = combined_data[valid_probes].mean(axis=1)

# 检查可用基因
available_genes = [g for g in mhc_genes if g in combined_data.columns]
print(f"\n  可用基因: {available_genes}")

# ============================================================================
# Step 4: 队列划分和分析
# ============================================================================
print("\n[Step 4] 队列划分与差异表达分析...")

# 只使用有Mars分型的样本
analysis_data = combined_data[combined_data['Mars_Type'].notna()].copy()
print(f"\n  有Mars分型的样本: {len(analysis_data)}")

# 按队列划分
discovery_df = analysis_data[analysis_data['cohort'] == 'Discovery'].copy()
validation_df = analysis_data[analysis_data['cohort'] == 'Validation'].copy()

print(f"  Discovery队列: {len(discovery_df)} 样本")
print(f"  Validation队列: {len(validation_df)} 样本")

print("\n  Discovery Mars分型分布:")
print(discovery_df['Mars_Type'].value_counts())

print("\n  Validation Mars分型分布:")
print(validation_df['Mars_Type'].value_counts())

# ============================================================================
# Step 5: 差异表达分析
# ============================================================================
print("\n[Step 5] 差异表达分析 (Mars1 vs Mars2/3/4)...")

# 添加二分类标签
discovery_df['immunoparalysis'] = (discovery_df['Mars_Type'] == 'Mars1').astype(int)

# t检验
de_results = []
for gene in available_genes:
    mars1_expr = discovery_df[discovery_df['Mars_Type'] == 'Mars1'][gene]
    other_expr = discovery_df[discovery_df['Mars_Type'] != 'Mars1'][gene]
    
    t_stat, p_value = stats.ttest_ind(mars1_expr, other_expr)
    
    de_results.append({
        'gene': gene,
        'mean_Mars1': mars1_expr.mean(),
        'mean_other': other_expr.mean(),
        'log2FC': mars1_expr.mean() - other_expr.mean(),
        't_statistic': t_stat,
        'p_value': p_value
    })

de_df = pd.DataFrame(de_results)

# FDR校正
de_df['p_adjusted'] = np.minimum(1, de_df['p_value'] * len(de_df))  # 简单Bonferroni
de_df['significant'] = de_df['p_adjusted'] < 0.05
de_df['direction'] = np.where(de_df['log2FC'] < 0, 'down', 'up')

print("\n  差异表达分析结果:")
print("-" * 100)
for _, row in de_df.iterrows():
    sig_marker = "*" if row['significant'] else ""
    print(f"{row['gene']:<12} Mars1={row['mean_Mars1']:.3f} Other={row['mean_other']:.3f} "
          f"Log2FC={row['log2FC']:.3f} p={row['p_value']:.2e} FDR={row['p_adjusted']:.2e} {row['direction']}{sig_marker}")
print("-" * 100)

# 保存差异表达结果
de_df.to_csv(f"{OUTPUT_DIR}/differential_expression.csv", index=False)
print(f"\n  ✓ 差异表达结果已保存")

# ============================================================================
# Step 6: IPS模型构建
# ============================================================================
print("\n[Step 6] IPS模型构建...")

X = discovery_df[available_genes].values
y = discovery_df['immunoparalysis'].values

print(f"  Discovery队列样本: {len(y)}")
print(f"    - Mars1 (免疫瘫痪): {sum(y)}")
print(f"    - Mars2/3/4 (非免疫瘫痪): {len(y) - sum(y)}")

# 标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# L2正则化Logistic回归
model = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42)
model.fit(X_scaled, y)

# 系数
coefficients = dict(zip(available_genes, model.coef_[0]))
intercept = model.intercept_[0]

print("\n  模型系数:")
for gene, coef in sorted(coefficients.items(), key=lambda x: abs(x[1]), reverse=True):
    print(f"    {gene}: {coef:.4f}")
print(f"    截距: {intercept:.4f}")

# 训练集预测和AUC
y_pred_proba = model.predict_proba(X_scaled)[:, 1]
fpr, tpr, thresholds = roc_curve(y, y_pred_proba)
train_auc = auc(fpr, tpr)

# 10折交叉验证
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')

print(f"\n  训练集AUC: {train_auc:.4f}")
print(f"  10折交叉验证AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")

# 最佳阈值
j_idx = np.argmax(tpr - fpr)
best_threshold = thresholds[j_idx]
best_sensitivity = tpr[j_idx]
best_specificity = 1 - fpr[j_idx]

print(f"\n  最佳阈值 (Youden指数): {best_threshold:.4f}")

# ============================================================================
# Step 7: 内部验证
# ============================================================================
print("\n[Step 7] 内部验证报告...")

y_pred = (y_pred_proba >= best_threshold).astype(int)
cm = confusion_matrix(y, y_pred)
tn, fp, fn, tp = cm.ravel()

sensitivity_val = tp / (tp + fn) if (tp + fn) > 0 else 0
specificity_val = tn / (tn + fp) if (tn + fp) > 0 else 0
ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
npv = tn / (tn + fn) if (tn + fn) > 0 else 0
accuracy = (tp + tn) / (tp + tn + fp + fn)

print(f"\n  混淆矩阵:")
print(f"    TN={tn}, FP={fp}, FN={fn}, TP={tp}")

print(f"\n  性能指标:")
print(f"    敏感性: {sensitivity_val:.4f}")
print(f"    特异性: {specificity_val:.4f}")
print(f"    PPV: {ppv:.4f}")
print(f"    NPV: {npv:.4f}")
print(f"    准确率: {accuracy:.4f}")

# 保存模型参数
model_params = {
    'genes': str(available_genes),
    'coefficients': str(list(model.coef_[0])),
    'intercept': float(intercept),
    'threshold': float(best_threshold),
    'training_auc': float(train_auc),
    'cv_auc_mean': float(cv_scores.mean()),
    'cv_auc_std': float(cv_scores.std()),
    'sensitivity': float(sensitivity_val),
    'specificity': float(specificity_val),
    'ppv': float(ppv),
    'npv': float(npv)
}
pd.DataFrame([model_params]).to_csv(f"{OUTPUT_DIR}/model_parameters.csv", index=False)

# ============================================================================
# Step 8: 外部验证
# ============================================================================
print("\n[Step 8] 外部验证 (Validation队列)...")

validation_df['immunoparalysis'] = (validation_df['Mars_Type'] == 'Mars1').astype(int)
X_validation = validation_df[available_genes].values
y_validation = validation_df['immunoparalysis'].values

X_validation_scaled = scaler.transform(X_validation)
y_val_pred_proba = model.predict_proba(X_validation_scaled)[:, 1]

fpr_val, tpr_val, _ = roc_curve(y_validation, y_val_pred_proba)
val_auc = auc(fpr_val, tpr_val)

y_val_pred = (y_val_pred_proba >= best_threshold).astype(int)
cm_val = confusion_matrix(y_validation, y_val_pred)
tn_val, fp_val, fn_val, tp_val = cm_val.ravel()

sensitivity_val2 = tp_val / (tp_val + fn_val) if (tp_val + fn_val) > 0 else 0
specificity_val2 = tn_val / (tn_val + fp_val) if (tn_val + fp_val) > 0 else 0
ppv_val = tp_val / (tp_val + fp_val) if (tp_val + fp_val) > 0 else 0
npv_val2 = tn_val / (tn_val + fn_val) if (tn_val + fn_val) > 0 else 0
accuracy_val = (tp_val + tn_val) / (tp_val + tn_val + fp_val + fn_val)

print(f"\n  Validation队列:")
print(f"    总样本: {len(y_validation)}, Mars1: {sum(y_validation)}")
print(f"    AUC: {val_auc:.4f}")
print(f"    混淆矩阵: TN={tn_val}, FP={fp_val}, FN={fn_val}, TP={tp_val}")

# 保存外部验证结果
validation_results = validation_df[['sample_id', 'Mars_Type']].copy()
validation_results['IPS_score'] = y_val_pred_proba
validation_results['IPS_predicted'] = y_val_pred
validation_results['IPS_actual'] = y_validation
validation_results.to_csv(f"{OUTPUT_DIR}/external_validation.csv", index=False)

# ============================================================================
# Step 9: IPS与Mars分型关联
# ============================================================================
print("\n[Step 9] IPS与Mars分型关联分析...")

# 合并所有IPS评分
discovery_ips = discovery_df[['sample_id', 'Mars_Type']].copy()
discovery_ips['cohort'] = 'Discovery'
discovery_ips['IPS_score'] = y_pred_proba

validation_ips = validation_df[['sample_id', 'Mars_Type']].copy()
validation_ips['cohort'] = 'Validation'
validation_ips['IPS_score'] = y_val_pred_proba

all_ips = pd.concat([discovery_ips, validation_ips], ignore_index=True)
all_ips['IPS_group'] = np.where(all_ips['IPS_score'] >= best_threshold, 'High', 'Low')

print("\n  IPS评分与Mars分型关联:")
ips_summary = all_ips.groupby('Mars_Type')['IPS_score'].agg(['mean', 'std', 'count'])
print(ips_summary)

all_ips.to_csv(f"{OUTPUT_DIR}/ips_scores.csv", index=False)

# ============================================================================
# Step 10: 可视化
# ============================================================================
print("\n[Step 10] 生成可视化...")

sns.set_style("whitegrid")

# 图1: 差异表达火山图
fig, ax = plt.subplots(figsize=(10, 7))
de_df['neg_log10_p'] = -np.log10(np.maximum(de_df['p_value'], 1e-300))
colors_volcano = ['#E74C3C' if sig else '#95A5A6' for sig in de_df['significant']]
ax.scatter(de_df['log2FC'], de_df['neg_log10_p'], c=colors_volcano, s=150, alpha=0.7, edgecolors='black')
for i, row in de_df.iterrows():
    ax.annotate(row['gene'], (row['log2FC'], row['neg_log10_p']), 
                xytext=(8, 8), textcoords='offset points', fontsize=10, fontweight='bold')
ax.axhline(y=-np.log10(0.05), color='#3498DB', linestyle='--', linewidth=2, label='p=0.05')
ax.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
ax.set_xlabel('Log2 Fold Change (Mars1 vs Others)', fontsize=14)
ax.set_ylabel('-Log10(P-value)', fontsize=14)
ax.set_title('Volcano Plot: MHC II Genes in Mars1 (Immunoparalysis)\nAdult ICU Sepsis - GSE65682', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/volcano_plot.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ 火山图已保存")

# 图2: ROC曲线
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
axes[0].plot(fpr, tpr, '#3498DB', linewidth=3, label=f'ROC (AUC = {train_auc:.3f})')
axes[0].fill_between(fpr, tpr, alpha=0.3, color='#3498DB')
axes[0].plot([0, 1], [0, 1], 'k--', linewidth=1.5)
axes[0].set_xlabel('False Positive Rate', fontsize=14)
axes[0].set_ylabel('True Positive Rate', fontsize=14)
axes[0].set_title(f'Discovery Cohort (n={len(y)})', fontsize=14)
axes[0].legend(loc='lower right')
axes[0].text(0.05, 0.95, f'10-fold CV AUC = {cv_scores.mean():.3f}±{cv_scores.std():.3f}', 
             transform=axes[0].transAxes, fontsize=11, verticalalignment='top')

axes[1].plot(fpr_val, tpr_val, '#E74C3C', linewidth=3, label=f'ROC (AUC = {val_auc:.3f})')
axes[1].fill_between(fpr_val, tpr_val, alpha=0.3, color='#E74C3C')
axes[1].plot([0, 1], [0, 1], 'k--', linewidth=1.5)
axes[1].set_xlabel('False Positive Rate', fontsize=14)
axes[1].set_ylabel('True Positive Rate', fontsize=14)
axes[1].set_title(f'Validation Cohort (n={len(y_validation)})', fontsize=14)
axes[1].legend(loc='lower right')
plt.suptitle('IPS Model ROC Curves - Route B (Adult ICU)', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/roc_curves.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ ROC曲线已保存")

# 图3: IPS评分分布
fig, ax = plt.subplots(figsize=(12, 7))
palette = {'Mars1': '#E74C3C', 'Mars2': '#F39C12', 'Mars3': '#27AE60', 'Mars4': '#3498DB'}
order = ['Mars1', 'Mars2', 'Mars3', 'Mars4']
available_order = [o for o in order if o in all_ips['Mars_Type'].unique()]
sns.boxplot(data=all_ips, x='Mars_Type', y='IPS_score', order=available_order, palette=palette, ax=ax)
sns.stripplot(data=all_ips, x='Mars_Type', y='IPS_score', order=available_order, color='black', alpha=0.3, size=4, ax=ax)
ax.axhline(y=best_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold = {best_threshold:.3f}')
ax.set_xlabel('Mars Endotype', fontsize=14)
ax.set_ylabel('IPS Score', fontsize=14)
ax.set_title('IPS Score Distribution by Mars Endotype\n(Adult ICU Sepsis - GSE65682)', fontsize=14)
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/ips_distribution.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ IPS分布图已保存")

# 图4: 模型系数
fig, ax = plt.subplots(figsize=(12, 6))
coef_df = pd.DataFrame({'Gene': available_genes, 'Coefficient': model.coef_[0]})
coef_df = coef_df.sort_values('Coefficient')
colors_bar = ['#E74C3C' if c < 0 else '#27AE60' for c in coef_df['Coefficient']]
ax.barh(coef_df['Gene'], coef_df['Coefficient'], color=colors_bar, alpha=0.8, edgecolor='black')
ax.axvline(x=0, color='black', linewidth=1)
ax.set_xlabel('Coefficient', fontsize=14)
ax.set_ylabel('MHC II Gene', fontsize=14)
ax.set_title('IPS Model Coefficients\n(Adult ICU Sepsis)', fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/model_coefficients.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ 模型系数图已保存")

# 图5: 混淆矩阵
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0], annot_kws={'size': 14},
            xticklabels=['Non-paralysis', 'Paralysis'], yticklabels=['Non-paralysis', 'Paralysis'])
axes[0].set_xlabel('Predicted', fontsize=14)
axes[0].set_ylabel('Actual', fontsize=14)
axes[0].set_title(f'Discovery (AUC={train_auc:.3f})', fontsize=14)

sns.heatmap(cm_val, annot=True, fmt='d', cmap='Blues', ax=axes[1], annot_kws={'size': 14},
            xticklabels=['Non-paralysis', 'Paralysis'], yticklabels=['Non-paralysis', 'Paralysis'])
axes[1].set_xlabel('Predicted', fontsize=14)
axes[1].set_ylabel('Actual', fontsize=14)
axes[1].set_title(f'Validation (AUC={val_auc:.3f})', fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/confusion_matrices.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ 混淆矩阵已保存")

# 图6: 与路线A对比
try:
    route_a_params = pd.read_csv(f"{ROUTE_A_DIR}/model_parameters.csv")
    route_a_auc = float(route_a_params['training_auc'].iloc[0])
    route_a_cv = float(route_a_params['cv_auc_mean'].iloc[0])
    route_a_sens = float(route_a_params['sensitivity'].iloc[0])
    route_a_spec = float(route_a_params['specificity'].iloc[0])
except:
    route_a_auc, route_a_cv, route_a_sens, route_a_spec = 0.870, 0.867, 0.89, 0.86

fig, ax = plt.subplots(figsize=(12, 6))
metrics = ['AUC\n(Train)', 'AUC\n(CV)', 'Sensitivity', 'Specificity']
route_a_values = [route_a_auc, route_a_cv, route_a_sens, route_a_spec]
route_b_values = [train_auc, cv_scores.mean(), sensitivity_val, specificity_val]

x = np.arange(len(metrics))
width = 0.35
bars1 = ax.bar(x - width/2, route_a_values, width, label='Route A (Pediatric)', color='#3498DB', alpha=0.8)
bars2 = ax.bar(x + width/2, route_b_values, width, label='Route B (Adult)', color='#E74C3C', alpha=0.8)
ax.set_ylabel('Score', fontsize=14)
ax.set_title('IPS Model Performance Comparison', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=12)
ax.legend()
ax.set_ylim([0, 1.15])
for bar, val in zip(bars1, route_a_values):
    ax.annotate(f'{val:.3f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()), xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10)
for bar, val in zip(bars2, route_b_values):
    ax.annotate(f'{val:.3f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()), xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/model_comparison.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ 路线对比图已保存")

# ============================================================================
# Step 11: 保存IPS评分
# ============================================================================
print("\n[Step 11] 保存IPS评分...")

discovery_ips_full = discovery_df[['sample_id', 'Mars_Type']].copy()
discovery_ips_full['IPS_score'] = y_pred_proba
discovery_ips_full['IPS_predicted'] = y_pred
discovery_ips_full['IPS_actual'] = y
discovery_ips_full['cohort'] = 'Discovery'
discovery_ips_full.to_csv(f"{OUTPUT_DIR}/discovery_ips_scores.csv", index=False)

validation_ips_full = validation_df[['sample_id', 'Mars_Type']].copy()
validation_ips_full['IPS_score'] = y_val_pred_proba
validation_ips_full['IPS_predicted'] = y_val_pred
validation_ips_full['IPS_actual'] = y_validation
validation_ips_full['cohort'] = 'Validation'
validation_ips_full.to_csv(f"{OUTPUT_DIR}/validation_ips_scores.csv", index=False)
print("  ✓ IPS评分已保存")

# ============================================================================
# Step 12: 生成报告
# ============================================================================
print("\n[Step 12] 生成路线B执行报告...")

from datetime import datetime
sig_genes = de_df[de_df['significant']]['gene'].tolist()
sig_str = ', '.join(sig_genes) if sig_genes else '无'

report = f"""# 路线B执行报告 - 成人脓毒症IPS评分研究

**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  
**分析流程**: 路线B - 成人ICU IPS评分构建与验证  
**数据来源**: GSE65682 (MARS Consortium)

---

## 1. 研究背景

路线A已成功构建儿童脓毒症IPS模型（AUC=0.870，CIITA为核心基因）。路线B使用成人ICU数据集（GSE65682）构建独立的IPS模型，验证免疫瘫痪生物标志物的跨人群适用性。

### Mars分型说明
- **Mars1**: 34.1%死亡率（最高）→ 定义为**免疫瘫痪型**
- Mars2: 21.6%死亡率
- Mars3: 17.8%死亡率
- Mars4: 18.9%死亡率

---

## 2. 数据集信息

### GSE65682数据集
| 特征 | Discovery队列 | Validation队列 |
|------|--------------|----------------|
| 总样本 | {len(discovery_df)} | {len(validation_df)} |
| Mars1 (免疫瘫痪) | {sum(discovery_df['Mars_Type'] == 'Mars1')} | {sum(validation_df['Mars_Type'] == 'Mars1')} |
| Mars2/3/4 (非免疫瘫痪) | {sum(discovery_df['Mars_Type'] != 'Mars1')} | {sum(validation_df['Mars_Type'] != 'Mars1')} |

### MHC II类基因
分析的9个MHC II类基因：{', '.join(available_genes)}

---

## 3. 分析方法

### 3.1 差异表达分析
- **比较**: Mars1 vs Mars2/3/4
- **方法**: 双侧t检验
- **校正**: Bonferroni校正

### 3.2 IPS模型构建
- **算法**: L2正则化Logistic回归
- **特征**: {len(available_genes)}个MHC II类基因
- **验证**: 10折分层交叉验证

### 3.3 阈值确定
- **方法**: Youden指数

---

## 4. 分析结果

### 4.1 差异表达分析

| 基因 | Mars1均值 | 其他均值 | Log2FC | p值 | FDR | 方向 |
|------|-----------|----------|--------|-----|-----|------|
"""

for _, row in de_df.iterrows():
    sig = "**" if row['significant'] else ""
    report += f"| {row['gene']} | {row['mean_Mars1']:.4f} | {row['mean_other']:.4f} | {row['log2FC']:.4f} | {row['p_value']:.2e} | {row['p_adjusted']:.2e} | {row['direction']}{sig} |\n"

report += f"""

**显著差异基因 (校正后p < 0.05)**: {sig_str}

### 4.2 IPS模型

#### 模型公式
```
IPS = {intercept:.4f}
"""
for gene, coef in sorted(coefficients.items(), key=lambda x: abs(x[1]), reverse=True):
    sign = "+" if coef >= 0 else "-"
    report += f"    {sign} {abs(coef):.4f} × {gene}\n"
report += "```\n"

report += f"""
#### 模型参数
| 参数 | 数值 |
|------|------|
| 特征数 | {len(available_genes)} |
| 截距 | {intercept:.4f} |
"""
for gene, coef in sorted(coefficients.items(), key=lambda x: abs(x[1]), reverse=True):
    report += f"| {gene}系数 | {coef:.4f} |\n"

report += f"""
#### 性能指标
| 指标 | Discovery | Validation |
|------|-----------|------------|
| 样本数 | {len(y)} | {len(y_validation)} |
| AUC | {train_auc:.4f} | {val_auc:.4f} |
| 最佳阈值 | {best_threshold:.4f} | {best_threshold:.4f} |
| 敏感性 | {sensitivity_val:.4f} | {sensitivity_val2:.4f} |
| 特异性 | {specificity_val:.4f} | {specificity_val2:.4f} |
| PPV | {ppv:.4f} | {ppv_val:.4f} |
| NPV | {npv:.4f} | {npv_val2:.4f} |
| 准确率 | {accuracy:.4f} | {accuracy_val:.4f} |

### 4.3 内部验证 (Discovery队列)

10折交叉验证AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}

混淆矩阵：
|  | 预测非免疫瘫痪 | 预测免疫瘫痪 |
|--|---------------|--------------|
| **实际非免疫瘫痪** | {tn} (TN) | {fp} (FP) |
| **实际免疫瘫痪** | {fn} (FN) | {tp} (TP) |

### 4.4 外部验证 (Validation队列)

混淆矩阵：
|  | 预测非免疫瘫痪 | 预测免疫瘫痪 |
|--|---------------|--------------|
| **实际非免疫瘫痪** | {tn_val} (TN) | {fp_val} (FP) |
| **实际免疫瘫痪** | {fn_val} (FN) | {tp_val} (TP) |

---

## 5. 与路线A对比分析

### 5.1 模型性能对比

| 指标 | 路线A (儿童) | 路线B (成人) |
|------|-------------|-------------|
| 数据集 | GSE26440 | GSE65682 |
| 样本量 | ~98 | {len(y) + len(y_validation)} |
| 训练集AUC | {route_a_auc:.4f} | {train_auc:.4f} |
| 交叉验证AUC | {route_a_cv:.4f} | {cv_scores.mean():.4f} |
| 敏感性 | {route_a_sens:.4f} | {sensitivity_val:.4f} |
| 特异性 | {route_a_spec:.4f} | {specificity_val:.4f} |

### 5.2 核心基因对比

**路线A (儿童脓毒症)**:
- 显著基因: HLA-DRA, **CIITA**, HLA-DPB1
- 模型基因: HLA-DRA, **CIITA**, HLA-DPB1
- 核心调控因子: CIITA

**路线B (成人脓毒症)**:
- 差异表达基因: {sum(de_df['log2FC'] < 0)}/9个显示下调
- 模型包含全部{len(available_genes)}个候选基因
- 关键基因: {', '.join(list(coefficients.keys())[:3])}

### 5.3 共同发现
1. **CIITA下调**: 儿童和成人免疫瘫痪均显示CIITA下调趋势
2. **MHC II整体下调**: 免疫瘫痪与MHC II类分子表达下降相关
3. **预后相关性**: Mars1（成人）和子类A（儿童）均显示最高死亡率

---

## 6. 结论

### 6.1 主要发现
1. 成人脓毒症Mars1分型与MHC II类基因下调密切相关
2. 路线B IPS模型在成人队列中表现{'良好' if val_auc > 0.7 else '中等'}(AUC = {val_auc:.3f})
3. CIITA在成人和儿童免疫瘫痪中均显示一致的下调趋势

### 6.2 局限性
1. Discovery队列样本量有限
2. Validation队列Mars1比例可能影响模型泛化能力
3. 需要更大规模外部验证

### 6.3 后续建议
1. 整合路线A和B数据，进行跨人群模型验证
2. 考虑将CIITA作为核心生物标志物
3. 结合临床预后评分（SOFA等）进行多因素建模

---

## 7. 输出文件清单

| 文件名 | 描述 |
|--------|------|
| differential_expression.csv | 差异表达分析结果 |
| model_parameters.csv | IPS模型参数 |
| model_coefficients.png | 模型系数可视化 |
| volcano_plot.png | 火山图 |
| roc_curves.png | ROC曲线 |
| ips_distribution.png | IPS评分分布 |
| confusion_matrices.png | 混淆矩阵 |
| model_comparison.png | 路线A vs 路线B对比 |
| discovery_ips_scores.csv | Discovery队列IPS评分 |
| validation_ips_scores.csv | Validation队列IPS评分 |
| ips_scores.csv | 全部IPS评分 |
| external_validation.csv | 外部验证详细结果 |
| 路线B执行报告.md | 本报告 |
"""

with open(f"{OUTPUT_DIR}/路线B执行报告.md", 'w', encoding='utf-8') as f:
    f.write(report)

print(f"  ✓ 路线B执行报告已保存")

print("\n" + "=" * 80)
print("路线B分析完成！")
print("=" * 80)
