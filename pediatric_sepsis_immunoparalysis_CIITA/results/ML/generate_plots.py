#!/usr/bin/env python3
"""
生成可视化图表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

BASE_DIR = '/app/data/所有对话/主对话/长期计划/儿童脓毒症免疫瘫痪研究'

# 读取数据
roc_data = pd.read_csv(f'{BASE_DIR}/results/ML/roc_curve_data.csv')
auc_data = pd.read_csv(f'{BASE_DIR}/results/ML/auc_values.csv')
ips_scores = pd.read_csv(f'{BASE_DIR}/results/ML/ips_scores_all.csv')
ips_weights = pd.read_csv(f'{BASE_DIR}/results/ML/ips_weights.csv')
feature_importance = pd.read_csv(f'{BASE_DIR}/results/ML/feature_importance.csv')
model_performance = pd.read_csv(f'{BASE_DIR}/results/ML/model_performance.csv')

# 读取探针到基因的映射
# 从DEGs文件中获取
degs = pd.read_csv(f'{BASE_DIR}/results/DEA/DEGs_significant_GSE26378.csv')
probe_to_gene = dict(zip(degs['ProbeID'], degs['GeneSymbol']))

# 转换探针ID为基因符号
ips_weights['GeneSymbol'] = ips_weights['Gene'].map(probe_to_gene)
ips_weights['GeneSymbol'] = ips_weights['GeneSymbol'].fillna(ips_weights['Gene'])

# ============= 图1: ROC曲线 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 训练集ROC
train_roc = roc_data[roc_data['Dataset'] == 'Training']
axes[0].plot(train_roc['FPR'], train_roc['TPR'], 'b-', linewidth=2, label=f'ROC (AUC = {auc_data[auc_data["Dataset"]=="Training"]["AUC"].values[0]:.4f})')
axes[0].plot([0, 1], [0, 1], 'r--', linewidth=1, label='Random')
axes[0].set_xlabel('False Positive Rate', fontsize=12)
axes[0].set_ylabel('True Positive Rate', fontsize=12)
axes[0].set_title('ROC Curve - Training Set (GSE26378)', fontsize=14)
axes[0].legend(loc='lower right')
axes[0].grid(True, alpha=0.3)

# 验证集ROC
val_roc = roc_data[roc_data['Dataset'] == 'Validation']
axes[1].plot(val_roc['FPR'], val_roc['TPR'], 'g-', linewidth=2, label=f'ROC (AUC = {auc_data[auc_data["Dataset"]=="Validation"]["AUC"].values[0]:.4f})')
axes[1].plot([0, 1], [0, 1], 'r--', linewidth=1, label='Random')
axes[1].set_xlabel('False Positive Rate', fontsize=12)
axes[1].set_ylabel('True Positive Rate', fontsize=12)
axes[1].set_title('ROC Curve - Validation Set (GSE26440)', fontsize=14)
axes[1].legend(loc='lower right')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/roc_curves.png', dpi=150, bbox_inches='tight')
plt.close()
print("ROC曲线已保存")

# ============= 图2: IPS分布 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 训练集IPS分布
train_ips = ips_scores[ips_scores['dataset'] == 'GSE26378']
for group in ['normal', 'sepsis']:
    data = train_ips[train_ips['group'] == group]['IPS_normalized']
    axes[0].hist(data, bins=20, alpha=0.6, label=f'{group} (n={len(data)})')
axes[0].axvline(x=0.0387, color='red', linestyle='--', linewidth=2, label='Threshold (0.0387)')
axes[0].set_xlabel('IPS Score (Normalized)', fontsize=12)
axes[0].set_ylabel('Frequency', fontsize=12)
axes[0].set_title('IPS Distribution - Training Set (GSE26378)', fontsize=14)
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# 验证集IPS分布
val_ips = ips_scores[ips_scores['dataset'] == 'GSE26440']
for group in ['normal', 'sepsis']:
    data = val_ips[val_ips['group'] == group]['IPS_normalized']
    axes[1].hist(data, bins=20, alpha=0.6, label=f'{group} (n={len(data)})')
axes[1].axvline(x=0.0387, color='red', linestyle='--', linewidth=2, label='Threshold (0.0387)')
axes[1].set_xlabel('IPS Score (Normalized)', fontsize=12)
axes[1].set_ylabel('Frequency', fontsize=12)
axes[1].set_title('IPS Distribution - Validation Set (GSE26440)', fontsize=14)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS分布图已保存")

# ============= 图3: 特征重要性 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# LASSO系数
top_lasso = feature_importance.nlargest(15, 'LASSO_coef')
genes_lasso = [probe_to_gene.get(g, g) for g in top_lasso.index]
axes[0].barh(range(len(top_lasso)), top_lasso['LASSO_coef'], color='steelblue')
axes[0].set_yticks(range(len(top_lasso)))
axes[0].set_yticklabels(genes_lasso, fontsize=9)
axes[0].set_xlabel('LASSO Coefficient', fontsize=12)
axes[0].set_title('Top 15 Features - LASSO', fontsize=14)
axes[0].invert_yaxis()
axes[0].grid(True, alpha=0.3, axis='x')

# 随机森林重要性
top_rf = feature_importance.nlargest(15, 'RF_importance')
genes_rf = [probe_to_gene.get(g, g) for g in top_rf.index]
axes[1].barh(range(len(top_rf)), top_rf['RF_importance'], color='forestgreen')
axes[1].set_yticks(range(len(top_rf)))
axes[1].set_yticklabels(genes_rf, fontsize=9)
axes[1].set_xlabel('Random Forest Importance', fontsize=12)
axes[1].set_title('Top 15 Features - Random Forest', fontsize=14)
axes[1].invert_yaxis()
axes[1].grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/feature_importance.png', dpi=150, bbox_inches='tight')
plt.close()
print("特征重要性图已保存")

# ============= 图4: IPS权重柱状图 =============
fig, ax = plt.subplots(figsize=(10, 6))
genes = ips_weights['GeneSymbol'].values
weights = ips_weights['Weight'].values
colors = ['steelblue' if w > 0 else 'coral' for w in weights]
ax.barh(range(len(genes)), weights, color=colors)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=10)
ax.set_xlabel('IPS Weight (Logistic Regression Coefficient)', fontsize=12)
ax.set_title('IPS Gene Weights', fontsize=14)
ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
ax.invert_yaxis()
ax.grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_weights.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS权重图已保存")

# ============= 图5: 模型性能对比 =============
fig, ax = plt.subplots(figsize=(10, 6))
metrics = ['Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1', 'AUC']
train_vals = [model_performance[model_performance['Model'].str.contains('Training')][m].values[0] for m in metrics]
val_vals = [model_performance[model_performance['Model'].str.contains('Validation')][m].values[0] for m in metrics]

x = np.arange(len(metrics))
width = 0.35
bars1 = ax.bar(x - width/2, train_vals, width, label='Training', color='steelblue')
bars2 = ax.bar(x + width/2, val_vals, width, label='Validation', color='forestgreen')

ax.set_ylabel('Score', fontsize=12)
ax.set_title('Model Performance Comparison', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(metrics)
ax.legend()
ax.set_ylim([0, 1.1])
ax.grid(True, alpha=0.3, axis='y')

# 添加数值标签
for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8)
for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/model_performance.png', dpi=150, bbox_inches='tight')
plt.close()
print("模型性能对比图已保存")

# ============= 图6: 箱线图比较 =============
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 训练集箱线图
train_ips = ips_scores[ips_scores['dataset'] == 'GSE26378']
data_normal = train_ips[train_ips['group'] == 'normal']['IPS_normalized']
data_sepsis = train_ips[train_ips['group'] == 'sepsis']['IPS_normalized']
bp1 = axes[0].boxplot([data_normal, data_sepsis], labels=['Normal', 'Sepsis'], patch_artist=True)
bp1['boxes'][0].set_facecolor('lightblue')
bp1['boxes'][1].set_facecolor('lightcoral')
axes[0].axhline(y=0.0387, color='red', linestyle='--', linewidth=1, label='Threshold')
axes[0].set_ylabel('IPS Score', fontsize=12)
axes[0].set_title('Training Set (GSE26378)', fontsize=14)
axes[0].legend()
axes[0].grid(True, alpha=0.3, axis='y')

# 验证集箱线图
val_ips = ips_scores[ips_scores['dataset'] == 'GSE26440']
data_normal_v = val_ips[val_ips['group'] == 'normal']['IPS_normalized']
data_sepsis_v = val_ips[val_ips['group'] == 'sepsis']['IPS_normalized']
bp2 = axes[1].boxplot([data_normal_v, data_sepsis_v], labels=['Normal', 'Sepsis'], patch_artist=True)
bp2['boxes'][0].set_facecolor('lightblue')
bp2['boxes'][1].set_facecolor('lightcoral')
axes[1].axhline(y=0.0387, color='red', linestyle='--', linewidth=1, label='Threshold')
axes[1].set_ylabel('IPS Score', fontsize=12)
axes[1].set_title('Validation Set (GSE26440)', fontsize=14)
axes[1].legend()
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_boxplot.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS箱线图已保存")

# ============= 保存带基因符号的IPS权重 =============
ips_weights.to_csv(f'{BASE_DIR}/results/ML/ips_weights_symbols.csv', index=False)

print("\n所有可视化图表已生成完成！")
