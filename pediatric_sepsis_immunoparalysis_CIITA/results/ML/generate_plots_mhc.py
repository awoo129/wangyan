#!/usr/bin/env python3
"""
生成可视化图表 - MHC基因版本
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
roc_data = pd.read_csv(f'{BASE_DIR}/results/ML/roc_curve_data_mhc.csv')
auc_data = pd.read_csv(f'{BASE_DIR}/results/ML/auc_values_mhc.csv')
ips_scores = pd.read_csv(f'{BASE_DIR}/results/ML/ips_scores_all_mhc.csv')
ips_weights = pd.read_csv(f'{BASE_DIR}/results/ML/ips_weights_mhc.csv')
feature_importance = pd.read_csv(f'{BASE_DIR}/results/ML/feature_importance_mhc.csv')
model_performance = pd.read_csv(f'{BASE_DIR}/results/ML/model_performance_mhc.csv')

# ============= 图1: ROC曲线 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 训练集ROC
train_roc = roc_data[roc_data['Dataset'] == 'Training']
auc_train = auc_data[auc_data['Dataset'] == 'Training']['AUC'].values[0]
axes[0].plot(train_roc['FPR'], train_roc['TPR'], 'b-', linewidth=2.5, label=f'ROC (AUC = {auc_train:.4f})')
axes[0].fill_between(train_roc['FPR'], train_roc['TPR'], alpha=0.2, color='blue')
axes[0].plot([0, 1], [0, 1], 'r--', linewidth=1.5, label='Random')
axes[0].set_xlabel('False Positive Rate', fontsize=12)
axes[0].set_ylabel('True Positive Rate', fontsize=12)
axes[0].set_title('ROC Curve - Training Set (GSE26378)', fontsize=14, fontweight='bold')
axes[0].legend(loc='lower right', fontsize=11)
axes[0].grid(True, alpha=0.3)

# 验证集ROC
val_roc = roc_data[roc_data['Dataset'] == 'Validation']
auc_val = auc_data[auc_data['Dataset'] == 'Validation']['AUC'].values[0]
axes[1].plot(val_roc['FPR'], val_roc['TPR'], 'g-', linewidth=2.5, label=f'ROC (AUC = {auc_val:.4f})')
axes[1].fill_between(val_roc['FPR'], val_roc['TPR'], alpha=0.2, color='green')
axes[1].plot([0, 1], [0, 1], 'r--', linewidth=1.5, label='Random')
axes[1].set_xlabel('False Positive Rate', fontsize=12)
axes[1].set_ylabel('True Positive Rate', fontsize=12)
axes[1].set_title('ROC Curve - Validation Set (GSE26440)', fontsize=14, fontweight='bold')
axes[1].legend(loc='lower right', fontsize=11)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/roc_curves_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("ROC曲线已保存: roc_curves_mhc.png")

# ============= 图2: IPS分布 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

threshold = 0.3636

# 训练集IPS分布
train_ips = ips_scores[ips_scores['dataset'] == 'GSE26378']
axes[0].hist(train_ips[train_ips['group'] == 'normal']['IPS_normalized'], bins=15, alpha=0.6, 
             color='steelblue', label=f'Normal (n={sum(train_ips["group"]=="normal")})', edgecolor='black')
axes[0].hist(train_ips[train_ips['group'] == 'sepsis']['IPS_normalized'], bins=15, alpha=0.6, 
             color='coral', label=f'Sepsis (n={sum(train_ips["group"]=="sepsis")})', edgecolor='black')
axes[0].axvline(x=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({threshold:.4f})')
axes[0].set_xlabel('IPS Score (Normalized)', fontsize=12)
axes[0].set_ylabel('Frequency', fontsize=12)
axes[0].set_title('IPS Distribution - Training Set (GSE26378)', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3)

# 验证集IPS分布
val_ips = ips_scores[ips_scores['dataset'] == 'GSE26440']
axes[1].hist(val_ips[val_ips['group'] == 'normal']['IPS_normalized'], bins=15, alpha=0.6, 
             color='steelblue', label=f'Normal (n={sum(val_ips["group"]=="normal")})', edgecolor='black')
axes[1].hist(val_ips[val_ips['group'] == 'sepsis']['IPS_normalized'], bins=15, alpha=0.6, 
             color='coral', label=f'Sepsis (n={sum(val_ips["group"]=="sepsis")})', edgecolor='black')
axes[1].axvline(x=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({threshold:.4f})')
axes[1].set_xlabel('IPS Score (Normalized)', fontsize=12)
axes[1].set_ylabel('Frequency', fontsize=12)
axes[1].set_title('IPS Distribution - Validation Set (GSE26440)', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_distribution_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS分布图已保存: ips_distribution_mhc.png")

# ============= 图3: 特征重要性 =============
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# LASSO系数
top_lasso = feature_importance.nlargest(5, 'LASSO_coef')
genes_lasso = top_lasso.index.tolist()
axes[0].barh(range(len(top_lasso)), top_lasso['LASSO_coef'], color='steelblue', edgecolor='black')
axes[0].set_yticks(range(len(top_lasso)))
axes[0].set_yticklabels(genes_lasso, fontsize=12)
axes[0].set_xlabel('LASSO Coefficient', fontsize=12)
axes[0].set_title('Feature Importance - LASSO', fontsize=14, fontweight='bold')
axes[0].invert_yaxis()
axes[0].grid(True, alpha=0.3, axis='x')

# 随机森林重要性
top_rf = feature_importance.nlargest(5, 'RF_importance')
genes_rf = top_rf.index.tolist()
axes[1].barh(range(len(top_rf)), top_rf['RF_importance'], color='forestgreen', edgecolor='black')
axes[1].set_yticks(range(len(top_rf)))
axes[1].set_yticklabels(genes_rf, fontsize=12)
axes[1].set_xlabel('Random Forest Importance', fontsize=12)
axes[1].set_title('Feature Importance - Random Forest', fontsize=14, fontweight='bold')
axes[1].invert_yaxis()
axes[1].grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/feature_importance_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("特征重要性图已保存: feature_importance_mhc.png")

# ============= 图4: IPS权重柱状图 =============
fig, ax = plt.subplots(figsize=(10, 6))
genes = ips_weights['Gene'].values
weights = ips_weights['Weight'].values
colors = ['steelblue' if w > 0 else 'coral' for w in weights]
bars = ax.barh(range(len(genes)), weights, color=colors, edgecolor='black')
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=12)
ax.set_xlabel('IPS Weight (Logistic Regression Coefficient)', fontsize=12)
ax.set_title('IPS Gene Weights', fontsize=14, fontweight='bold')
ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
ax.invert_yaxis()
ax.grid(True, alpha=0.3, axis='x')

# 添加数值标签
for bar, w in zip(bars, weights):
    width = bar.get_width()
    ax.annotate(f'{w:.3f}', xy=(width, bar.get_y() + bar.get_height()/2),
                xytext=(5 if width > 0 else -5, 0), textcoords="offset points",
                ha='left' if width > 0 else 'right', va='center', fontsize=10)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_weights_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS权重图已保存: ips_weights_mhc.png")

# ============= 图5: 模型性能对比 =============
fig, ax = plt.subplots(figsize=(10, 6))
metrics = ['Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1', 'AUC']
train_vals = [model_performance[model_performance['Model'].str.contains('Training')][m].values[0] for m in metrics]
val_vals = [model_performance[model_performance['Model'].str.contains('Validation')][m].values[0] for m in metrics]

x = np.arange(len(metrics))
width = 0.35
bars1 = ax.bar(x - width/2, train_vals, width, label='Training', color='steelblue', edgecolor='black')
bars2 = ax.bar(x + width/2, val_vals, width, label='Validation', color='forestgreen', edgecolor='black')

ax.set_ylabel('Score', fontsize=12)
ax.set_title('Model Performance Comparison', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=11)
ax.legend(fontsize=11)
ax.set_ylim([0, 1.15])
ax.grid(True, alpha=0.3, axis='y')

# 添加数值标签
for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/model_performance_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("模型性能对比图已保存: model_performance_mhc.png")

# ============= 图6: 箱线图比较 =============
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
threshold = 0.3636

# 训练集箱线图
train_ips = ips_scores[ips_scores['dataset'] == 'GSE26378']
data_normal = train_ips[train_ips['group'] == 'normal']['IPS_normalized']
data_sepsis = train_ips[train_ips['group'] == 'sepsis']['IPS_normalized']
bp1 = axes[0].boxplot([data_normal, data_sepsis], labels=['Normal', 'Sepsis'], patch_artist=True, widths=0.6)
bp1['boxes'][0].set_facecolor('lightblue')
bp1['boxes'][1].set_facecolor('lightcoral')
for box in bp1['boxes']:
    box.set_edgecolor('black')
    box.set_linewidth(1.5)
axes[0].axhline(y=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({threshold:.4f})')
axes[0].set_ylabel('IPS Score', fontsize=12)
axes[0].set_title('Training Set (GSE26378)', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3, axis='y')

# 验证集箱线图
val_ips = ips_scores[ips_scores['dataset'] == 'GSE26440']
data_normal_v = val_ips[val_ips['group'] == 'normal']['IPS_normalized']
data_sepsis_v = val_ips[val_ips['group'] == 'sepsis']['IPS_normalized']
bp2 = axes[1].boxplot([data_normal_v, data_sepsis_v], labels=['Normal', 'Sepsis'], patch_artist=True, widths=0.6)
bp2['boxes'][0].set_facecolor('lightblue')
bp2['boxes'][1].set_facecolor('lightcoral')
for box in bp2['boxes']:
    box.set_edgecolor('black')
    box.set_linewidth(1.5)
axes[1].axhline(y=threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({threshold:.4f})')
axes[1].set_ylabel('IPS Score', fontsize=12)
axes[1].set_title('Validation Set (GSE26440)', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/ips_boxplot_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("IPS箱线图已保存: ips_boxplot_mhc.png")

# ============= 图7: 综合结果汇总图 =============
fig = plt.figure(figsize=(16, 12))

# 1. ROC曲线
ax1 = fig.add_subplot(2, 3, 1)
ax1.plot(train_roc['FPR'], train_roc['TPR'], 'b-', linewidth=2, label=f'Train (AUC={auc_train:.4f})')
ax1.plot(val_roc['FPR'], val_roc['TPR'], 'g-', linewidth=2, label=f'Val (AUC={auc_val:.4f})')
ax1.plot([0, 1], [0, 1], 'r--', linewidth=1)
ax1.set_xlabel('FPR')
ax1.set_ylabel('TPR')
ax1.set_title('ROC Curves', fontweight='bold')
ax1.legend(loc='lower right', fontsize=9)
ax1.grid(True, alpha=0.3)

# 2. IPS分布-训练集
ax2 = fig.add_subplot(2, 3, 2)
ax2.hist(train_ips[train_ips['group']=='normal']['IPS_normalized'], bins=12, alpha=0.6, color='steelblue', label='Normal')
ax2.hist(train_ips[train_ips['group']=='sepsis']['IPS_normalized'], bins=12, alpha=0.6, color='coral', label='Sepsis')
ax2.axvline(x=threshold, color='red', linestyle='--', linewidth=2)
ax2.set_xlabel('IPS Score')
ax2.set_ylabel('Count')
ax2.set_title('Training Set IPS Distribution', fontweight='bold')
ax2.legend(fontsize=9)

# 3. IPS分布-验证集
ax3 = fig.add_subplot(2, 3, 3)
ax3.hist(val_ips[val_ips['group']=='normal']['IPS_normalized'], bins=12, alpha=0.6, color='steelblue', label='Normal')
ax3.hist(val_ips[val_ips['group']=='sepsis']['IPS_normalized'], bins=12, alpha=0.6, color='coral', label='Sepsis')
ax3.axvline(x=threshold, color='red', linestyle='--', linewidth=2)
ax3.set_xlabel('IPS Score')
ax3.set_ylabel('Count')
ax3.set_title('Validation Set IPS Distribution', fontweight='bold')
ax3.legend(fontsize=9)

# 4. 基因权重
ax4 = fig.add_subplot(2, 3, 4)
genes = ips_weights['Gene'].values
weights = ips_weights['Weight'].values
colors = ['steelblue' if w > 0 else 'coral' for w in weights]
ax4.barh(range(len(genes)), weights, color=colors)
ax4.set_yticks(range(len(genes)))
ax4.set_yticklabels(genes, fontsize=10)
ax4.axvline(x=0, color='black', linestyle='-')
ax4.set_xlabel('Weight')
ax4.set_title('IPS Gene Weights', fontweight='bold')
ax4.invert_yaxis()

# 5. 特征重要性
ax5 = fig.add_subplot(2, 3, 5)
top_features = feature_importance.nlargest(5, 'RF_importance')
ax5.barh(range(len(top_features)), top_features['RF_importance'], color='forestgreen')
ax5.set_yticks(range(len(top_features)))
ax5.set_yticklabels(top_features.index.tolist(), fontsize=10)
ax5.set_xlabel('Importance')
ax5.set_title('RF Feature Importance', fontweight='bold')
ax5.invert_yaxis()

# 6. 模型性能
ax6 = fig.add_subplot(2, 3, 6)
metrics_short = ['Acc', 'Sens', 'Spec', 'AUC']
train_vals_short = [train_vals[0], train_vals[1], train_vals[2], train_vals[5]]
val_vals_short = [val_vals[0], val_vals[1], val_vals[2], val_vals[5]]
x = np.arange(len(metrics_short))
ax6.bar(x - 0.2, train_vals_short, 0.4, label='Train', color='steelblue')
ax6.bar(x + 0.2, val_vals_short, 0.4, label='Val', color='forestgreen')
ax6.set_xticks(x)
ax6.set_xticklabels(metrics_short)
ax6.set_ylabel('Score')
ax6.set_title('Model Performance', fontweight='bold')
ax6.legend()
ax6.set_ylim([0, 1.15])

plt.tight_layout()
plt.savefig(f'{BASE_DIR}/results/ML/summary_figures_mhc.png', dpi=150, bbox_inches='tight')
plt.close()
print("综合汇总图已保存: summary_figures_mhc.png")

print("\n所有可视化图表已生成完成！")
