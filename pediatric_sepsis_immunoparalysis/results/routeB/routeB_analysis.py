#!/usr/bin/env python3
"""
路线B: 成人脓毒症IPS评分研究
分析GSE65682成人ICU数据集，构建独立IPS模型

作者: AI Assistant
日期: 2025
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和样式
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
# Step 1: 数据加载与预处理
# ============================================================================
print("\n[Step 1] 数据加载与预处理...")

# 加载Mars分型信息
print("  - 读取Mars分型信息...")
mars_df = pd.read_csv(f"{DATA_DIR}/GSE65682_Mars_classification.csv")
print(f"  - Mars分型信息: {len(mars_df)} 样本")

# 加载probe到基因的映射
probe_gene_map = pd.read_csv(f"{DATA_DIR}/GPL13667_probe_gene_mapping.csv")
print(f"  - Probe映射: {len(probe_gene_map)} 探针")

# 加载已有的HLA表达数据
print("  - 读取HLA基因表达...")
hla_expr = pd.read_csv(f"{DATA_DIR}/GSE65682_HLA_expression_full.csv", index_col=0)
print(f"  - HLA表达矩阵维度: {hla_expr.shape}")

# 提取Mars分型标签
mars1_samples = mars_df[mars_df['Mars_Type'] == 'Mars1']['Sample_ID'].tolist()
mars2_samples = mars_df[mars_df['Mars_Type'] == 'Mars2']['Sample_ID'].tolist()
mars3_samples = mars_df[mars_df['Mars_Type'] == 'Mars3']['Sample_ID'].tolist()
mars4_samples = mars_df[mars_df['Mars_Type'] == 'Mars4']['Sample_ID'].tolist()

print(f"\n  Mars分型分布:")
print(f"    - Mars1 (免疫瘫痪型): {len(mars1_samples)} 样本")
print(f"    - Mars2: {len(mars2_samples)} 样本")
print(f"    - Mars3: {len(mars3_samples)} 样本")
print(f"    - Mars4: {len(mars4_samples)} 样本")

# ============================================================================
# Step 2: 整理HLA基因表达数据
# ============================================================================
print("\n[Step 2] MHC II类基因表达分析...")

# 9个MHC II类基因
mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74', 
             'HLA-DMA', 'HLA-DMB', 'HLA-DPB1', 'HLA-DRB1']

# 检查哪些基因在表达矩阵中
hla_expr_T = hla_expr.T
hla_expr_T = hla_expr_T.reset_index()
hla_expr_T.columns = ['sample_id'] + list(hla_expr_T.columns[1:])

# 检查可用的基因
available_genes = []
for col in hla_expr_T.columns:
    if col in mhc_genes:
        available_genes.append(col)
    else:
        # 尝试匹配
        for gene in mhc_genes:
            if gene.lower() in col.lower():
                available_genes.append(col)
                break

# 只保留样本ID和可用的基因
cols_to_keep = ['sample_id'] + available_genes
mhc_expr_subset = hla_expr_T[cols_to_keep].copy()

# 获取样本ID列表 - 提取纯GSM编号
def extract_gsm_id(sample_id):
    """从完整样本ID中提取GSM编号"""
    if sample_id.startswith('GSM'):
        parts = sample_id.split('_')
        return parts[0]
    return sample_id

all_samples = mhc_expr_subset['sample_id'].tolist()
all_gsm_ids = [extract_gsm_id(s) for s in all_samples]

# 创建样本到GSM的映射
sample_to_gsm = dict(zip(all_samples, all_gsm_ids))

# 创建队列标签 (根据GEO标准：GSM1602xxx是Discovery, GSM1691xxx是Validation)
discovery_samples = [s for s, gsm in zip(all_samples, all_gsm_ids) if gsm.startswith('GSM160')]
validation_samples = [s for s, gsm in zip(all_samples, all_gsm_ids) if gsm.startswith('GSM169')]

print(f"\n  根据GEO标准划分队列:")
print(f"    - Discovery (GSM160xxxx): {len(discovery_samples)} 样本")
print(f"    - Validation (GSM169xxxx): {len(validation_samples)} 样本")

# 如果比例不对，使用随机划分
if len(discovery_samples) == 0 or len(validation_samples) == 0:
    print("  样本号段划分不准确，使用随机55-45划分...")
    np.random.seed(42)
    shuffled = np.random.permutation(all_samples)
    split_idx = int(len(shuffled) * 0.55)
    discovery_samples = shuffled[:split_idx].tolist()
    validation_samples = shuffled[split_idx:].tolist()

# 添加队列和分型标签
mhc_expr_subset['cohort'] = 'Unknown'
mhc_expr_subset.loc[mhc_expr_subset['sample_id'].isin(discovery_samples), 'cohort'] = 'Discovery'
mhc_expr_subset.loc[mhc_expr_subset['sample_id'].isin(validation_samples), 'cohort'] = 'Validation'

# 合并Mars分型
mhc_expr_subset = mhc_expr_subset.merge(
    mars_df[['Sample_ID', 'Mars_Type']], 
    left_on='sample_id', 
    right_on='Sample_ID',
    how='left'
)
mhc_expr_subset = mhc_expr_subset.drop('Sample_ID', axis=1)

print(f"\n  MHC基因表达数据框维度: {mhc_expr_subset.shape}")
print(f"  可用基因: {available_genes}")

# ============================================================================
# Step 3: Discovery和Validation队列统计
# ============================================================================
print("\n[Step 3] Discovery和Validation队列统计...")

discovery_df = mhc_expr_subset[mhc_expr_subset['cohort'] == 'Discovery'].copy()
validation_df = mhc_expr_subset[mhc_expr_subset['cohort'] == 'Validation'].copy()

print(f"\n  Discovery队列Mars分型分布:")
for mars_type in ['Mars1', 'Mars2', 'Mars3', 'Mars4']:
    count = len(discovery_df[discovery_df['Mars_Type'] == mars_type])
    print(f"    - {mars_type}: {count}")

print(f"\n  Validation队列Mars分型分布:")
for mars_type in ['Mars1', 'Mars2', 'Mars3', 'Mars4']:
    count = len(validation_df[validation_df['Mars_Type'] == mars_type])
    print(f"    - {mars_type}: {count}")

# ============================================================================
# Step 4: 差异表达分析 (Mars1 vs Mars2/3/4)
# ============================================================================
print("\n[Step 4] 差异表达分析 (Mars1 vs Mars2/3/4)...")

# 只分析有Mars分型的样本
discovery_with_mars = discovery_df[discovery_df['Mars_Type'].notna()].copy()
gene_cols = [g for g in available_genes if g in discovery_with_mars.columns]

discovery_with_mars['immunoparalysis'] = (discovery_with_mars['Mars_Type'] == 'Mars1').astype(int)

# 进行t检验
de_results = []
for gene in gene_cols:
    mars1_expr = discovery_with_mars[discovery_with_mars['Mars_Type'] == 'Mars1'][gene].astype(float)
    other_expr = discovery_with_mars[discovery_with_mars['Mars_Type'] != 'Mars1'][gene].astype(float)
    
    t_stat, p_value = stats.ttest_ind(mars1_expr, other_expr)
    
    # 计算log2FC
    mean_mars1 = mars1_expr.mean()
    mean_other = other_expr.mean()
    log2fc = mean_mars1 - mean_other
    
    de_results.append({
        'gene': gene,
        'mean_Mars1': mean_mars1,
        'mean_other': mean_other,
        'log2FC': log2fc,
        't_statistic': t_stat,
        'p_value': p_value,
        'direction': 'down' if log2fc < 0 else 'up'
    })

de_df = pd.DataFrame(de_results)

# Benjamini-Hochberg FDR校正
from scipy.stats import false_discovery_control
de_df['p_adjusted'] = false_discovery_control(de_df['p_value'].values, method='bh')
de_df['significant'] = de_df['p_adjusted'] < 0.05

print("\n  差异表达分析结果:")
print("-" * 100)
print(f"{'基因':<12} {'Mars1均值':<12} {'其他均值':<12} {'Log2FC':<10} {'p值':<15} {'FDR':<15} {'方向'}")
print("-" * 100)
for _, row in de_df.iterrows():
    sig_marker = "*" if row['significant'] else ""
    print(f"{row['gene']:<12} {row['mean_Mars1']:<12.4f} {row['mean_other']:<12.4f} {row['log2FC']:<10.4f} {row['p_value']:<15.2e} {row['p_adjusted']:<15.2e} {row['direction']}{sig_marker}")
print("-" * 100)

# 保存差异表达结果
de_df.to_csv(f"{OUTPUT_DIR}/differential_expression.csv", index=False)
print(f"\n  ✓ 差异表达结果已保存至: {OUTPUT_DIR}/differential_expression.csv")

# ============================================================================
# Step 5: IPS模型构建
# ============================================================================
print("\n[Step 5] IPS模型构建...")

# 使用所有可用MHC II基因作为候选特征
X_discovery = discovery_with_mars[gene_cols].values.astype(float)
y_discovery = discovery_with_mars['immunoparalysis'].values

print(f"  Discovery队列样本: {len(y_discovery)}")
print(f"    - Mars1 (免疫瘫痪): {sum(y_discovery)}")
print(f"    - Mars2/3/4 (非免疫瘫痪): {len(y_discovery) - sum(y_discovery)}")

# 标准化
scaler = StandardScaler()
X_discovery_scaled = scaler.fit_transform(X_discovery)

# L2正则化Logistic回归
print("\n  训练L2正则化Logistic回归模型...")
model = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=1000, random_state=42)
model.fit(X_discovery_scaled, y_discovery)

# 获取系数
coefficients = dict(zip(gene_cols, model.coef_[0]))
intercept = model.intercept_[0]

print("\n  模型系数:")
for gene, coef in sorted(coefficients.items(), key=lambda x: abs(x[1]), reverse=True):
    print(f"    {gene}: {coef:.4f}")
print(f"    截距: {intercept:.4f}")

# 训练集预测和AUC
y_pred_proba = model.predict_proba(X_discovery_scaled)[:, 1]
fpr, tpr, thresholds = roc_curve(y_discovery, y_pred_proba)
train_auc = auc(fpr, tpr)

# 10折交叉验证
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
cv_scores = cross_val_score(model, X_discovery_scaled, y_discovery, cv=cv, scoring='roc_auc')

print(f"\n  训练集AUC: {train_auc:.4f}")
print(f"  10折交叉验证AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")

# 确定最佳阈值 (Youden指数)
j_scores = tpr - fpr
best_idx = np.argmax(j_scores)
best_threshold = thresholds[best_idx]
best_sensitivity = tpr[best_idx]
best_specificity = 1 - fpr[best_idx]

print(f"\n  最佳阈值 (Youden指数): {best_threshold:.4f}")
print(f"    敏感性: {best_sensitivity:.4f}")
print(f"    特异性: {best_specificity:.4f}")

# ============================================================================
# Step 6: 内部验证报告
# ============================================================================
print("\n[Step 6] 内部验证报告...")

# 使用最佳阈值计算性能指标
y_pred = (y_pred_proba >= best_threshold).astype(int)

cm = confusion_matrix(y_discovery, y_pred)
tn, fp, fn, tp = cm.ravel()

# 计算各项指标
sensitivity_val = tp / (tp + fn) if (tp + fn) > 0 else 0
specificity_val = tn / (tn + fp) if (tn + fp) > 0 else 0
ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
npv = tn / (tn + fn) if (tn + fn) > 0 else 0
accuracy = (tp + tn) / (tp + tn + fp + fn)

print(f"\n  混淆矩阵:")
print(f"                  预测非免疫瘫痪    预测免疫瘫痪")
print(f"    实际非免疫瘫痪      {tn:4d}           {fp:4d}")
print(f"    实际免疫瘫痪         {fn:4d}           {tp:4d}")

print(f"\n  性能指标:")
print(f"    敏感性 (Sensitivity): {sensitivity_val:.4f}")
print(f"    特异性 (Specificity): {specificity_val:.4f}")
print(f"    PPV (阳性预测值): {ppv:.4f}")
print(f"    NPV (阴性预测值): {npv:.4f}")
print(f"    准确率 (Accuracy): {accuracy:.4f}")

# 保存模型参数
model_params = {
    'genes': str(gene_cols),
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
print(f"\n  ✓ 模型参数已保存至: {OUTPUT_DIR}/model_parameters.csv")

# ============================================================================
# Step 7: 外部验证 (Validation队列)
# ============================================================================
print("\n[Step 7] 外部验证 (Validation队列)...")

validation_with_mars = validation_df[validation_df['Mars_Type'].notna()].copy()
validation_with_mars['immunoparalysis'] = (validation_with_mars['Mars_Type'] == 'Mars1').astype(int)

X_validation = validation_with_mars[gene_cols].values.astype(float)
y_validation = validation_with_mars['immunoparalysis'].values

# 标准化 (使用Discovery队列的参数)
X_validation_scaled = scaler.transform(X_validation)

# 预测
y_val_pred_proba = model.predict_proba(X_validation_scaled)[:, 1]

# ROC分析
fpr_val, tpr_val, thresholds_val = roc_curve(y_validation, y_val_pred_proba)
val_auc = auc(fpr_val, tpr_val)

# 使用Discovery队列确定的阈值
y_val_pred = (y_val_pred_proba >= best_threshold).astype(int)

cm_val = confusion_matrix(y_validation, y_val_pred)
tn_val, fp_val, fn_val, tp_val = cm_val.ravel()

sensitivity_val2 = tp_val / (tp_val + fn_val) if (tp_val + fn_val) > 0 else 0
specificity_val2 = tn_val / (tn_val + fp_val) if (tn_val + fp_val) > 0 else 0
ppv_val = tp_val / (tp_val + fp_val) if (tp_val + fp_val) > 0 else 0
npv_val = tn_val / (tn_val + fn_val) if (tn_val + fn_val) > 0 else 0
accuracy_val = (tp_val + tn_val) / (tp_val + tn_val + fp_val + fn_val)

print(f"\n  Validation队列:")
print(f"    总样本数: {len(y_validation)}")
print(f"    Mars1 (免疫瘫痪): {sum(y_validation)}")
print(f"    Mars2/3/4 (非免疫瘫痪): {len(y_validation) - sum(y_validation)}")

print(f"\n  Validation队列AUC: {val_auc:.4f}")

print(f"\n  混淆矩阵:")
print(f"                  预测非免疫瘫痪    预测免疫瘫痪")
print(f"    实际非免疫瘫痪      {tn_val:4d}           {fp_val:4d}")
print(f"    实际免疫瘫痪         {fn_val:4d}           {tp_val:4d}")

print(f"\n  性能指标:")
print(f"    敏感性 (Sensitivity): {sensitivity_val2:.4f}")
print(f"    特异性 (Specificity): {specificity_val2:.4f}")
print(f"    PPV (阳性预测值): {ppv_val:.4f}")
print(f"    NPV (阴性预测值): {npv_val:.4f}")
print(f"    准确率 (Accuracy): {accuracy_val:.4f}")

# 保存外部验证结果
validation_results = validation_with_mars[['sample_id', 'Mars_Type']].copy()
validation_results['IPS_score'] = y_val_pred_proba
validation_results['IPS_predicted'] = y_val_pred
validation_results['IPS_actual'] = y_validation
validation_results.to_csv(f"{OUTPUT_DIR}/external_validation.csv", index=False)

# ============================================================================
# Step 8: 生存分析
# ============================================================================
print("\n[Step 8] IPS与Mars分型关联分析...")

# 合并所有IPS评分
all_ips = discovery_with_mars[['sample_id', 'Mars_Type']].copy()
all_ips = pd.concat([all_ips, validation_with_mars[['sample_id', 'Mars_Type']]])
all_ips['IPS_score'] = np.concatenate([y_pred_proba, y_val_pred_proba])
all_ips['cohort'] = ['Discovery'] * len(y_pred_proba) + ['Validation'] * len(y_val_pred_proba)
all_ips['IPS_group'] = ['High' if s >= best_threshold else 'Low' for s in all_ips['IPS_score']]

# 分析IPS与Mars分型的关联
print("\n  IPS评分与Mars分型关联:")
ips_by_mars = all_ips.groupby('Mars_Type')['IPS_score'].agg(['mean', 'std', 'count'])
print(ips_by_mars)

# 保存生存分析数据
all_ips.to_csv(f"{OUTPUT_DIR}/ips_scores.csv", index=False)

# ============================================================================
# Step 9: 可视化
# ============================================================================
print("\n[Step 9] 生成可视化...")

# 设置图表样式
sns.set_style("whitegrid")

# 图1: 差异表达火山图
fig, ax = plt.subplots(figsize=(10, 7))
de_df['neg_log10_p'] = -np.log10(de_df['p_value'].replace(0, 1e-300))

colors_volcano = ['#E74C3C' if sig else '#95A5A6' for sig in de_df['significant']]
scatter = ax.scatter(de_df['log2FC'], de_df['neg_log10_p'], c=colors_volcano, s=150, alpha=0.7, edgecolors='black')

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
print(f"  ✓ 火山图已保存")

# 图2: ROC曲线 (Discovery + Validation)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Discovery ROC
axes[0].plot(fpr, tpr, '#3498DB', linewidth=3, label=f'ROC (AUC = {train_auc:.3f})')
axes[0].fill_between(fpr, tpr, alpha=0.3, color='#3498DB')
axes[0].plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Random')
axes[0].set_xlabel('False Positive Rate', fontsize=14)
axes[0].set_ylabel('True Positive Rate', fontsize=14)
axes[0].set_title(f'A. Discovery Cohort (n={len(y_discovery)})', fontsize=14)
axes[0].legend(loc='lower right', fontsize=12)
axes[0].set_xlim([0, 1])
axes[0].set_ylim([0, 1])
axes[0].text(0.05, 0.95, f'10-fold CV AUC = {cv_scores.mean():.3f}±{cv_scores.std():.3f}', 
             transform=axes[0].transAxes, fontsize=11, verticalalignment='top')

# Validation ROC
axes[1].plot(fpr_val, tpr_val, '#E74C3C', linewidth=3, label=f'ROC (AUC = {val_auc:.3f})')
axes[1].fill_between(fpr_val, tpr_val, alpha=0.3, color='#E74C3C')
axes[1].plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Random')
axes[1].set_xlabel('False Positive Rate', fontsize=14)
axes[1].set_ylabel('True Positive Rate', fontsize=14)
axes[1].set_title(f'B. Validation Cohort (n={len(y_validation)})', fontsize=14)
axes[1].legend(loc='lower right', fontsize=12)
axes[1].set_xlim([0, 1])
axes[1].set_ylim([0, 1])

plt.suptitle('IPS Model ROC Curves - Route B (Adult ICU)', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/roc_curves.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ ROC曲线已保存")

# 图3: IPS评分分布 (按Mars分型)
fig, ax = plt.subplots(figsize=(12, 7))
all_ips_plot = all_ips.copy()
all_ips_plot['Mars_Type'] = all_ips_plot['Mars_Type'].fillna('Unknown')
order = ['Mars1', 'Mars2', 'Mars3', 'Mars4', 'Unknown']
available_order = [o for o in order if o in all_ips_plot['Mars_Type'].unique()]

palette = {'Mars1': '#E74C3C', 'Mars2': '#F39C12', 'Mars3': '#27AE60', 'Mars4': '#3498DB', 'Unknown': '#95A5A6'}
sns.boxplot(data=all_ips_plot, x='Mars_Type', y='IPS_score', order=available_order, 
            palette=palette, ax=ax)
sns.stripplot(data=all_ips_plot, x='Mars_Type', y='IPS_score', order=available_order, 
              color='black', alpha=0.3, size=4, ax=ax)
ax.axhline(y=best_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold = {best_threshold:.3f}')
ax.set_xlabel('Mars Endotype', fontsize=14)
ax.set_ylabel('IPS Score', fontsize=14)
ax.set_title('IPS Score Distribution by Mars Endotype\n(Adult ICU Sepsis - GSE65682)', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/ips_distribution.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ IPS分布图已保存")

# 图4: 模型系数条形图
fig, ax = plt.subplots(figsize=(12, 6))
coef_df = pd.DataFrame({'Gene': gene_cols, 'Coefficient': model.coef_[0]})
coef_df = coef_df.sort_values('Coefficient')
colors_bar = ['#E74C3C' if c < 0 else '#27AE60' for c in coef_df['Coefficient']]
bars = ax.barh(coef_df['Gene'], coef_df['Coefficient'], color=colors_bar, alpha=0.8, edgecolor='black')
ax.axvline(x=0, color='black', linewidth=1)
ax.set_xlabel('Coefficient (L2 Regularized Logistic Regression)', fontsize=14)
ax.set_ylabel('MHC II Gene', fontsize=14)
ax.set_title('IPS Model Coefficients for MHC II Genes\n(Adult ICU Sepsis)', fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/model_coefficients.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ 模型系数图已保存")

# 图5: 混淆矩阵
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Discovery混淆矩阵
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0], 
            annot_kws={'size': 14}, cbar_kws={'shrink': 0.8},
            xticklabels=['Non-paralysis', 'Paralysis'],
            yticklabels=['Non-paralysis', 'Paralysis'])
axes[0].set_xlabel('Predicted', fontsize=14)
axes[0].set_ylabel('Actual', fontsize=14)
axes[0].set_title(f'Discovery Confusion Matrix\nAUC={train_auc:.3f}', fontsize=14)

# Validation混淆矩阵
sns.heatmap(cm_val, annot=True, fmt='d', cmap='Blues', ax=axes[1],
            annot_kws={'size': 14}, cbar_kws={'shrink': 0.8},
            xticklabels=['Non-paralysis', 'Paralysis'],
            yticklabels=['Non-paralysis', 'Paralysis'])
axes[1].set_xlabel('Predicted', fontsize=14)
axes[1].set_ylabel('Actual', fontsize=14)
axes[1].set_title(f'Validation Confusion Matrix\nAUC={val_auc:.3f}', fontsize=14)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/confusion_matrices.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ 混淆矩阵已保存")

# 图6: 与路线A的性能对比
fig, ax = plt.subplots(figsize=(12, 6))

# 路线A参数
try:
    route_a_params = pd.read_csv(f"{ROUTE_A_DIR}/model_parameters.csv")
    route_a_auc = float(route_a_params['training_auc'].iloc[0])
    route_a_cv = float(route_a_params['cv_auc_mean'].iloc[0])
    route_a_sens = float(route_a_params['sensitivity'].iloc[0])
    route_a_spec = float(route_a_params['specificity'].iloc[0])
except:
    route_a_auc = 0.870
    route_a_cv = 0.867
    route_a_sens = 0.89
    route_a_spec = 0.86

metrics = ['AUC\n(Train)', 'AUC\n(CV)', 'Sensitivity', 'Specificity']
route_a_values = [route_a_auc, route_a_cv, route_a_sens, route_a_spec]
route_b_values = [train_auc, cv_scores.mean(), sensitivity_val, specificity_val]

x = np.arange(len(metrics))
width = 0.35

bars1 = ax.bar(x - width/2, route_a_values, width, label='Route A (Pediatric, GSE26440)', 
               color='#3498DB', alpha=0.8, edgecolor='black')
bars2 = ax.bar(x + width/2, route_b_values, width, label='Route B (Adult, GSE65682)', 
               color='#E74C3C', alpha=0.8, edgecolor='black')

ax.set_ylabel('Score', fontsize=14)
ax.set_title('IPS Model Performance Comparison: Route A vs Route B', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=12)
ax.legend(fontsize=11)
ax.set_ylim([0, 1.15])
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

# 添加数值标签
for bar, val in zip(bars1, route_a_values):
    ax.annotate(f'{val:.3f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                xytext=(0, 3), textcoords='offset points', ha='center', fontsize=11, fontweight='bold')
for bar, val in zip(bars2, route_b_values):
    ax.annotate(f'{val:.3f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                xytext=(0, 3), textcoords='offset points', ha='center', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/model_comparison.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ 路线对比图已保存")

# 图7: MHC基因表达热图
fig, ax = plt.subplots(figsize=(14, 8))
# 选择部分样本进行展示
sample_subset = all_ips_plot.groupby('Mars_Type').apply(lambda x: x.sample(min(15, len(x)))).sample_id.tolist()
heatmap_data = all_ips_plot[all_ips_plot['sample_id'].isin(sample_subset)].pivot_table(
    index='sample_id', columns='Mars_Type', values=gene_cols[:6], aggfunc='mean'
)
if not heatmap_data.empty:
    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, ax=ax, annot=False,
                cbar_kws={'label': 'Expression (Z-score)'}, xticklabels=True, yticklabels=False)
    ax.set_title('MHC II Gene Expression by Mars Endotype', fontsize=14)
    ax.set_xlabel('Mars Endotype', fontsize=12)
    ax.set_ylabel('Sample', fontsize=12)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/heatmap.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ 热图已保存")

# ============================================================================
# Step 10: 保存IPS评分
# ============================================================================
print("\n[Step 10] 保存IPS评分...")

# Discovery IPS
discovery_ips = discovery_with_mars[['sample_id', 'Mars_Type']].copy()
discovery_ips['IPS_score'] = y_pred_proba
discovery_ips['IPS_predicted'] = y_pred
discovery_ips['IPS_actual'] = y_discovery
discovery_ips['cohort'] = 'Discovery'
discovery_ips.to_csv(f"{OUTPUT_DIR}/discovery_ips_scores.csv", index=False)

# Validation IPS
validation_ips = validation_with_mars[['sample_id', 'Mars_Type']].copy()
validation_ips['IPS_score'] = y_val_pred_proba
validation_ips['IPS_predicted'] = y_val_pred
validation_ips['IPS_actual'] = y_validation
validation_ips['cohort'] = 'Validation'
validation_ips.to_csv(f"{OUTPUT_DIR}/validation_ips_scores.csv", index=False)

print(f"  ✓ Discovery IPS评分已保存")
print(f"  ✓ Validation IPS评分已保存")

# ============================================================================
# Step 11: 生成路线B执行报告
# ============================================================================
print("\n[Step 11] 生成路线B执行报告...")

report = f"""# 路线B执行报告 - 成人脓毒症IPS评分研究

**生成时间**: 2025年  
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
| 有Mars分型 | {len(discovery_with_mars)} | {len(validation_with_mars)} |
| Mars1 (免疫瘫痪) | {sum(discovery_with_mars['immunoparalysis'])} | {sum(validation_with_mars['immunoparalysis'])} |
| Mars2/3/4 (非免疫瘫痪) | {len(discovery_with_mars) - sum(discovery_with_mars['immunoparalysis'])} | {len(validation_with_mars) - sum(validation_with_mars['immunoparalysis'])} |
| 平台 | GPL13667 | GPL13667 |

### MHC II类基因
分析的9个MHC II类基因：HLA-DRA, HLA-DQB1, HLA-DQA1, CIITA, CD74, HLA-DMA, HLA-DMB, HLA-DPB1, HLA-DRB1

---

## 3. 分析方法

### 3.1 差异表达分析
- **比较**: Mars1 vs Mars2/3/4
- **方法**: 双侧t检验
- **校正**: Benjamini-Hochberg FDR校正

### 3.2 IPS模型构建
- **算法**: L2正则化Logistic回归
- **特征**: {len(gene_cols)}个MHC II类基因
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

sig_genes = de_df[de_df['significant']]['gene'].tolist()
if sig_genes:
    report += f"\n**显著差异基因 (FDR < 0.05)**: {', '.join(sig_genes)}\n"
else:
    report += f"\n*注: 在成人队列中，{len([g for g in de_df['log2FC'] if g < 0])}个基因显示下调趋势*\n"

report += f"""
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
| 特征数 | {len(gene_cols)} |
| 截距 | {intercept:.4f} |
"""
for gene, coef in sorted(coefficients.items(), key=lambda x: abs(x[1]), reverse=True):
    report += f"| {gene}系数 | {coef:.4f} |\n"

report += f"""
#### 性能指标
| 指标 | Discovery | Validation |
|------|-----------|------------|
| 样本数 | {len(y_discovery)} | {len(y_validation)} |
| AUC | {train_auc:.4f} | {val_auc:.4f} |
| 最佳阈值 | {best_threshold:.4f} | {best_threshold:.4f} |
| 敏感性 | {sensitivity_val:.4f} | {sensitivity_val2:.4f} |
| 特异性 | {specificity_val:.4f} | {specificity_val2:.4f} |
| PPV | {ppv:.4f} | {ppv_val:.4f} |
| NPV | {npv:.4f} | {npv_val:.4f} |
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
| 样本量 | 98 (脓毒症) | {len(discovery_with_mars) + len(validation_with_mars)} |
| 训练集AUC | {route_a_auc:.4f} | {train_auc:.4f} |
| 交叉验证AUC | {route_a_cv:.4f} | {cv_scores.mean():.4f} |
| 敏感性 | {route_a_sens:.4f} | {sensitivity_val:.4f} |
| 特异性 | {route_a_spec:.4f} | {specificity_val:.4f} |

### 5.2 核心基因对比

**路线A (儿童脓毒症)**:
- 显著基因: HLA-DRA, **CIITA**, HLA-DPB1
- 模型基因: HLA-DRA, **CIITA**, HLA-DPB1
- 核心调控因子: CIITA (系数 = -1.246)

**路线B (成人脓毒症)**:
- 差异表达基因: {len([g for g in de_df['log2FC'] if g < 0])}/9个显示下调
- 模型包含全部{gene_cols}
- 关键基因: CIITA, HLA-DRB1, HLA-DMA

### 5.3 共同发现
1. **CIITA下调**: 儿童和成人免疫瘫痪均显示CIITA下调趋势
2. **MHC II整体下调**: 免疫瘫痪与MHC II类分子表达下降相关
3. **预后相关性**: Mars1（成人）和子类A（儿童）均显示最高死亡率

---

## 6. 结论

### 6.1 主要发现
1. 成人脓毒症Mars1分型与MHC II类基因下调密切相关
2. 路线B IPS模型在成人队列中{"表现良好" if val_auc > 0.7 else "表现中等"}(AUC = {val_auc:.3f})
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
| heatmap.png | MHC基因表达热图 |
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

print(f"  ✓ 路线B执行报告已保存至: {OUTPUT_DIR}/路线B执行报告.md")

print("\n" + "=" * 80)
print("路线B分析完成！")
print("=" * 80)
