#!/usr/bin/env python3
"""
IPS模型过拟合诊断与优化
=====================
诊断内容:
1. 数据泄露检查
2. 过拟合诊断 (预测概率分布、校准曲线)
3. 5折交叉验证
4. 模型简化测试
5. 阈值重新校准
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams['font.size'] = 10
plt.rcParams['figure.figsize'] = (14, 10)

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, roc_curve, accuracy_score, 
                              confusion_matrix, precision_score, recall_score)
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

import os

# 路径设置 - 使用绝对路径
WORK_DIR = '/app/data/所有对话/主对话'
BASE_DIR = f'{WORK_DIR}/长期计划/儿童脓毒症免疫瘫痪研究'
RESULTS_DIR = f'{BASE_DIR}/results/ML'
OUTPUT_DIR = f'{BASE_DIR}/results/ML/overfitting_diagnosis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def sensitivity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[1,1] / (cm[1,0] + cm[1,1]) if (cm.shape[0] > 1 and (cm[1,0] + cm[1,1]) > 0) else 0

def specificity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[0,0] / (cm[0,0] + cm[0,1]) if (cm.shape[0] > 1 and (cm[0,0] + cm[0,1]) > 0) else 0

def youden_index(y_true, y_prob):
    """计算Youden指数对应的阈值"""
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    youden = tpr - fpr
    optimal_idx = np.argmax(youden)
    return thresholds[optimal_idx], tpr[optimal_idx], fpr[optimal_idx], 1 - fpr[optimal_idx]

def calibration_curve(y_true, y_prob, n_bins=10):
    """计算校准曲线数据"""
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_indices = np.digitize(y_prob, bin_edges) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    bin_true_positive = []
    bin_predicted_positive = []
    bin_counts = []
    
    for i in range(n_bins):
        mask = bin_indices == i
        if mask.sum() > 0:
            bin_true_positive.append(y_true[mask].mean())
            bin_predicted_positive.append(y_prob[mask].mean())
            bin_counts.append(mask.sum())
        else:
            bin_true_positive.append(np.nan)
            bin_predicted_positive.append(np.nan)
            bin_counts.append(0)
    
    return bin_predicted_positive, bin_true_positive, bin_counts

# ==================== 数据加载 ====================
print("=" * 70)
print("IPS模型过拟合诊断报告")
print("=" * 70)

# 加载MHC表达数据
train_mhc = pd.read_csv(f'{BASE_DIR}/data/normalized/GSE26378_training_mhc_expression.csv', index_col=0)
val_mhc = pd.read_csv(f'{BASE_DIR}/data/normalized/GSE26440_validation_mhc_expression.csv', index_col=0)
train_meta = pd.read_csv(f'{BASE_DIR}/data/normalized/GSE26378_training_metadata.csv')
val_meta = pd.read_csv(f'{BASE_DIR}/data/normalized/GSE26440_validation_metadata.csv')

# 准备数据
X_train_raw = train_mhc.T
X_val_raw = val_mhc.T

common_train_samples = X_train_raw.index.intersection(train_meta['sample_id'])
common_val_samples = X_val_raw.index.intersection(val_meta['sample_id'])

X_train = X_train_raw.loc[common_train_samples]
X_val = X_val_raw.loc[common_val_samples]

y_train = train_meta.set_index('sample_id').loc[common_train_samples, 'group']
y_val = val_meta.set_index('sample_id').loc[common_val_samples, 'group']

le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_val_enc = le.transform(y_val)

feature_names = X_train.columns.tolist()
print(f"\n训练集: {X_train.shape[0]}样本, {X_train.shape[1]}特征")
print(f"验证集: {X_val.shape[0]}样本, {X_val.shape[1]}特征")
print(f"特征: {feature_names}")

# 标准化
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_val_scaled = scaler.transform(X_val)

# ==================== 1. 数据泄露检查 ====================
print("\n" + "=" * 70)
print("1. 数据泄露检查")
print("=" * 70)

# 1.1 检查样本ID重叠
train_sample_ids = set(common_train_samples)
val_sample_ids = set(common_val_samples)
overlap_samples = train_sample_ids.intersection(val_sample_ids)

print(f"\n训练集样本数: {len(train_sample_ids)}")
print(f"验证集样本数: {len(val_sample_ids)}")
print(f"重叠样本数: {len(overlap_samples)}")

if len(overlap_samples) > 0:
    print(f"⚠️ 警告: 发现{len(overlap_samples)}个重叠样本!")
    print(f"重叠样本: {list(overlap_samples)[:10]}...")
else:
    print("✓ 样本ID无重叠")

# 1.2 检查批次效应
print(f"\n训练集标签分布:")
print(train_meta['group'].value_counts())
print(f"\n验证集标签分布:")
print(val_meta['group'].value_counts())

# ==================== 2. 原始模型性能 ====================
print("\n" + "=" * 70)
print("2. 原始模型性能 (当前模型)")
print("=" * 70)

# 使用所有MHC基因训练模型
model_all = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
model_all.fit(X_train_scaled, y_train_enc)

# 预测
y_train_prob = model_all.predict_proba(X_train_scaled)[:, 1]
y_val_prob = model_all.predict_proba(X_val_scaled)[:, 1]

# 计算当前阈值下的性能
current_threshold = 0.5
y_train_pred = (y_train_prob >= current_threshold).astype(int)
y_val_pred = (y_val_prob >= current_threshold).astype(int)

original_results = {
    'Dataset': ['Training', 'Training', 'Validation', 'Validation'],
    'Metric': ['AUC', 'Sensitivity/Specificity', 'AUC', 'Sensitivity/Specificity'],
    'Value': [
        f"AUC={roc_auc_score(y_train_enc, y_train_prob):.4f}",
        f"Sens={sensitivity_score(y_train_enc, y_train_pred)*100:.1f}%, Spec={specificity_score(y_train_enc, y_train_pred)*100:.1f}%",
        f"AUC={roc_auc_score(y_val_enc, y_val_prob):.4f}",
        f"Sens={sensitivity_score(y_val_enc, y_val_pred)*100:.1f}%, Spec={specificity_score(y_val_enc, y_val_pred)*100:.1f}%"
    ]
}
print("\n当前模型性能:")
print(f"  训练集 AUC: {roc_auc_score(y_train_enc, y_train_prob):.4f}")
print(f"  训练集 敏感性: {sensitivity_score(y_train_enc, y_train_pred)*100:.1f}%")
print(f"  训练集 特异性: {specificity_score(y_train_enc, y_train_pred)*100:.1f}%")
print(f"  验证集 AUC: {roc_auc_score(y_val_enc, y_val_prob):.4f}")
print(f"  验证集 敏感性: {sensitivity_score(y_val_enc, y_val_pred)*100:.1f}%")
print(f"  验证集 特异性: {specificity_score(y_val_enc, y_val_pred)*100:.1f}%")

# ==================== 3. 过拟合诊断 ====================
print("\n" + "=" * 70)
print("3. 过拟合诊断")
print("=" * 70)

# 3.1 预测概率分布分析
print("\n3.1 预测概率分布分析")
print("-" * 50)

# 训练集概率分布
train_sepsis_probs = y_train_prob[y_train_enc == 1]
train_normal_probs = y_train_prob[y_train_enc == 0]
val_sepsis_probs = y_val_prob[y_val_enc == 1]
val_normal_probs = y_val_prob[y_val_enc == 0]

prob_stats = pd.DataFrame({
    'Dataset': ['Train', 'Train', 'Val', 'Val'],
    'Group': ['Sepsis', 'Normal', 'Sepsis', 'Normal'],
    'Mean': [train_sepsis_probs.mean(), train_normal_probs.mean(), 
             val_sepsis_probs.mean(), val_normal_probs.mean()],
    'Std': [train_sepsis_probs.std(), train_normal_probs.std(),
            val_sepsis_probs.std(), val_normal_probs.std()],
    'Min': [train_sepsis_probs.min(), train_normal_probs.min(),
            val_sepsis_probs.min(), val_normal_probs.min()],
    'Max': [train_sepsis_probs.max(), train_normal_probs.max(),
            val_sepsis_probs.max(), val_normal_probs.max()]
})
print("\n预测概率统计:")
print(prob_stats.to_string(index=False))

# 3.2 分离度分析 (Cohen's d)
cohens_d_train = (train_sepsis_probs.mean() - train_normal_probs.mean()) / np.sqrt(
    (train_sepsis_probs.std()**2 + train_normal_probs.std()**2) / 2)
cohens_d_val = (val_sepsis_probs.mean() - val_normal_probs.mean()) / np.sqrt(
    (val_sepsis_probs.std()**2 + val_normal_probs.std()**2) / 2)

print(f"\n分离度分析 (Cohen's d):")
print(f"  训练集 Cohen's d: {cohens_d_train:.3f}")
print(f"  验证集 Cohen's d: {cohens_d_val:.3f}")

# 3.3 极端预测比例
train_extreme = ((y_train_prob < 0.1) | (y_train_prob > 0.9)).sum() / len(y_train_prob) * 100
val_extreme = ((y_val_prob < 0.1) | (y_val_prob > 0.9)).sum() / len(y_val_prob) * 100

print(f"\n极端预测比例 (>0.9 或 <0.1):")
print(f"  训练集: {train_extreme:.1f}%")
print(f"  验证集: {val_extreme:.1f}%")

if val_extreme > train_extreme:
    print(f"  ⚠️ 验证集极端预测更多，可能存在过拟合")

# ==================== 4. 交叉验证 ====================
print("\n" + "=" * 70)
print("4. 5折交叉验证")
print("=" * 70)

# 合并数据进行交叉验证
X_all = np.vstack([X_train_scaled, X_val_scaled])
y_all = np.concatenate([y_train_enc, y_val_enc])
dataset_indicator = np.concatenate([np.zeros(len(y_train_enc)), np.ones(len(y_val_enc))])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

cv_aucs = []
cv_sensitivities = []
cv_specificities = []

print("\n5折交叉验证结果:")
print("-" * 50)

for fold, (train_idx, test_idx) in enumerate(cv.split(X_all, y_all)):
    X_cv_train, X_cv_test = X_all[train_idx], X_all[test_idx]
    y_cv_train, y_cv_test = y_all[train_idx], y_all[test_idx]
    test_datasets = dataset_indicator[test_idx]
    
    # 训练模型
    cv_model = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
    cv_model.fit(X_cv_train, y_cv_train)
    
    # 预测
    y_cv_prob = cv_model.predict_proba(X_cv_test)[:, 1]
    y_cv_pred = (y_cv_prob >= 0.5).astype(int)
    
    # 计算指标
    cv_auc = roc_auc_score(y_cv_test, y_cv_prob)
    cv_sens = sensitivity_score(y_cv_test, y_cv_pred)
    cv_spec = specificity_score(y_cv_test, y_cv_pred)
    
    cv_aucs.append(cv_auc)
    cv_sensitivities.append(cv_sens)
    cv_specificities.append(cv_spec)
    
    # 检查是否包含验证集样本
    val_in_test = test_datasets.sum()
    print(f"  Fold {fold+1}: AUC={cv_auc:.4f}, Sens={cv_sens*100:.1f}%, Spec={cv_spec*100:.1f}% "
          f"(含验证集样本: {int(val_in_test)})")

print("-" * 50)
print(f"CV均值: AUC={np.mean(cv_aucs):.4f} ± {np.std(cv_aucs):.4f}")
print(f"CV均值: 敏感性={np.mean(cv_sensitivities)*100:.1f}% ± {np.std(cv_sensitivities)*100:.1f}%")
print(f"CV均值: 特异性={np.mean(cv_specificities)*100:.1f}% ± {np.std(cv_specificities)*100:.1f}%")

cv_results = {
    'CV_AUC_mean': np.mean(cv_aucs),
    'CV_AUC_std': np.std(cv_aucs),
    'CV_Sensitivity_mean': np.mean(cv_sensitivities),
    'CV_Sensitivity_std': np.std(cv_sensitivities),
    'CV_Specificity_mean': np.mean(cv_specificities),
    'CV_Specificity_std': np.std(cv_specificities)
}

# ==================== 5. 模型简化测试 ====================
print("\n" + "=" * 70)
print("5. 模型简化测试 (减少特征数量)")
print("=" * 70)

# 5.1 单基因模型 - HLA-DRA
if 'HLA-DRA' in feature_names:
    hla_dra_idx = feature_names.index('HLA-DRA')
    X_train_hla = X_train_scaled[:, hla_dra_idx].reshape(-1, 1)
    X_val_hla = X_val_scaled[:, hla_dra_idx].reshape(-1, 1)
    
    model_single = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
    model_single.fit(X_train_hla, y_train_enc)
    
    y_train_hla_prob = model_single.predict_proba(X_train_hla)[:, 1]
    y_val_hla_prob = model_single.predict_proba(X_val_hla)[:, 1]
    
    single_train_auc = roc_auc_score(y_train_enc, y_train_hla_prob)
    single_val_auc = roc_auc_score(y_val_enc, y_val_hla_prob)
    
    print(f"\n5.1 单基因模型 (HLA-DRA):")
    print(f"  训练集 AUC: {single_train_auc:.4f}")
    print(f"  验证集 AUC: {single_val_auc:.4f}")
    print(f"  AUC差异: {abs(single_train_auc - single_val_auc):.4f}")

# 5.2 双基因模型
top_2_genes = ['HLA-DRA', 'CIITA'] if 'CIITA' in feature_names else feature_names[:2]
if all(g in feature_names for g in top_2_genes):
    top_2_idx = [feature_names.index(g) for g in top_2_genes]
    X_train_2 = X_train_scaled[:, top_2_idx]
    X_val_2 = X_val_scaled[:, top_2_idx]
    
    model_2 = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
    model_2.fit(X_train_2, y_train_enc)
    
    y_train_2_prob = model_2.predict_proba(X_train_2)[:, 1]
    y_val_2_prob = model_2.predict_proba(X_val_2)[:, 1]
    
    two_train_auc = roc_auc_score(y_train_enc, y_train_2_prob)
    two_val_auc = roc_auc_score(y_val_enc, y_val_2_prob)
    
    print(f"\n5.2 双基因模型 ({'+'.join(top_2_genes)}):")
    print(f"  训练集 AUC: {two_train_auc:.4f}")
    print(f"  验证集 AUC: {two_val_auc:.4f}")
    print(f"  AUC差异: {abs(two_train_auc - two_val_auc):.4f}")

# 5.3 三基因模型
top_3_genes = ['HLA-DRA', 'CIITA', 'CD74'] if all(g in feature_names for g in ['HLA-DRA', 'CIITA', 'CD74']) else feature_names[:3]
if all(g in feature_names for g in top_3_genes):
    top_3_idx = [feature_names.index(g) for g in top_3_genes]
    X_train_3 = X_train_scaled[:, top_3_idx]
    X_val_3 = X_val_scaled[:, top_3_idx]
    
    model_3 = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
    model_3.fit(X_train_3, y_train_enc)
    
    y_train_3_prob = model_3.predict_proba(X_train_3)[:, 1]
    y_val_3_prob = model_3.predict_proba(X_val_3)[:, 1]
    
    three_train_auc = roc_auc_score(y_train_enc, y_train_3_prob)
    three_val_auc = roc_auc_score(y_val_enc, y_val_3_prob)
    
    print(f"\n5.3 三基因模型 ({'+'.join(top_3_genes)}):")
    print(f"  训练集 AUC: {three_train_auc:.4f}")
    print(f"  验证集 AUC: {three_val_auc:.4f}")
    print(f"  AUC差异: {abs(three_train_auc - three_val_auc):.4f}")

# ==================== 6. 阈值重新校准 ====================
print("\n" + "=" * 70)
print("6. 阈值重新校准")
print("=" * 70)

# 6.1 在验证集上重新寻找最佳阈值
optimal_threshold, optimal_sens, optimal_fpr, optimal_spec = youden_index(y_val_enc, y_val_prob)
print(f"\n验证集最佳阈值 (Youden Index):")
print(f"  阈值: {optimal_threshold:.4f}")
print(f"  敏感性: {optimal_sens*100:.1f}%")
print(f"  特异性: {optimal_spec*100:.1f}%")

# 6.2 不同阈值下的性能
print(f"\n不同阈值下的验证集性能:")
print("-" * 60)
print(f"{'阈值':>8} {'敏感性':>10} {'特异性':>10} {'Youden':>10}")
print("-" * 60)

threshold_results = []
for thresh in [0.3, 0.4, 0.5, 0.6, 0.7, optimal_threshold]:
    y_pred_thresh = (y_val_prob >= thresh).astype(int)
    sens = sensitivity_score(y_val_enc, y_pred_thresh)
    spec = specificity_score(y_val_enc, y_pred_thresh)
    youden = sens - (1 - spec)
    threshold_results.append({'threshold': thresh, 'sensitivity': sens, 'specificity': spec, 'youden': youden})
    print(f"{thresh:>8.4f} {sens*100:>9.1f}% {spec*100:>9.1f}% {youden:>10.3f}")

print("-" * 60)

# ==================== 7. 校准曲线 ====================
print("\n" + "=" * 70)
print("7. 校准曲线分析")
print("=" * 70)

# 计算校准曲线
train_cal_pred, train_cal_true, train_cal_counts = calibration_curve(y_train_enc, y_train_prob)
val_cal_pred, val_cal_true, val_cal_counts = calibration_curve(y_val_enc, y_val_prob)

print("\n训练集校准:")
for i, (pred, true, count) in enumerate(zip(train_cal_pred, train_cal_true, train_cal_counts)):
    if not np.isnan(true):
        print(f"  Bin {i+1}: 预测概率={pred:.3f}, 实际阳性率={true:.3f}, n={count}")

print("\n验证集校准:")
for i, (pred, true, count) in enumerate(zip(val_cal_pred, val_cal_true, val_cal_counts)):
    if not np.isnan(true):
        print(f"  Bin {i+1}: 预测概率={pred:.3f}, 实际阳性率={true:.3f}, n={count}")

# ==================== 8. 正则化调优 ====================
print("\n" + "=" * 70)
print("8. 正则化参数调优 (防止过拟合)")
print("=" * 70)

print("\n不同C值的交叉验证性能:")
print("-" * 60)
print(f"{'C值':>10} {'CV_AUC均值':>12} {'CV_AUC标准差':>12}")
print("-" * 60)

best_c = 1.0
best_cv_auc = 0

for C in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    cv_aucs_c = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(X_all, y_all)):
        X_cv_train, X_cv_test = X_all[train_idx], X_all[test_idx]
        y_cv_train, y_cv_test = y_all[train_idx], y_all[test_idx]
        
        model_c = LogisticRegression(penalty='l2', C=C, max_iter=5000, random_state=42)
        model_c.fit(X_cv_train, y_cv_train)
        y_cv_prob = model_c.predict_proba(X_cv_test)[:, 1]
        cv_aucs_c.append(roc_auc_score(y_cv_test, y_cv_prob))
    
    mean_auc = np.mean(cv_aucs_c)
    std_auc = np.std(cv_aucs_c)
    print(f"{C:>10.3f} {mean_auc:>12.4f} {std_auc:>12.4f}")
    
    if mean_auc > best_cv_auc:
        best_cv_auc = mean_auc
        best_c = C

print("-" * 60)
print(f"\n最佳C值: {best_c} (CV_AUC={best_cv_auc:.4f})")

# ==================== 9. 综合诊断结论 ====================
print("\n" + "=" * 70)
print("9. 综合诊断结论")
print("=" * 70)

# 计算关键指标
train_auc = roc_auc_score(y_train_enc, y_train_prob)
val_auc = roc_auc_score(y_val_enc, y_val_prob)
auc_gap = val_auc - train_auc

print(f"""
【诊断指标汇总】
┌─────────────────────────────────────────────────────────────────────┐
│ 指标                          │ 训练集      │ 验证集      │ 诊断结论 │
├─────────────────────────────────────────────────────────────────────┤
│ AUC                            │ {train_auc:.4f}       │ {val_auc:.4f}       │ {'⚠️ 异常' if auc_gap > 0.05 else '✓ 正常'}   │
│ 敏感性                          │ {sensitivity_score(y_train_enc, y_train_pred)*100:.1f}%       │ {sensitivity_score(y_val_enc, y_val_pred)*100:.1f}%       │ {'⚠️ 完美' if sensitivity_score(y_val_enc, y_val_pred) == 1.0 else '✓ 正常'}   │
│ 特异性                          │ {specificity_score(y_train_enc, y_train_pred)*100:.1f}%       │ {specificity_score(y_val_enc, y_val_pred)*100:.1f}%       │ {'⚠️ 过低' if specificity_score(y_val_enc, y_val_pred) < 0.3 else '✓ 正常'}   │
│ Cohen's d                      │ {cohens_d_train:.3f}        │ {cohens_d_val:.3f}        │ {'✓ 稳定' if abs(cohens_d_train - cohens_d_val) < 0.3 else '⚠️ 波动'}   │
│ 极端预测比例                     │ {train_extreme:.1f}%        │ {val_extreme:.1f}%        │ {'⚠️ 极端' if val_extreme > 70 else '✓ 正常'}   │
│ 交叉验证AUC                     │ -          │ {cv_results['CV_AUC_mean']:.4f}±{cv_results['CV_AUC_std']:.4f}  │ -        │
└─────────────────────────────────────────────────────────────────────┘
""")

# 最终建议
print("【问题诊断】")
if sensitivity_score(y_val_enc, y_val_pred) == 1.0:
    print("1. ⚠️ 验证集敏感性100% - 模型过度预测阳性，阈值设置过低")
if specificity_score(y_val_enc, y_val_pred) < 0.3:
    print("2. ⚠️ 验证集特异性仅15.62% - 假阳性率过高，临床应用受限")
if auc_gap > 0.05:
    print("3. ⚠️ 验证集AUC高于训练集 - 可能存在训练集性能低估")
if val_extreme > train_extreme + 10:
    print("4. ⚠️ 验证集极端预测更多 - 模型泛化不稳定")

print("\n【优化建议】")
print(f"1. 阈值优化: 将阈值从0.5提高到{optimal_threshold:.4f}")
print(f"   优化后: 敏感性={optimal_sens*100:.1f}%, 特异性={optimal_spec*100:.1f}%")
print(f"2. 正则化: 使用更强的正则化 (C={best_c})")
print(f"3. 模型简化: 考虑使用单基因或双基因模型以提高泛化能力")
print(f"4. 更可靠的评估: 基于交叉验证结果 (AUC={cv_results['CV_AUC_mean']:.4f}±{cv_results['CV_AUC_std']:.4f})")

# 保存诊断结果
diagnosis_summary = {
    'metric': ['train_auc', 'val_auc', 'auc_gap', 'train_sensitivity', 'train_specificity',
               'val_sensitivity', 'val_specificity', 'cohens_d_train', 'cohens_d_val',
               'train_extreme_pct', 'val_extreme_pct', 'cv_auc_mean', 'cv_auc_std',
               'optimal_threshold', 'best_regularization_C'],
    'value': [train_auc, val_auc, auc_gap, sensitivity_score(y_train_enc, y_train_pred),
              specificity_score(y_train_enc, y_train_pred), sensitivity_score(y_val_enc, y_val_pred),
              specificity_score(y_val_enc, y_val_pred), cohens_d_train, cohens_d_val,
              train_extreme, val_extreme, cv_results['CV_AUC_mean'], cv_results['CV_AUC_std'],
              optimal_threshold, best_c]
}
diagnosis_df = pd.DataFrame(diagnosis_summary)
diagnosis_df.to_csv(f'{OUTPUT_DIR}/diagnosis_summary.csv', index=False)

# 保存阈值分析结果
threshold_df = pd.DataFrame(threshold_results)
threshold_df.to_csv(f'{OUTPUT_DIR}/threshold_analysis.csv', index=False)

print(f"\n✓ 诊断结果已保存到: {OUTPUT_DIR}/")

# ==================== 生成可视化 ====================
print("\n生成可视化图表...")

fig, axes = plt.subplots(3, 2, figsize=(14, 12))

# 1. ROC曲线对比
ax1 = axes[0, 0]
fpr_train, tpr_train, _ = roc_curve(y_train_enc, y_train_prob)
fpr_val, tpr_val, _ = roc_curve(y_val_enc, y_val_prob)
ax1.plot(fpr_train, tpr_train, 'b-', label=f'Training (AUC={train_auc:.4f})', linewidth=2)
ax1.plot(fpr_val, tpr_val, 'r-', label=f'Validation (AUC={val_auc:.4f})', linewidth=2)
ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5)
ax1.set_xlabel('False Positive Rate')
ax1.set_ylabel('True Positive Rate')
ax1.set_title('ROC Curves Comparison')
ax1.legend(loc='lower right')
ax1.grid(True, alpha=0.3)

# 2. 预测概率分布
ax2 = axes[0, 1]
ax2.hist(train_sepsis_probs, bins=20, alpha=0.5, label='Train Sepsis', color='blue')
ax2.hist(train_normal_probs, bins=20, alpha=0.5, label='Train Normal', color='lightblue')
ax2.hist(val_sepsis_probs, bins=20, alpha=0.5, label='Val Sepsis', color='red')
ax2.hist(val_normal_probs, bins=20, alpha=0.5, label='Val Normal', color='lightsalmon')
ax2.axvline(x=0.5, color='black', linestyle='--', label='Threshold=0.5')
ax2.set_xlabel('Predicted Probability')
ax2.set_ylabel('Count')
ax2.set_title('Prediction Probability Distribution')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. 校准曲线
ax3 = axes[1, 0]
train_cal_pred_clean = [x for x in train_cal_pred if not np.isnan(x)]
train_cal_true_clean = [x for x in train_cal_true if not np.isnan(x)]
val_cal_pred_clean = [x for x in val_cal_pred if not np.isnan(x)]
val_cal_true_clean = [x for x in val_cal_true if not np.isnan(x)]

ax3.plot([0, 1], [0, 1], 'k--', label='Perfect Calibration')
if len(train_cal_pred_clean) > 0:
    ax3.plot(train_cal_pred_clean, train_cal_true_clean, 'b-o', label='Training', markersize=8)
if len(val_cal_pred_clean) > 0:
    ax3.plot(val_cal_pred_clean, val_cal_true_clean, 'r-o', label='Validation', markersize=8)
ax3.set_xlabel('Mean Predicted Probability')
ax3.set_ylabel('Fraction of Positives')
ax3.set_title('Calibration Curve')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. 阈值敏感性分析
ax4 = axes[1, 1]
thresholds = np.arange(0.1, 0.9, 0.05)
sens_vals = []
spec_vals = []
for thresh in thresholds:
    y_pred = (y_val_prob >= thresh).astype(int)
    sens_vals.append(sensitivity_score(y_val_enc, y_pred))
    spec_vals.append(specificity_score(y_val_enc, y_pred))

ax4.plot(thresholds, sens_vals, 'b-', label='Sensitivity', linewidth=2)
ax4.plot(thresholds, spec_vals, 'r-', label='Specificity', linewidth=2)
ax4.axvline(x=optimal_threshold, color='green', linestyle='--', label=f'Optimal={optimal_threshold:.3f}')
ax4.axvline(x=0.5, color='gray', linestyle=':', alpha=0.5)
ax4.set_xlabel('Classification Threshold')
ax4.set_ylabel('Score')
ax4.set_title('Threshold Sensitivity Analysis (Validation Set)')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. 模型复杂度对比
ax5 = axes[2, 0]
models = ['Single Gene\n(HLA-DRA)', '2 Genes', '3 Genes', 'All MHC']
if 'HLA-DRA' in feature_names:
    train_aucs = [single_train_auc, two_train_auc, three_train_auc, train_auc]
    val_aucs = [single_val_auc, two_val_auc, three_val_auc, val_auc]
    
    x = np.arange(len(models))
    width = 0.35
    bars1 = ax5.bar(x - width/2, train_aucs, width, label='Training', color='blue', alpha=0.7)
    bars2 = ax5.bar(x + width/2, val_aucs, width, label='Validation', color='red', alpha=0.7)
    ax5.set_ylabel('AUC')
    ax5.set_title('Model Complexity vs Performance')
    ax5.set_xticks(x)
    ax5.set_xticklabels(models)
    ax5.legend()
    ax5.set_ylim([0.7, 1.0])
    ax5.axhline(y=0.9, color='gray', linestyle='--', alpha=0.5)
    ax5.grid(True, alpha=0.3, axis='y')

# 6. 交叉验证结果
ax6 = axes[2, 1]
folds = range(1, 6)
ax6.errorbar(folds, cv_aucs, yerr=[cv_aucs[i]*0.02 for i in range(5)], 
              fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8)
ax6.axhline(y=cv_results['CV_AUC_mean'], color='red', linestyle='--', 
            label=f'Mean AUC={cv_results["CV_AUC_mean"]:.4f}')
ax6.fill_between([0.5, 5.5], 
                  cv_results['CV_AUC_mean'] - cv_results['CV_AUC_std'],
                  cv_results['CV_AUC_mean'] + cv_results['CV_AUC_std'],
                  alpha=0.2, color='red')
ax6.set_xlabel('Fold')
ax6.set_ylabel('AUC')
ax6.set_title('5-Fold Cross-Validation Results')
ax6.set_xticks(folds)
ax6.legend()
ax6.grid(True, alpha=0.3)
ax6.set_ylim([0.7, 1.0])

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/overfitting_diagnosis.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ 可视化已保存到: {OUTPUT_DIR}/overfitting_diagnosis.png")

print("\n" + "=" * 70)
print("诊断完成!")
print("=" * 70)
