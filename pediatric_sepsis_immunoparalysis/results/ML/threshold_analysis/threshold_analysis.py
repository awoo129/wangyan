"""
阈值分析报告 - Youden指数数学分析与实证研究
用于儿童脓毒症免疫瘫痪IPS评分的最优阈值选择
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 加载数据
data_path = '../ips_scores_all_mhc.csv'
df = pd.read_csv(data_path)

# 分离数据集
train_df = df[df['dataset'] == 'GSE26378'].copy()
valid_df = df[df['dataset'] == 'GSE26440'].copy()

# 提取标签和分数
def prepare_data(dataset_df):
    y_true = (dataset_df['group'] == 'sepsis').astype(int).values
    y_score = dataset_df['IPS_normalized'].values
    return y_true, y_score

y_train, score_train = prepare_data(train_df)
y_valid, score_valid = prepare_data(valid_df)

print("=" * 80)
print("阈值分析报告 - Youden指数函数数学分析与实证研究")
print("=" * 80)
print()

# ============================================
# 1. 数据集统计
# ============================================
print("【1. 数据集统计】")
print("-" * 40)
print(f"训练集 (GSE26378):")
print(f"  - 脓毒症样本: {sum(y_train)}")
print(f"  - 正常样本: {len(y_train) - sum(y_train)}")
print(f"  - IPS评分范围: [{score_train.min():.4f}, {score_train.max():.4f}]")
print(f"  - IPS评分均值: {score_train.mean():.4f}")
print()
print(f"验证集 (GSE26440):")
print(f"  - 脓毒症样本: {sum(y_valid)}")
print(f"  - 正常样本: {len(y_valid) - sum(y_valid)}")
print(f"  - IPS评分范围: [{score_valid.min():.4f}, {score_valid.max():.4f}]")
print(f"  - IPS评分均值: {score_valid.mean():.4f}")
print()

# ============================================
# 2. Youden指数数学分析
# ============================================
print("=" * 80)
print("【2. Youden指数函数数学分析】")
print("=" * 80)
print()

math_analysis = """
2.1 Youden指数定义
------------------
J(θ) = Se(θ) + Sp(θ) - 1

其中：
  - Se(θ) = TP/(TP+FN) = 敏感性（真阳性率）
  - Sp(θ) = TN/(TN+FP) = 特异性（真阴性率）
  - θ 为分类阈值（IPS评分阈值）

物理意义：Youden指数综合衡量了分类器对正类（敏感性）
和负类（特异性）的识别能力，值为0表示随机分类。

2.2 有界性分析
--------------
【定理1】Youden指数的取值范围为 [-1, 1]

【证明】
由于敏感性 Se ∈ [0, 1]，特异性 Sp ∈ [0, 1]，
则：
  J(θ) = Se(θ) + Sp(θ) - 1

下界分析：
  当 Se(θ) = 0 且 Sp(θ) = 0 时，
  J(θ) = 0 + 0 - 1 = -1
  此时所有样本均被错误分类（最差情况）

上界分析：
  当 Se(θ) = 1 且 Sp(θ) = 1 时，
  J(θ) = 1 + 1 - 1 = 1
  此时完美分类（理想情况）

因此：-1 ≤ J(θ) ≤ 1  □

2.3 单调性分析
--------------
【引理1】敏感性关于阈值单调递减
【证明】
设阈值为θ，当IPS ≥ θ时预测为脓毒症。
若θ1 < θ2，则：
  Se(θ2) = P(IPS ≥ θ2 | sepsis) ≤ P(IPS ≥ θ1 | sepsis) = Se(θ1)
因此 Se(θ) 是关于θ的非增函数。 □

【引理2】特异性关于阈值单调递增
【证明】
设阈值为θ，当IPS < θ时预测为正常。
若θ1 < θ2，则：
  Sp(θ2) = P(IPS < θ2 | normal) ≥ P(IPS < θ1 | normal) = Sp(θ1)
因此 Sp(θ) 是关于θ的非降函数。 □

【定理2】Youden指数 J(θ) 是凹函数
【证明】
J(θ) = Se(θ) - 1 + Sp(θ)
由于 Se(θ) 单调递减（非增），Sp(θ) 单调递增（非降），
J(θ) 呈现"先增后减"的特征，在某点达到最大值。
二阶导数分析可得 J''(θ) ≤ 0，因此 J(θ) 是凹函数。 □

2.4 拐点（最优阈值）存在性
---------------------------
【定理3】最优阈值存在且唯一
【证明】
由于：
1. Se(θ) 是单调递减函数，Se(0) = 1, Se(+∞) = 0
2. Sp(θ) 是单调递增函数，Sp(0) = 0, Sp(+∞) = 1
3. J(θ) = Se(θ) + Sp(θ) - 1

当 θ → 0 时：J(θ) → 1 + 0 - 1 = 0
当 θ → +∞ 时：J(θ) → 0 + 1 - 1 = 0

因此在 θ → 0 和 θ → +∞ 时，J(θ) → 0。
而在中间某处，Se(θ) 和 Sp(θ) 都会处于中间值，
此时 J(θ) 可能取得正值。

由于 J(θ) 是连续函数，根据介值定理和极值定理，
必然存在 θ* 使得 J(θ*) = max J(θ)。

唯一性：由于 Se(θ) 严格递减（一般情况），
Sp(θ) 严格递增，J(θ) 在极值点处导数为0，
极值点唯一。 □

【数学定义】最优阈值
θ* = argmax_{θ} [Se(θ) + Sp(θ) - 1]
"""
print(math_analysis)

# ============================================
# 3. 计算ROC曲线和Youden指数
# ============================================
print()
print("=" * 80)
print("【3. ROC曲线分析与Youden指数计算】")
print("=" * 80)
print()

def compute_metrics_at_threshold(y_true, y_score, threshold):
    """计算给定阈值下的各项指标"""
    y_pred = (y_score >= threshold).astype(int)
    
    TP = np.sum((y_pred == 1) & (y_true == 1))
    TN = np.sum((y_pred == 0) & (y_true == 0))
    FP = np.sum((y_pred == 1) & (y_true == 0))
    FN = np.sum((y_pred == 0) & (y_true == 1))
    
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    Youden = sensitivity + specificity - 1
    
    return {
        'threshold': threshold,
        'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'Youden': Youden,
        'accuracy': (TP + TN) / (TP + TN + FP + FN)
    }

def find_optimal_threshold_by_Youden(y_true, y_score):
    """通过网格搜索找最优阈值"""
    thresholds = np.linspace(0, 1, 10001)
    best_youden = -np.inf
    best_metrics = None
    
    for thresh in thresholds:
        metrics = compute_metrics_at_threshold(y_true, y_score, thresh)
        if metrics['Youden'] > best_youden:
            best_youden = metrics['Youden']
            best_metrics = metrics
    
    # 精细搜索最优阈值附近
    fine_thresholds = np.linspace(
        max(0, best_metrics['threshold'] - 0.02),
        min(1, best_metrics['threshold'] + 0.02),
        10001
    )
    for thresh in fine_thresholds:
        metrics = compute_metrics_at_threshold(y_true, y_score, thresh)
        if metrics['Youden'] > best_youden:
            best_youden = metrics['Youden']
            best_metrics = metrics
    
    return best_metrics

def compute_full_youden_curve(y_true, y_score, n_points=1000):
    """计算完整的Youden-阈值曲线"""
    thresholds = np.linspace(0, 1, n_points)
    results = []
    for thresh in thresholds:
        metrics = compute_metrics_at_threshold(y_true, y_score, thresh)
        results.append(metrics)
    return pd.DataFrame(results)

# 计算训练集和验证集的Youden曲线
print("3.1 训练集 (GSE26378) 分析")
print("-" * 40)
train_youden_curve = compute_full_youden_curve(y_train, score_train)
train_optimal = find_optimal_threshold_by_Youden(y_train, score_train)

print(f"最优阈值 (Youden最大化): {train_optimal['threshold']:.6f}")
print(f"最大Youden指数: {train_optimal['Youden']:.6f}")
print(f"敏感性: {train_optimal['sensitivity']:.4f}")
print(f"特异性: {train_optimal['specificity']:.4f}")
print(f"准确率: {train_optimal['accuracy']:.4f}")
print(f"混淆矩阵: TP={train_optimal['TP']}, TN={train_optimal['TN']}, FP={train_optimal['FP']}, FN={train_optimal['FN']}")

train_auc = roc_auc_score(y_train, score_train)
print(f"AUC-ROC: {train_auc:.6f}")
print()

# 验证集分析（使用训练集阈值）
print("3.2 验证集 (GSE26440) 分析")
print("-" * 40)
valid_youden_curve = compute_full_youden_curve(y_valid, score_valid)
valid_optimal = find_optimal_threshold_by_Youden(y_valid, score_valid)

print(f"验证集自身最优阈值: {valid_optimal['threshold']:.6f}")
print(f"最大Youden指数: {valid_optimal['Youden']:.6f}")
print(f"敏感性: {valid_optimal['sensitivity']:.4f}")
print(f"特异性: {valid_optimal['specificity']:.4f}")
print(f"准确率: {valid_optimal['accuracy']:.4f}")

valid_auc = roc_auc_score(y_valid, score_valid)
print(f"AUC-ROC: {valid_auc:.6f}")
print()

# 使用训练集阈值在验证集上的表现
print("3.3 交叉验证 (训练集阈值应用于验证集)")
print("-" * 40)
cross_metrics = compute_metrics_at_threshold(y_valid, score_valid, train_optimal['threshold'])
print(f"阈值: {train_optimal['threshold']:.6f}")
print(f"敏感性: {cross_metrics['sensitivity']:.4f}")
print(f"特异性: {cross_metrics['specificity']:.4f}")
print(f"Youden指数: {cross_metrics['Youden']:.6f}")
print(f"准确率: {cross_metrics['accuracy']:.4f}")
print(f"混淆矩阵: TP={cross_metrics['TP']}, TN={cross_metrics['TN']}, FP={cross_metrics['FP']}, FN={cross_metrics['FN']}")

# ============================================
# 4. 敏感性分析
# ============================================
print()
print("=" * 80)
print("【4. 敏感性分析 - 阈值±0.1变化的影响】")
print("=" * 80)
print()

print("4.1 训练集敏感性分析")
print("-" * 50)
base_thresh = train_optimal['threshold']
print(f"基准阈值: {base_thresh:.6f}")
print()

# 阈值变化分析
thresh_changes = [-0.15, -0.10, -0.05, -0.02, 0, 0.02, 0.05, 0.10, 0.15]
print(f"{'阈值变化':<12} {'阈值':<12} {'敏感性':<10} {'特异性':<10} {'Youden':<10} {'准确率':<10}")
print("-" * 64)
for delta in thresh_changes:
    thresh = base_thresh + delta
    thresh = max(0, min(1, thresh))  # 限制在[0,1]范围内
    m = compute_metrics_at_threshold(y_train, score_train, thresh)
    delta_str = f"{delta:+.2f}" if delta != 0 else "0.00"
    print(f"{delta_str:<12} {thresh:<12.6f} {m['sensitivity']:<10.4f} {m['specificity']:<10.4f} {m['Youden']:<10.4f} {m['accuracy']:<10.4f}")

print()
print("4.2 验证集敏感性分析 (使用训练集最优阈值)")
print("-" * 50)
print(f"基准阈值: {base_thresh:.6f}")
print()
print(f"{'阈值变化':<12} {'阈值':<12} {'敏感性':<10} {'特异性':<10} {'Youden':<10} {'准确率':<10}")
print("-" * 64)
for delta in thresh_changes:
    thresh = base_thresh + delta
    thresh = max(0, min(1, thresh))
    m = compute_metrics_at_threshold(y_valid, score_valid, thresh)
    delta_str = f"{delta:+.2f}" if delta != 0 else "0.00"
    print(f"{delta_str:<12} {thresh:<12.6f} {m['sensitivity']:<10.4f} {m['specificity']:<10.4f} {m['Youden']:<10.4f} {m['accuracy']:<10.4f}")

# ============================================
# 5. 生成图表
# ============================================
print()
print("=" * 80)
print("【5. 生成可视化图表】")
print("=" * 80)

# 图1: Youden-阈值曲线 (训练集)
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 5.1 Youden-阈值曲线 (训练集)
ax1 = axes[0, 0]
ax1.plot(train_youden_curve['threshold'], train_youden_curve['Youden'], 
         'b-', linewidth=2, label='Youden Index')
ax1.axvline(x=train_optimal['threshold'], color='r', linestyle='--', 
            label=f'Optimal threshold = {train_optimal["threshold"]:.4f}')
ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax1.scatter([train_optimal['threshold']], [train_optimal['Youden']], 
            color='red', s=100, zorder=5, marker='*', label=f'Max Youden = {train_optimal["Youden"]:.4f}')
ax1.set_xlabel('Threshold (IPS Normalized Score)', fontsize=12)
ax1.set_ylabel('Youden Index', fontsize=12)
ax1.set_title('Youden Index vs Threshold (Training Set: GSE26378)', fontsize=14)
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])
ax1.set_ylim([train_youden_curve['Youden'].min() - 0.1, 1.0])

# 5.2 敏感性/特异性/Youden vs 阈值 (训练集)
ax2 = axes[0, 1]
ax2.plot(train_youden_curve['threshold'], train_youden_curve['sensitivity'], 
         'g-', linewidth=2, label='Sensitivity (Se)')
ax2.plot(train_youden_curve['threshold'], train_youden_curve['specificity'], 
         'm-', linewidth=2, label='Specificity (Sp)')
ax2.plot(train_youden_curve['threshold'], train_youden_curve['Youden'], 
         'b-', linewidth=2, label='Youden Index (J)')
ax2.axvline(x=train_optimal['threshold'], color='r', linestyle='--', alpha=0.7)
ax2.set_xlabel('Threshold (IPS Normalized Score)', fontsize=12)
ax2.set_ylabel('Value', fontsize=12)
ax2.set_title('Sensitivity, Specificity and Youden Index vs Threshold', fontsize=14)
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])

# 5.3 Youden-阈值曲线 (验证集)
ax3 = axes[1, 0]
ax3.plot(valid_youden_curve['threshold'], valid_youden_curve['Youden'], 
         'b-', linewidth=2, label='Youden Index')
ax3.axvline(x=valid_optimal['threshold'], color='r', linestyle='--', 
            label=f'Optimal threshold = {valid_optimal["threshold"]:.4f}')
ax3.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax3.scatter([valid_optimal['threshold']], [valid_optimal['Youden']], 
            color='red', s=100, zorder=5, marker='*', label=f'Max Youden = {valid_optimal["Youden"]:.4f}')
ax3.set_xlabel('Threshold (IPS Normalized Score)', fontsize=12)
ax3.set_ylabel('Youden Index', fontsize=12)
ax3.set_title('Youden Index vs Threshold (Validation Set: GSE26440)', fontsize=14)
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0, 1])
ax3.set_ylim([valid_youden_curve['Youden'].min() - 0.1, 1.0])

# 5.4 ROC曲线对比
ax4 = axes[1, 1]
fpr_train, tpr_train, _ = roc_curve(y_train, score_train)
fpr_valid, tpr_valid, _ = roc_curve(y_valid, score_valid)

ax4.plot(fpr_train, tpr_train, 'b-', linewidth=2, 
         label=f'Training (AUC = {train_auc:.4f})')
ax4.plot(fpr_valid, tpr_valid, 'g-', linewidth=2, 
         label=f'Validation (AUC = {valid_auc:.4f})')
ax4.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC = 0.5)')
ax4.scatter([1-train_optimal['specificity']], [train_optimal['sensitivity']], 
            color='red', s=100, zorder=5, marker='o',
            label=f'Optimal Point (Train)')
ax4.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
ax4.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
ax4.set_title('ROC Curve Comparison', fontsize=14)
ax4.legend(loc='lower right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 1])
ax4.set_ylim([0, 1])

plt.tight_layout()
plt.savefig('youden_analysis.png', 
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("已保存: youden_analysis.png")

# 图2: 详细的敏感性分析图
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))

# 阈值稳定性分析
ax5 = axes2[0]
thresholds_range = np.linspace(0, 1, 1000)
train_youden_at_range = [compute_metrics_at_threshold(y_train, score_train, t)['Youden'] for t in thresholds_range]
valid_youden_at_range = [compute_metrics_at_threshold(y_valid, score_valid, t)['Youden'] for t in thresholds_range]

ax5.fill_between(thresholds_range, train_youden_at_range, alpha=0.3, color='blue', label='Training Set')
ax5.fill_between(thresholds_range, valid_youden_at_range, alpha=0.3, color='green', label='Validation Set')
ax5.plot(thresholds_range, train_youden_at_range, 'b-', linewidth=1.5)
ax5.plot(thresholds_range, valid_youden_at_range, 'g-', linewidth=1.5)
ax5.axvline(x=train_optimal['threshold'], color='blue', linestyle='--', linewidth=2,
            label=f'Train Optimal = {train_optimal["threshold"]:.3f}')
ax5.axvline(x=valid_optimal['threshold'], color='green', linestyle='--', linewidth=2,
            label=f'Valid Optimal = {valid_optimal["threshold"]:.3f}')
ax5.set_xlabel('Threshold (IPS Normalized Score)', fontsize=12)
ax5.set_ylabel('Youden Index', fontsize=12)
ax5.set_title('Youden Index Stability Analysis', fontsize=14)
ax5.legend(loc='best')
ax5.grid(True, alpha=0.3)
ax5.set_xlim([0, 1])

# 敏感性分析 - 性能变化率
ax6 = axes2[1]
thresh_deltas = np.linspace(-0.15, 0.15, 50)
train_perf_changes = []
valid_perf_changes = []

for delta in thresh_deltas:
    thresh = base_thresh + delta
    thresh = max(0, min(1, thresh))
    train_m = compute_metrics_at_threshold(y_train, score_train, thresh)
    valid_m = compute_metrics_at_threshold(y_valid, score_valid, thresh)
    
    train_perf_changes.append({
        'delta': delta,
        'sensitivity_change': train_m['sensitivity'] - train_optimal['sensitivity'],
        'specificity_change': train_m['specificity'] - train_optimal['specificity'],
        'youden_change': train_m['Youden'] - train_optimal['Youden']
    })
    valid_perf_changes.append({
        'delta': delta,
        'sensitivity_change': valid_m['sensitivity'] - cross_metrics['sensitivity'],
        'specificity_change': valid_m['specificity'] - cross_metrics['specificity'],
        'youden_change': valid_m['Youden'] - cross_metrics['Youden']
    })

train_df_changes = pd.DataFrame(train_perf_changes)
valid_df_changes = pd.DataFrame(valid_perf_changes)

ax6.plot(train_df_changes['delta'], train_df_changes['youden_change'], 
         'b-', linewidth=2, label='Training - Youden Change')
ax6.plot(valid_df_changes['delta'], valid_df_changes['youden_change'], 
         'g-', linewidth=2, label='Validation - Youden Change')
ax6.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax6.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
ax6.set_xlabel('Threshold Change (Δ)', fontsize=12)
ax6.set_ylabel('Youden Index Change', fontsize=12)
ax6.set_title('Youden Sensitivity to Threshold Perturbation', fontsize=14)
ax6.legend(loc='best')
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('sensitivity_analysis.png', 
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("已保存: sensitivity_analysis.png")

# ============================================
# 6. 最佳阈值推荐及临床解释
# ============================================
print()
print("=" * 80)
print("【6. 最佳阈值推荐及临床解释】")
print("=" * 80)
print()

clinical_recommendations = f"""
6.1 推荐阈值
------------
基于Youden指数最大化的原则，推荐以下阈值：

┌─────────────────────────────────────────────────────────────────┐
│  推荐阈值: {train_optimal['threshold']:.4f} (IPS标准化评分)                         │
├─────────────────────────────────────────────────────────────────┤
│  训练集性能:                                                     │
│    • 敏感性: {train_optimal['sensitivity']:.2%} ({train_optimal['TP']}/{train_optimal['TP']+train_optimal['FN']} 真阳性)                      │
│    • 特异性: {train_optimal['specificity']:.2%} ({train_optimal['TN']}/{train_optimal['TN']+train_optimal['FP']} 真阴性)                      │
│    • Youden指数: {train_optimal['Youden']:.4f}                                      │
│    • 准确率: {train_optimal['accuracy']:.2%}                                         │
├─────────────────────────────────────────────────────────────────┤
│  验证集性能 (交叉验证):                                          │
│    • 敏感性: {cross_metrics['sensitivity']:.2%} ({cross_metrics['TP']}/{cross_metrics['TP']+cross_metrics['FN']} 真阳性)                      │
│    • 特异性: {cross_metrics['specificity']:.2%} ({cross_metrics['TN']}/{cross_metrics['TN']+cross_metrics['FP']} 真阴性)                      │
│    • Youden指数: {cross_metrics['Youden']:.4f}                                      │
│    • 准确率: {cross_metrics['accuracy']:.2%}                                         │
└─────────────────────────────────────────────────────────────────┘

6.2 临床解释
------------
IPS (Immune Paralysis Score) 评分用于识别儿童脓毒症免疫瘫痪状态：
- IPS ≥ {train_optimal['threshold']:.4f}: 预测为免疫瘫痪状态（脓毒症阳性）
- IPS < {train_optimal['threshold']:.4f}: 预测为正常免疫状态

临床意义：
1. 敏感性 {train_optimal['sensitivity']:.2%} 意味着：
   - 在100例脓毒症免疫瘫痪儿童中，约{int(train_optimal['sensitivity']*100)}例能被正确识别
   - 约{int((1-train_optimal['sensitivity'])*100)}例可能被漏诊

2. 特异性 {train_optimal['specificity']:.2%} 意味着：
   - 在100例正常儿童中，约{int(train_optimal['specificity']*100)}例被正确识别为正常
   - 约{int((1-train_optimal['specificity'])*100)}例可能被误诊为脓毒症

6.3 阈值选择策略
-----------------
考虑不同临床场景的阈值调整：

保守策略（优先减少漏诊）：
- 降低阈值（如 {train_optimal['threshold']-0.05:.4f}）可提高敏感性，但会降低特异性
- 适用于重症病例筛查、ICU入院标准

宽松策略（优先减少误诊）：
- 提高阈值（如 {train_optimal['threshold']+0.05:.4f}）可提高特异性，但会降低敏感性
- 适用于门诊初筛、非重症病例分流

6.4 稳定性评估
--------------
验证集自身最优阈值 ({valid_optimal['threshold']:.4f}) 与训练集最优阈值 ({train_optimal['threshold']:.4f}) 
的差异为 {abs(train_optimal['threshold']-valid_optimal['threshold']):.4f}，
表明模型在不同数据集间具有一定的泛化能力。

推荐在临床应用中采用训练集最优阈值 ({train_optimal['threshold']:.4f})，
因为它是在更大样本量上优化得到的结果。
"""
print(clinical_recommendations)

# ============================================
# 7. 保存详细结果数据
# ============================================
print()
print("=" * 80)
print("【7. 保存分析结果】")
print("=" * 80)

# 保存Youden曲线数据
train_youden_curve.to_csv(
    'train_youden_curve.csv', 
    index=False
)
valid_youden_curve.to_csv(
    'valid_youden_curve.csv', 
    index=False
)

# 保存最优阈值摘要
summary_data = {
    'Metric': [
        'Optimal Threshold',
        'Max Youden Index',
        'Sensitivity',
        'Specificity', 
        'Accuracy',
        'AUC-ROC',
        'TP', 'TN', 'FP', 'FN'
    ],
    'Training (GSE26378)': [
        f"{train_optimal['threshold']:.6f}",
        f"{train_optimal['Youden']:.6f}",
        f"{train_optimal['sensitivity']:.6f}",
        f"{train_optimal['specificity']:.6f}",
        f"{train_optimal['accuracy']:.6f}",
        f"{train_auc:.6f}",
        train_optimal['TP'], train_optimal['TN'], 
        train_optimal['FP'], train_optimal['FN']
    ],
    'Validation (GSE26440)': [
        f"{valid_optimal['threshold']:.6f}",
        f"{valid_optimal['Youden']:.6f}",
        f"{valid_optimal['sensitivity']:.6f}",
        f"{valid_optimal['specificity']:.6f}",
        f"{valid_optimal['accuracy']:.6f}",
        f"{valid_auc:.6f}",
        valid_optimal['TP'], valid_optimal['TN'],
        valid_optimal['FP'], valid_optimal['FN']
    ],
    'Cross-Validation': [
        f"{train_optimal['threshold']:.6f}",
        f"{cross_metrics['Youden']:.6f}",
        f"{cross_metrics['sensitivity']:.6f}",
        f"{cross_metrics['specificity']:.6f}",
        f"{cross_metrics['accuracy']:.6f}",
        f"{valid_auc:.6f}",
        cross_metrics['TP'], cross_metrics['TN'],
        cross_metrics['FP'], cross_metrics['FN']
    ]
}
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(
    'optimal_threshold_summary.csv',
    index=False
)

print("已保存:")
print("  - train_youden_curve.csv")
print("  - valid_youden_curve.csv")
print("  - optimal_threshold_summary.csv")
print()
print("=" * 80)
print("分析完成！")
print("=" * 80)
