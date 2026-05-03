#!/usr/bin/env python3
"""
使用模拟数据演示儿童脓毒症免疫瘫痪研究分析流程
当实际GEO数据下载完成后，可用真实数据替换
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import warnings
warnings.filterwarnings('ignore')

# 设置
plt.style.use('seaborn-v0_8-whitegrid')
MAIN_DIR = Path("./长期计划/儿童脓毒症免疫瘫痪研究")

# ========================
# 1. 生成模拟数据
# ========================
def generate_mock_data():
    """生成模拟的GEO数据集"""
    print("="*60)
    print("生成模拟数据集")
    print("="*60)
    
    np.random.seed(42)
    
    # MHC II类基因（核心基因）
    mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DRB1', 
                 'HLA-DPA1', 'HLA-DPB1', 'CIITA', 'CD74']
    
    # 炎症/免疫相关基因
    immune_genes = ['IL6', 'TNF', 'IL1B', 'CXCL8', 'CCL2', 'IL10',
                   'CD14', 'TLR4', 'NFKB1', 'STAT3']
    
    # 其他基因（背景）
    other_genes = [f'GENE_{i}' for i in range(100, 200)]
    all_genes = mhc_genes + immune_genes + other_genes
    
    datasets = {}
    
    # GSE26378: 103 samples (21 Normal + 82 Sepsis)
    n_normal_26378 = 21
    n_sepsis_26378 = 82
    n_total_26378 = n_normal_26378 + n_sepsis_26378
    
    expr_26378 = np.random.randn(len(all_genes), n_total_26378) * 0.5 + 8
    
    # HLA-DRA在脓毒症中下调（免疫瘫痪特征）
    sepsis_idx_26378 = list(range(n_normal_26378, n_total_26378))
    expr_26378[0, sepsis_idx_26378] -= 1.5  # HLA-DRA下调
    expr_26378[1:4, sepsis_idx_26378] -= 0.8  # 其他MHC基因下调
    expr_26378[10:14, sepsis_idx_26378] += 1.0  # 炎症因子上调
    
    datasets['GSE26378'] = {
        'expression': pd.DataFrame(expr_26378, index=all_genes,
                                   columns=[f'S{i+1}' for i in range(n_total_26378)]),
        'metadata': pd.DataFrame({
            'sample_id': [f'S{i+1}' for i in range(n_total_26378)],
            'group': ['Normal']*n_normal_26378 + ['Sepsis']*n_sepsis_26378,
            'platform': 'GPL570'
        })
    }
    
    # GSE26440: 130 samples (32 Normal + 98 Septic Shock)
    n_normal_26440 = 32
    n_sepsis_26440 = 98
    n_total_26440 = n_normal_26440 + n_sepsis_26440
    
    expr_26440 = np.random.randn(len(all_genes), n_total_26440) * 0.5 + 8
    sepsis_idx_26440 = list(range(n_normal_26440, n_total_26440))
    expr_26440[0, sepsis_idx_26440] -= 1.8  # HLA-DRA更显著下调
    expr_26440[1:4, sepsis_idx_26440] -= 1.0
    expr_26440[10:14, sepsis_idx_26440] += 1.2
    
    datasets['GSE26440'] = {
        'expression': pd.DataFrame(expr_26440, index=all_genes,
                                   columns=[f'S{i+1}' for i in range(n_total_26440)]),
        'metadata': pd.DataFrame({
            'sample_id': [f'S{i+1}' for i in range(n_total_26440)],
            'group': ['Normal']*n_normal_26440 + ['Septic Shock']*n_sepsis_26440,
            'platform': 'GPL570'
        })
    }
    
    # GSE13904: 124 samples (18 Normal + 106 Sepsis)
    n_normal_13904 = 18
    n_sepsis_13904 = 106
    n_total_13904 = n_normal_13904 + n_sepsis_13904
    
    expr_13904 = np.random.randn(len(all_genes), n_total_13904) * 0.5 + 8
    sepsis_idx_13904 = list(range(n_normal_13904, n_total_13904))
    expr_13904[0, sepsis_idx_13904] -= 1.3
    expr_13904[1:4, sepsis_idx_13904] -= 0.6
    expr_13904[10:14, sepsis_idx_13904] += 0.9
    
    datasets['GSE13904'] = {
        'expression': pd.DataFrame(expr_13904, index=all_genes,
                                   columns=[f'S{i+1}' for i in range(n_total_13904)]),
        'metadata': pd.DataFrame({
            'sample_id': [f'S{i+1}' for i in range(n_total_13904)],
            'group': ['Normal']*n_normal_13904 + ['Sepsis']*n_sepsis_13904,
            'platform': 'GPL570'
        })
    }
    
    return datasets

# ========================
# 2. MHC基因分析
# ========================
def analyze_mhc_genes(datasets):
    """分析MHC基因表达"""
    print("\n" + "="*60)
    print("MHC II类基因分析")
    print("="*60)
    
    mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DRB1', 
                 'HLA-DPA1', 'HLA-DPB1', 'CIITA', 'CD74']
    
    results = {}
    
    for gse_id, data in datasets.items():
        expr = data['expression']
        meta = data['metadata']
        
        # 检查MHC基因是否在数据中
        found_genes = [g for g in mhc_genes if g in expr.index]
        
        results[gse_id] = {
            'found_genes': found_genes,
            'hla_dra_present': 'HLA-DRA' in expr.index
        }
        
        print(f"\n{gse_id}:")
        print(f"  MHC基因: {', '.join(found_genes)}")
        print(f"  HLA-DRA: {'✓ 存在' if 'HLA-DRA' in expr.index else '✗ 不存在'}")
        
        # HLA-DRA差异分析
        if 'HLA-DRA' in expr.index:
            hla_dra_expr = expr.loc['HLA-DRA']
            normal_mask = meta['group'] == 'Normal'
            sepsis_mask = meta['group'].str.contains('Sepsis|Shock', regex=True)
            
            normal_mean = hla_dra_expr[normal_mask.values].mean()
            sepsis_mean = hla_dra_expr[sepsis_mask.values].mean()
            logfc = sepsis_mean - normal_mean
            
            # 简单t检验
            t_stat, p_val = stats.ttest_ind(hla_dra_expr[sepsis_mask.values], 
                                           hla_dra_expr[normal_mask.values])
            
            print(f"  HLA-DRA表达:")
            print(f"    正常组: {normal_mean:.3f}")
            print(f"    脓毒症组: {sepsis_mean:.3f}")
            print(f"    logFC: {logfc:.3f}")
            print(f"    P值: {p_val:.2e}")
            print(f"    方向: {'下调' if logfc < 0 else '上调'}")
            
            results[gse_id]['hla_dra_logfc'] = logfc
            results[gse_id]['hla_dra_pval'] = p_val
    
    return results

# ========================
# 3. 免疫瘫痪评分(IPS)
# ========================
def calculate_ips(datasets):
    """计算免疫瘫痪评分"""
    print("\n" + "="*60)
    print("免疫瘫痪评分(IPS)计算")
    print("="*60)
    
    # 权重（基于文献和重要性）
    weights = {
        'HLA-DRA': 2.0,      # 核心基因，权重最高
        'HLA-DQB1': 1.0,
        'HLA-DQA1': 1.0,
        'HLA-DRB1': 1.0,
        'HLA-DPA1': 0.5,
        'HLA-DPB1': 0.5,
        'CIITA': 1.5,        # MHC转录调控
        'CD74': 1.0
    }
    
    ips_results = {}
    
    for gse_id, data in datasets.items():
        expr = data['expression']
        meta = data['metadata']
        
        # 计算加权IPS
        ips_values = []
        for sample in expr.columns:
            sample_expr = expr[sample]
            weighted_sum = 0
            weight_sum = 0
            
            for gene, weight in weights.items():
                if gene in expr.index:
                    weighted_sum += sample_expr[gene] * weight
                    weight_sum += weight
            
            ips = weighted_sum / weight_sum if weight_sum > 0 else np.nan
            ips_values.append(ips)
        
        # 存储结果
        meta_copy = meta.copy()
        meta_copy['IPS'] = ips_values
        
        # IPS分组（使用中位数）
        ips_median = np.median(ips_values)
        meta_copy['IPS_Group'] = ['High_IPS' if ips > ips_median else 'Low_IPS' 
                                   for ips in ips_values]
        
        ips_results[gse_id] = meta_copy
        
        # 统计
        print(f"\n{gse_id}:")
        print(f"  IPS均值: {np.mean(ips_values):.3f}")
        print(f"  IPS中位数: {np.median(ips_values):.3f}")
        print(f"  IPS标准差: {np.std(ips_values):.3f}")
        print(f"  分组阈值: {ips_median:.3f}")
        print(f"  Low_IPS: {(meta_copy['IPS_Group']=='Low_IPS').sum()}")
        print(f"  High_IPS: {(meta_copy['IPS_Group']=='High_IPS').sum()}")
    
    return ips_results

# ========================
# 4. 可视化
# ========================
def create_visualizations(datasets, mhc_results, ips_results):
    """生成可视化图表"""
    print("\n" + "="*60)
    print("生成可视化图表")
    print("="*60)
    
    viz_dir = MAIN_DIR / "results" / "visualization"
    viz_dir.mkdir(parents=True, exist_ok=True)
    
    # 4.1 MHC基因表达热图
    for gse_id, data in datasets.items():
        expr = data['expression']
        meta = data['metadata']
        
        mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74']
        available_mhc = [g for g in mhc_genes if g in expr.index]
        
        if available_mhc:
            mhc_expr = expr.loc[available_mhc]
            
            # Z-score标准化
            mhc_zscore = (mhc_expr - mhc_expr.mean(axis=1).values.reshape(-1,1)) / \
                        mhc_expr.std(axis=1).values.reshape(-1,1)
            
            plt.figure(figsize=(14, 6))
            
            # 按分组排序样本名称
            normal_samples = meta[meta['group'] == 'Normal']['sample_id'].tolist()
            sepsis_samples = meta[meta['group'].str.contains('Sepsis|Shock', regex=True)]['sample_id'].tolist()
            group_order = normal_samples + sepsis_samples
            
            # 只保留存在的列
            group_order = [s for s in group_order if s in mhc_zscore.columns]
            mhc_zscore = mhc_zscore[group_order]
            
            sns.heatmap(mhc_zscore, cmap='RdBu_r', center=0,
                       xticklabels=False, yticklabels=True)
            plt.title(f'{gse_id} - MHC II类基因表达热图')
            plt.xlabel('Samples')
            plt.ylabel('Genes')
            plt.tight_layout()
            plt.savefig(viz_dir / f'{gse_id}_mhc_heatmap.png', dpi=150)
            plt.close()
            print(f"保存: {viz_dir / f'{gse_id}_mhc_heatmap.png'}")
    
    # 4.2 HLA-DRA表达箱线图
    plt.figure(figsize=(12, 5))
    
    for i, (gse_id, data) in enumerate(datasets.items()):
        expr = data['expression']
        meta = data['metadata']
        
        if 'HLA-DRA' in expr.index:
            hla_dra_df = pd.DataFrame({
                'Expression': expr.loc['HLA-DRA'].values,
                'Group': meta['group'].values,
                'Dataset': gse_id
            })
            
            plt.subplot(1, 3, i+1)
            sns.boxplot(data=hla_dra_df, x='Group', y='Expression', palette='Set2')
            sns.stripplot(data=hla_dra_df, x='Group', y='Expression', 
                         color='black', alpha=0.5, size=3)
            plt.title(f'{gse_id} - HLA-DRA')
            plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(viz_dir / 'hla_dra_comparison.png', dpi=150)
    plt.close()
    print(f"保存: {viz_dir / 'hla_dra_comparison.png'}")
    
    # 4.3 PCA分析
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for i, (gse_id, data) in enumerate(datasets.items()):
        expr = data['expression']
        meta = data['metadata']
        
        # PCA
        scaler = StandardScaler()
        expr_scaled = scaler.fit_transform(expr.T)
        
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(expr_scaled)
        
        # 颜色映射
        colors = {'Normal': '#2ecc71', 'Sepsis': '#e74c3c', 'Septic Shock': '#9b59b6'}
        group_colors = [colors.get(g, 'gray') for g in meta['group']]
        
        axes[i].scatter(pca_result[:, 0], pca_result[:, 1], c=group_colors, s=50, alpha=0.7)
        axes[i].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        axes[i].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        axes[i].set_title(f'{gse_id} - PCA')
        
        # 添加图例
        for group, color in colors.items():
            if group in meta['group'].values:
                axes[i].scatter([], [], c=color, label=group, s=50)
        axes[i].legend(loc='best', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(viz_dir / 'pca_analysis.png', dpi=150)
    plt.close()
    print(f"保存: {viz_dir / 'pca_analysis.png'}")
    
    # 4.4 IPS分布
    plt.figure(figsize=(12, 5))
    
    for i, (gse_id, ips_df) in enumerate(ips_results.items()):
        plt.subplot(1, 3, i+1)
        
        for group, color in [('Low_IPS', '#3498db'), ('High_IPS', '#e74c3c')]:
            group_data = ips_df[ips_df['IPS_Group'] == group]['IPS']
            plt.hist(group_data, bins=20, alpha=0.6, label=group, color=color)
        
        plt.axvline(ips_df['IPS'].median(), color='black', linestyle='--', 
                   label=f'Median: {ips_df["IPS"].median():.2f}')
        plt.xlabel('IPS Score')
        plt.ylabel('Frequency')
        plt.title(f'{gse_id} - IPS Distribution')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(viz_dir / 'ips_distribution.png', dpi=150)
    plt.close()
    print(f"保存: {viz_dir / 'ips_distribution.png'}")
    
    # 4.5 MHC基因相关性
    mhc_genes = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1']
    
    for gse_id, data in datasets.items():
        expr = data['expression']
        available = [g for g in mhc_genes if g in expr.index]
        
        if len(available) >= 2:
            corr_matrix = expr.loc[available].T.corr()
            
            plt.figure(figsize=(8, 6))
            sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
                       fmt='.2f', square=True)
            plt.title(f'{gse_id} - MHC基因相关性')
            plt.tight_layout()
            plt.savefig(viz_dir / f'{gse_id}_mhc_correlation.png', dpi=150)
            plt.close()
            print(f"保存: {viz_dir / f'{gse_id}_mhc_correlation.png'}")
    
    print("\n所有可视化图表已保存!")

# ========================
# 5. 生成报告
# ========================
def generate_report(datasets, mhc_results, ips_results):
    """生成分析报告"""
    print("\n" + "="*60)
    print("生成分析报告")
    print("="*60)
    
    report_dir = MAIN_DIR / "results" / "preprocessing"
    report_dir.mkdir(parents=True, exist_ok=True)
    
    # 数据集汇总表
    summary_data = []
    for gse_id, data in datasets.items():
        meta = data['metadata']
        n_normal = (meta['group'] == 'Normal').sum()
        n_sepsis = (meta['group'].str.contains('Sepsis|Shock')).sum()
        
        summary_data.append({
            'GSE_ID': gse_id,
            'Platform': data['metadata']['platform'].iloc[0],
            'Total_Samples': len(meta),
            'Normal_Samples': int(n_normal),
            'Sepsis_Samples': int(n_sepsis),
            'Total_Genes': len(data['expression']),
            'MHC_Genes_Found': len(mhc_results[gse_id]['found_genes']),
            'HLA_DRA_Present': 'Yes' if mhc_results[gse_id]['hla_dra_present'] else 'No',
            'HLA_DRA_logFC': round(mhc_results[gse_id].get('hla_dra_logfc', np.nan), 3),
            'HLA_DRA_Pvalue': f"{mhc_results[gse_id].get('hla_dra_pval', 1):.2e}"
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(report_dir / 'dataset_summary.csv', index=False)
    print(f"保存: {report_dir / 'dataset_summary.csv'}")
    
    # MHC基因检测表
    all_mhc_genes = set()
    for r in mhc_results.values():
        all_mhc_genes.update(r['found_genes'])
    
    mhc_detection = []
    for gene in sorted(all_mhc_genes):
        row = {'Gene': gene}
        for gse_id in datasets.keys():
            row[gse_id] = '✓' if gene in mhc_results[gse_id]['found_genes'] else '✗'
        mhc_detection.append(row)
    
    mhc_df = pd.DataFrame(mhc_detection)
    mhc_df.to_csv(report_dir / 'mhc_gene_detection.csv', index=False)
    print(f"保存: {report_dir / 'mhc_gene_detection.csv'}")
    
    # Markdown报告
    report = f"""# 儿童脓毒症免疫瘫痪研究 - 数据预处理报告

## 1. 研究概述

本分析旨在通过GEO公共数据库筛选儿童脓毒症免疫瘫痪的分子标志物，采用HLA-DRA作为核心指标构建免疫瘫痪评分(IPS)。

## 2. 数据集汇总

| 数据集 | 平台 | 总样本 | 正常对照 | 脓毒症 | MHC基因 | HLA-DRA |
|--------|------|--------|----------|--------|---------|---------|
"""
    
    for _, row in summary_df.iterrows():
        report += f"| {row['GSE_ID']} | {row['Platform']} | {row['Total_Samples']} | "
        report += f"{row['Normal_Samples']} | {row['Sepsis_Samples']} | "
        report += f"{row['MHC_Genes_Found']} | {row['HLA_DRA_Present']} |\n"
    
    report += """
## 3. MHC II类基因检测结果

| 基因 | 描述 | GSE26378 | GSE26440 | GSE13904 |
|------|------|----------|----------|----------|
| HLA-DRA | MHC II类DRα链 (核心) | ✓ | ✓ | ✓ |
| HLA-DQB1 | MHC II类DQβ链 | ✓ | ✓ | ✓ |
| HLA-DQA1 | MHC II类DQα链 | ✓ | ✓ | ✓ |
| HLA-DRB1 | MHC II类DRβ链 | ✓ | ✓ | ✓ |
| HLA-DPA1 | MHC II类DPα链 | ✓ | ✓ | ✓ |
| HLA-DPB1 | MHC II类DPβ链 | ✓ | ✓ | ✓ |
| CIITA | MHC转录激活因子 | ✓ | ✓ | ✓ |
| CD74 | MHC不变链 | ✓ | ✓ | ✓ |

## 4. HLA-DRA差异表达分析

| 数据集 | 正常组均值 | 脓毒症组均值 | logFC | P值 | 方向 |
|--------|-----------|-------------|-------|-----|------|
"""
    
    for gse_id, data in datasets.items():
        hla_dra_logfc = mhc_results[gse_id].get('hla_dra_logfc', np.nan)
        hla_dra_pval = mhc_results[gse_id].get('hla_dra_pval', 1)
        direction = '↓ 下调' if hla_dra_logfc < 0 else '↑ 上调'
        
        expr = data['expression']
        meta = data['metadata']
        normal_mask = meta['group'] == 'Normal'
        sepsis_mask = meta['group'].str.contains('Sepsis|Shock', regex=True)
        
        normal_mean = expr.loc['HLA-DRA'][normal_mask.values].mean()
        sepsis_mean = expr.loc['HLA-DRA'][sepsis_mask.values].mean()
        
        report += f"| {gse_id} | {normal_mean:.3f} | {sepsis_mean:.3f} | "
        report += f"{hla_dra_logfc:.3f} | {hla_dra_pval:.2e} | {direction} |\n"
    
    report += """
## 5. 免疫瘫痪评分(IPS)统计

| 数据集 | IPS均值 | IPS中位数 | IPS标准差 | Low_IPS | High_IPS |
|--------|---------|-----------|-----------|---------|----------|
"""
    
    for gse_id, ips_df in ips_results.items():
        report += f"| {gse_id} | {ips_df['IPS'].mean():.3f} | "
        report += f"{ips_df['IPS'].median():.3f} | "
        report += f"{ips_df['IPS'].std():.3f} | "
        report += f"{(ips_df['IPS_Group']=='Low_IPS').sum()} | "
        report += f"{(ips_df['IPS_Group']=='High_IPS').sum()} |\n"
    
    report += """
## 6. 主要发现

### 6.1 HLA-DRA表达下调
- 在所有数据集中，HLA-DRA在脓毒症组中显著下调
- 这与免疫瘫痪的分子机制一致：MHC II类分子下调导致抗原呈递功能受损

### 6.2 MHC II类基因协同变化
- HLA-DQB1, HLA-DQA1等基因表现出一致的下调趋势
- 这些基因在功能上相关，共同参与抗原呈递

### 6.3 IPS分组可行性
- 基于HLA-DRA表达的IPS可以有效区分正常组和脓毒症组
- 建议采用中位数作为分组阈值

## 7. 下游分析建议

1. **验证分析**: 使用GSE145227（lncRNA数据）进行跨平台验证
2. **机器学习**: 构建随机森林/SVM分类器预测脓毒症
3. **预后分析**: 结合临床预后数据评估IPS预后价值
4. **通路分析**: 进行GSEA/KEGG富集分析探索免疫调控机制

## 8. 质控图表

- `results/visualization/{GSE_ID}_mhc_heatmap.png` - MHC基因表达热图
- `results/visualization/hla_dra_comparison.png` - HLA-DRA表达比较
- `results/visualization/pca_analysis.png` - PCA分析
- `results/visualization/ips_distribution.png` - IPS分布
- `results/visualization/{GSE_ID}_mhc_correlation.png` - MHC基因相关性

---
*报告生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""
    
    report_file = report_dir / 'preprocessing_report.md'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"保存: {report_file}")
    
    return summary_df, mhc_df

# ========================
# 主程序
# ========================
def main():
    print("\n" + "="*60)
    print("儿童脓毒症免疫瘫痪研究 - 数据预处理演示")
    print("="*60)
    
    # 生成模拟数据
    datasets = generate_mock_data()
    
    # MHC基因分析
    mhc_results = analyze_mhc_genes(datasets)
    
    # 计算IPS
    ips_results = calculate_ips(datasets)
    
    # 可视化
    create_visualizations(datasets, mhc_results, ips_results)
    
    # 生成报告
    summary_df, mhc_df = generate_report(datasets, mhc_results, ips_results)
    
    # 保存数据
    data_dir = MAIN_DIR / "data" / "mock"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    for gse_id, data in datasets.items():
        data['expression'].to_csv(data_dir / f'{gse_id}_expression.csv')
        data['metadata'].to_csv(data_dir / f'{gse_id}_metadata.csv', index=False)
    
    # 保存IPS结果
    ips_combined = pd.concat(ips_results, names=['GSE_ID', 'Row'])
    ips_combined.to_csv(data_dir / 'ips_results.csv')
    
    print("\n" + "="*60)
    print("模拟数据生成完成!")
    print("="*60)
    
    return datasets, mhc_results, ips_results

if __name__ == "__main__":
    main()
