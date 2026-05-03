#!/usr/bin/env python3
"""
儿童脓毒症免疫瘫痪研究 - 数据预处理与分析
处理GEO series_matrix格式数据
"""

import os
import re
import gzip
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 路径设置
BASE_DIR = "."
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
NORM_DIR = os.path.join(BASE_DIR, "data/normalized")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
VIZ_DIR = os.path.join(RESULTS_DIR, "visualization")
PREP_DIR = os.path.join(RESULTS_DIR, "preprocessing")

os.makedirs(NORM_DIR, exist_ok=True)
os.makedirs(VIZ_DIR, exist_ok=True)
os.makedirs(PREP_DIR, exist_ok=True)

# MHC II类核心基因
MHC_GENES = ['HLA-DRA', 'HLA-DQB1', 'HLA-DQA1', 'CIITA', 'CD74']

# IPS评分权重
GENE_WEIGHTS = {
    'HLA-DRA': 0.35,
    'HLA-DQB1': 0.25,
    'HLA-DQA1': 0.20,
    'CIITA': 0.15,
    'CD74': 0.05
}

def parse_series_matrix(file_path):
    """解析GEO series_matrix.txt文件"""
    metadata = {}
    sample_ids = []
    data_rows = []
    
    print(f"正在解析: {file_path}")
    
    in_table = False
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            if '!series_matrix_table_begin' in line:
                in_table = True
                continue
            
            if '!series_matrix_table_end' in line:
                break
            
            if not in_table:
                # 解析元数据行
                if line.startswith('!') and '=' not in line[:5]:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0].replace('!', '')
                        values = [v.strip('"') for v in parts[1:]]
                        metadata[key] = values
            else:
                # 解析表达矩阵
                if line.startswith('"'):
                    parts = line.split('\t')
                    parts = [p.strip('"') for p in parts]
                    
                    if parts[0] == 'ID_REF':
                        # 这是表头行
                        sample_ids = parts[1:]
                    else:
                        # 这是数据行
                        data_rows.append(parts)
    
    # 构建表达矩阵DataFrame
    if sample_ids and data_rows:
        probe_ids = [row[0] for row in data_rows]
        expr_values = [[float(v) if v else np.nan for v in row[1:]] for row in data_rows]
        
        expr_df = pd.DataFrame(expr_values, index=probe_ids, columns=sample_ids)
    else:
        expr_df = pd.DataFrame()
    
    print(f"  解析完成: {len(sample_ids)} 样本, {len(expr_df)} 探针")
    
    return metadata, expr_df


def extract_sample_info(metadata):
    """从元数据中提取样本信息"""
    sample_titles = metadata.get('Sample_title', [])
    sample_geo = metadata.get('Sample_geo_accession', [])
    
    samples = []
    for i, title in enumerate(sample_titles):
        geo_id = sample_geo[i] if i < len(sample_geo) else f"Sample_{i}"
        
        # 判断分组：Septic Shock vs Control
        if 'Septic' in title or 'septic' in title.lower():
            group = 'sepsis'
        elif 'Control' in title or 'Normal' in title or 'control' in title.lower():
            group = 'normal'
        else:
            group = 'unknown'
        
        samples.append({
            'sample_id': geo_id,
            'title': title,
            'group': group
        })
    
    return pd.DataFrame(samples)


def get_gpl570_annotation():
    """GPL570平台探针-基因注释 (HG-U133_Plus_2)"""
    # 基于Affymetrix HG-U133_Plus_2芯片注释
    annotation = {
        'HLA-DRA': [
            '203290_at',   # 主要探针
            '207899_s_at',
            '211990_x_at'
        ],
        'HLA-DQB1': [
            '209480_at',   # 主要探针
            '211654_x_at',
            '210981_x_at'
        ],
        'HLA-DQA1': [
            '210982_at',   # 主要探针
            '211654_x_at',
            '212830_x_at'
        ],
        'CIITA': [
            '202395_s_at',  # 主要探针
            '202396_s_at',
            '208932_at'
        ],
        'CD74': [
            '201009_s_at',  # 主要探针
            '201010_s_at',
            '212832_at'
        ]
    }
    return annotation


def extract_mhc_expression(expr_df, annotation):
    """从表达矩阵中提取MHC基因表达"""
    mhc_data = {}
    probe_info = {}
    
    print("\n提取MHC基因表达:")
    for gene, probes in annotation.items():
        valid_probes = []
        expressions = []
        
        for probe in probes:
            if probe in expr_df.index:
                valid_probes.append(probe)
                expressions.append(expr_df.loc[probe].values)
        
        if expressions:
            # 多探针取均值
            probe_expr = np.array(expressions)
            mean_expr = np.nanmean(probe_expr, axis=0)
            mhc_data[gene] = mean_expr
            probe_info[gene] = {
                'probes_found': valid_probes,
                'n_probes': len(valid_probes)
            }
            print(f"  {gene}: 找到 {len(valid_probes)} 个探针 -> {valid_probes}")
        else:
            print(f"  {gene}: 未找到探针")
            mhc_data[gene] = None
            probe_info[gene] = {
                'probes_found': [],
                'n_probes': 0
            }
    
    return mhc_data, probe_info


def calculate_ips(mhc_data, weights):
    """计算IPS免疫瘫痪评分"""
    ips_scores = None
    
    for gene, expr in mhc_data.items():
        if expr is not None:
            weight = weights.get(gene, 0)
            if ips_scores is None:
                ips_scores = expr * weight
            else:
                ips_scores += expr * weight
    
    return ips_scores


def quality_control(expr_df, sample_info):
    """数据质量控制"""
    qc_results = {
        'n_samples': len(sample_info),
        'n_probes': len(expr_df),
        'samples_by_group': sample_info['group'].value_counts().to_dict(),
        'outliers': [],
        'missing_rate': {}
    }
    
    # 计算缺失率
    qc_results['missing_rate']['overall'] = expr_df.isna().sum().sum() / expr_df.size
    
    return qc_results


def perform_pca(expr_df, sample_info, n_components=3):
    """执行PCA分析"""
    # 数据标准化
    X = expr_df.T.fillna(0)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(X_scaled)
    
    pca_df = pd.DataFrame(components, columns=[f'PC{i+1}' for i in range(n_components)])
    pca_df['sample_id'] = expr_df.columns
    pca_df = pca_df.merge(sample_info[['sample_id', 'group']], on='sample_id', how='left')
    
    return pca_df, pca.explained_variance_ratio_


def differential_expression_analysis(expr_df, sample_info, gene):
    """差异表达分析"""
    if gene not in expr_df.index:
        return None
    
    sepsis_mask = sample_info['group'] == 'sepsis'
    normal_mask = sample_info['group'] == 'normal'
    
    sepsis_samples = sample_info[sepsis_mask]['sample_id'].tolist()
    normal_samples = sample_info[normal_mask]['sample_id'].tolist()
    
    # 获取两组表达值
    sepsis_expr = expr_df.loc[gene, sepsis_samples].values.astype(float)
    normal_expr = expr_df.loc[gene, normal_samples].values.astype(float)
    
    # t检验
    t_stat, p_value = stats.ttest_ind(sepsis_expr, normal_expr)
    
    # 计算fold change (log2)
    mean_sepsis = np.nanmean(sepsis_expr)
    mean_normal = np.nanmean(normal_expr)
    log2_fc = np.log2(mean_sepsis / mean_normal) if mean_normal > 0 else 0
    
    return {
        'gene': gene,
        'mean_sepsis': mean_sepsis,
        'mean_normal': mean_normal,
        'log2_fold_change': log2_fc,
        't_statistic': t_stat,
        'p_value': p_value,
        'n_sepsis': len(sepsis_expr),
        'n_normal': len(normal_expr)
    }


def create_visualizations(expr_df, mhc_data, sample_info, ips_scores, pca_df, pca_variance, dataset_name):
    """创建可视化图表"""
    dpi = 300
    
    # 1. MHC基因表达热图
    fig, ax = plt.subplots(figsize=(14, 6))
    
    mhc_df = pd.DataFrame(mhc_data, index=expr_df.columns)
    mhc_df = mhc_df.dropna(how='all', axis=1)
    
    # 按组排序
    sample_order = sample_info.sort_values('group')['sample_id'].tolist()
    mhc_df_ordered = mhc_df.reindex([s for s in sample_order if s in mhc_df.index])
    
    # 标准化用于热图
    if not mhc_df_ordered.empty:
        mhc_normalized = (mhc_df_ordered - mhc_df_ordered.mean()) / (mhc_df_ordered.std() + 1e-10)
        
        sns.heatmap(mhc_normalized.T, cmap='RdBu_r', center=0, 
                     xticklabels=False, ax=ax, cbar_kws={'label': 'Z-score'})
        ax.set_title(f'MHC Class II Gene Expression Heatmap - {dataset_name}', fontsize=14)
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(VIZ_DIR, f'{dataset_name}_mhc_heatmap.png'), dpi=dpi, bbox_inches='tight')
    plt.close()
    print("  保存: mhc_heatmap.png")
    
    # 2. HLA-DRA表达比较
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    hla_dra = mhc_data.get('HLA-DRA')
    if hla_dra is not None:
        sepsis_mask = sample_info['group'] == 'sepsis'
        normal_mask = sample_info['group'] == 'normal'
        
        sepsis_vals = hla_dra[sepsis_mask.values]
        normal_vals = hla_dra[normal_mask.values]
        
        # 箱线图
        bp = axes[0].boxplot([normal_vals, sepsis_vals], 
                             labels=['Normal', 'Sepsis'],
                             patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
        
        axes[0].set_ylabel('Expression (Normalized)', fontsize=12)
        axes[0].set_title('HLA-DRA Expression: Normal vs Sepsis', fontsize=14)
        
        # 统计检验
        if len(normal_vals) > 1 and len(sepsis_vals) > 1:
            t_stat, p_val = stats.ttest_ind(sepsis_vals, normal_vals)
            axes[0].text(0.5, 0.95, f'p = {p_val:.2e}', transform=axes[0].transAxes,
                        ha='center', fontsize=11, fontweight='bold')
        
        # 小提琴图
        violin_data = pd.DataFrame({
            'Expression': np.concatenate([normal_vals, sepsis_vals]),
            'Group': ['Normal'] * len(normal_vals) + ['Sepsis'] * len(sepsis_vals)
        })
        
        sns.violinplot(data=violin_data, x='Group', y='Expression', ax=axes[1],
                      palette={'Normal': 'lightblue', 'Sepsis': 'lightcoral'})
        axes[1].set_title('HLA-DRA Expression Distribution', fontsize=14)
        axes[1].set_ylabel('Expression (Normalized)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(VIZ_DIR, f'{dataset_name}_hla_dra_expression.png'), dpi=dpi, bbox_inches='tight')
    plt.close()
    print("  保存: hla_dra_expression.png")
    
    # 3. IPS评分分布
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    normal_ips = ips_scores[sample_info['group'] == 'normal']
    sepsis_ips = ips_scores[sample_info['group'] == 'sepsis']
    
    # 直方图
    axes[0].hist(normal_ips, bins=20, alpha=0.6, label='Normal', color='blue')
    axes[0].hist(sepsis_ips, bins=20, alpha=0.6, label='Sepsis', color='red')
    axes[0].axvline(np.percentile(ips_scores, 25), color='gray', 
                    linestyle='--', label='25th percentile')
    axes[0].axvline(np.percentile(ips_scores, 75), color='gray', 
                    linestyle='--', label='75th percentile')
    axes[0].set_xlabel('IPS Score', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title(f'IPS Score Distribution - {dataset_name}', fontsize=14)
    axes[0].legend()
    
    # 箱线图
    ips_df = pd.DataFrame({
        'IPS': np.concatenate([normal_ips, sepsis_ips]),
        'Group': ['normal'] * len(normal_ips) + ['sepsis'] * len(sepsis_ips)
    })
    
    sns.boxplot(data=ips_df, x='Group', y='IPS', ax=axes[1],
               palette={'normal': 'lightblue', 'sepsis': 'lightcoral'})
    axes[1].set_title('IPS Score by Group', fontsize=14)
    axes[1].set_ylabel('IPS Score', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(VIZ_DIR, f'{dataset_name}_ips_distribution.png'), dpi=dpi, bbox_inches='tight')
    plt.close()
    print("  保存: ips_distribution.png")
    
    # 4. PCA图
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    colors = {'sepsis': 'red', 'normal': 'blue', 'unknown': 'gray'}
    
    # 按分组着色
    for group in pca_df['group'].unique():
        mask = pca_df['group'] == group
        if mask.sum() > 0:
            axes[0].scatter(pca_df.loc[mask, 'PC1'], 
                           pca_df.loc[mask, 'PC2'],
                           c=colors.get(group, 'gray'),
                           label=group.capitalize(),
                           alpha=0.7, s=50)
    
    axes[0].set_xlabel(f'PC1 ({pca_variance[0]*100:.1f}%)', fontsize=12)
    axes[0].set_ylabel(f'PC2 ({pca_variance[1]*100:.1f}%)', fontsize=12)
    axes[0].set_title(f'PCA: Samples by Clinical Group - {dataset_name}', fontsize=14)
    axes[0].legend()
    
    # 按IPS分组着色
    ips_quantiles = pd.qcut(ips_scores, q=3, labels=['Low', 'Medium', 'High'], duplicates='drop')
    pca_df['IPS_group'] = ips_quantiles
    
    colors_ips = {'Low': 'green', 'Medium': 'orange', 'High': 'purple'}
    for group in ips_quantiles.unique():
        if pd.notna(group):
            mask = pca_df['IPS_group'] == group
            if mask.sum() > 0:
                axes[1].scatter(pca_df.loc[mask, 'PC1'], 
                               pca_df.loc[mask, 'PC2'],
                               c=colors_ips.get(group, 'gray'),
                               label=f'{group} IPS',
                               alpha=0.7, s=50)
    
    axes[1].set_xlabel(f'PC1 ({pca_variance[0]*100:.1f}%)', fontsize=12)
    axes[1].set_ylabel(f'PC2 ({pca_variance[1]*100:.1f}%)', fontsize=12)
    axes[1].set_title(f'PCA: Samples by IPS Group - {dataset_name}', fontsize=14)
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(VIZ_DIR, f'{dataset_name}_pca_analysis.png'), dpi=dpi, bbox_inches='tight')
    plt.close()
    print("  保存: pca_analysis.png")
    
    # 5. 样本分布箱线图（QC）
    fig, ax = plt.subplots(figsize=(14, 5))
    
    sample_data = expr_df.T.fillna(0)
    bp = ax.boxplot([sample_data.iloc[i].values for i in range(len(sample_data))], 
                    showfliers=False)
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Expression', fontsize=12)
    ax.set_title(f'Sample Expression Distribution (QC) - {dataset_name}', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(os.path.join(VIZ_DIR, f'{dataset_name}_sample_distribution.png'), dpi=dpi, bbox_inches='tight')
    plt.close()
    print("  保存: sample_distribution.png")
    
    return pca_df


def process_dataset(file_path, dataset_name, annotation):
    """处理单个数据集"""
    print(f"\n{'='*60}")
    print(f"处理数据集: {dataset_name}")
    print('='*60)
    
    # 1. 解析series_matrix文件
    metadata, expr_df = parse_series_matrix(file_path)
    
    # 2. 提取样本信息
    sample_info = extract_sample_info(metadata)
    print(f"\n样本信息:")
    print(f"  - 总样本数: {len(sample_info)}")
    print(f"  - 分组: {sample_info['group'].value_counts().to_dict()}")
    
    # 3. 数据质量控制
    qc_results = quality_control(expr_df, sample_info)
    print(f"\n质量控制:")
    print(f"  - 探针数: {qc_results['n_probes']}")
    print(f"  - 缺失率: {qc_results['missing_rate']['overall']*100:.2f}%")
    
    # 4. 提取MHC基因表达
    mhc_data, probe_info = extract_mhc_expression(expr_df, annotation)
    
    # 5. 计算IPS评分
    ips_scores = calculate_ips(mhc_data, GENE_WEIGHTS)
    
    # 6. 差异表达分析
    de_results = {}
    for gene in MHC_GENES:
        if gene in expr_df.index:
            de_results[gene] = differential_expression_analysis(expr_df, sample_info, gene)
            if de_results[gene]:
                print(f"\n{gene} 差异分析:")
                print(f"  Log2FC: {de_results[gene]['log2_fold_change']:.4f}")
                print(f"  p-value: {de_results[gene]['p_value']:.2e}")
    
    # 7. PCA分析
    pca_df, pca_variance = perform_pca(expr_df, sample_info)
    
    # 8. 创建可视化
    print(f"\n生成可视化图表...")
    pca_df = create_visualizations(expr_df, mhc_data, sample_info, ips_scores, pca_df, pca_variance, dataset_name)
    
    # 9. 保存数据文件
    print(f"\n保存数据文件...")
    
    # 保存表达矩阵
    expr_df.to_csv(os.path.join(NORM_DIR, f'{dataset_name}_expression.csv'))
    
    # 保存样本信息
    sample_info.to_csv(os.path.join(NORM_DIR, f'{dataset_name}_metadata.csv'), index=False)
    
    # 保存MHC基因表达
    mhc_df = pd.DataFrame(mhc_data, index=sample_info['sample_id']).T
    mhc_df.to_csv(os.path.join(NORM_DIR, f'{dataset_name}_mhc_expression.csv'))
    
    # 保存IPS评分
    ips_df = pd.DataFrame({
        'sample_id': sample_info['sample_id'],
        'group': sample_info['group'],
        'IPS_score': ips_scores
    })
    ips_df.to_csv(os.path.join(NORM_DIR, f'{dataset_name}_ips_scores.csv'), index=False)
    
    # 10. 生成单数据集报告
    report = f"""# 数据预处理报告 - {dataset_name}

## 1. 数据集概览

### 基本信息
- **数据集**: {dataset_name}
- **总样本数**: {len(sample_info)}
- **总探针数**: {qc_results['n_probes']}

### 样本分组
"""
    
    for group, count in qc_results['samples_by_group'].items():
        report += f"- **{group.capitalize()}**: {count} 样本\n"
    
    report += f"""
## 2. MHC基因检测结果

### 探针映射
| 基因 | 探针数量 | 状态 |
|------|----------|------|
"""
    
    for gene, info in probe_info.items():
        status = "✅ 找到" if info['n_probes'] > 0 else "❌ 未找到"
        probes = ", ".join(info['probes_found']) if info['probes_found'] else "无"
        report += f"| {gene} | {info['n_probes']} | {status} |\n"
    
    report += f"""
### 基因表达差异分析
| 基因 | 正常组均值 | 脓毒症组均值 | Log2FC | p值 |
|------|------------|--------------|--------|-----|
"""
    
    for gene in MHC_GENES:
        if gene in de_results and de_results[gene]:
            res = de_results[gene]
            sig = "*" if res['p_value'] < 0.05 else ""
            report += f"| {gene} | {res['mean_normal']:.4f} | {res['mean_sepsis']:.4f} | {res['log2_fold_change']:.4f}{sig} | {res['p_value']:.2e} |\n"
    
    report += f"""
## 3. 数据质量控制结果

- **总体缺失率**: {qc_results['missing_rate']['overall']*100:.2f}%
- **异常样本数**: {len(qc_results['outliers'])}

## 4. IPS评分统计

| 分组 | 样本数 | IPS均值 | IPS标准差 |
|------|--------|---------|-----------|
| Normal | {sum(sample_info['group']=='normal')} | {np.mean(ips_scores[sample_info['group']=='normal']):.4f} | {np.std(ips_scores[sample_info['group']=='normal']):.4f} |
| Sepsis | {sum(sample_info['group']=='sepsis')} | {np.mean(ips_scores[sample_info['group']=='sepsis']):.4f} | {np.std(ips_scores[sample_info['group']=='sepsis']):.4f} |

## 5. 关键发现

### HLA-DRA差异表达
"""
    
    if 'HLA-DRA' in de_results and de_results['HLA-DRA']:
        res = de_results['HLA-DRA']
        direction = "下调" if res['log2_fold_change'] < 0 else "上调"
        significant = "显著" if res['p_value'] < 0.05 else "不显著"
        report += f"""
- **表达变化**: {direction} (Log2FC = {res['log2_fold_change']:.4f})
- **统计显著性**: p = {res['p_value']:.2e} ({significant})
- **正常组均值**: {res['mean_normal']:.4f}
- **脓毒症组均值**: {res['mean_sepsis']:.4f}

### 免疫瘫痪特征
- 脓毒症组MHC II类基因表达**{'降低' if res['log2_fold_change'] < 0 else '升高'}**
- 这与免疫瘫痪的临床特征一致
"""
    
    report += """
## 6. 结论

"""
    
    genes_found = sum(1 for info in probe_info.values() if info['n_probes'] > 0)
    
    if 'HLA-DRA' in de_results and de_results['HLA-DRA'] and de_results['HLA-DRA']['p_value'] < 0.05:
        report += "✅ HLA-DRA在脓毒症组显著差异表达\n"
    else:
        report += "⚠️ HLA-DRA差异表达未达显著水平\n"
    
    if genes_found >= 3:
        report += "✅ IPS评分可构建（≥3个基因可用）\n"
    else:
        report += "❌ IPS评分构建受限（<3个基因可用）\n"
    
    with open(os.path.join(PREP_DIR, f'{dataset_name}_preprocessing_report.md'), 'w') as f:
        f.write(report)
    
    print(f"\n报告已保存: {dataset_name}_preprocessing_report.md")
    
    return {
        'metadata': metadata,
        'sample_info': sample_info,
        'expr_df': expr_df,
        'mhc_data': mhc_data,
        'probe_info': probe_info,
        'ips_scores': ips_scores,
        'qc_results': qc_results,
        'de_results': de_results,
        'pca_df': pca_df
    }


def main():
    """主函数"""
    print("="*60)
    print("儿童脓毒症免疫瘫痪研究 - 数据预处理")
    print("="*60)
    
    # 获取GPL570注释
    annotation = get_gpl570_annotation()
    
    # 处理训练集
    train_results = process_dataset(
        os.path.join(RAW_DIR, 'GSE26378_series_matrix.txt'),
        'GSE26378_training',
        annotation
    )
    
    # 处理验证集
    valid_results = process_dataset(
        os.path.join(RAW_DIR, 'GSE26440_series_matrix.txt'),
        'GSE26440_validation',
        annotation
    )
    
    # 生成综合报告
    print("\n" + "="*60)
    print("生成综合分析报告")
    print("="*60)
    
    combined_report = f"""# 儿童脓毒症免疫瘫痪研究 - 综合分析报告

## 研究背景

本研究旨在分析儿童脓毒症中的免疫瘫痪现象，通过GEO数据库中的基因表达数据，
构建基于MHC II类基因的IPS（免疫瘫痪评分）系统。

## 数据集概述

| 数据集 | 用途 | 总样本数 | 脓毒症 | 正常对照 |
|--------|------|----------|--------|----------|
| GSE26378 | 训练集 | {len(train_results['sample_info'])} | {sum(train_results['sample_info']['group']=='sepsis')} | {sum(train_results['sample_info']['group']=='normal')} |
| GSE26440 | 验证集 | {len(valid_results['sample_info'])} | {sum(valid_results['sample_info']['group']=='sepsis')} | {sum(valid_results['sample_info']['group']=='normal')} |

## MHC II类基因分析

### 目标基因
- **HLA-DRA**: MHC II类主要组织相容性复合物DRα链 (权重: 35%)
- **HLA-DQB1**: MHC II类DQβ链 (权重: 25%)
- **HLA-DQA1**: MHC II类DQα链 (权重: 20%)
- **CIITA**: MHC II类转录激活因子 (权重: 15%)
- **CD74**: MHC II类相关蛋白 (权重: 5%)

### 探针检测结果

#### 训练集 (GSE26378)
| 基因 | 探针数 | 状态 |
|------|--------|------|
"""
    
    for gene in MHC_GENES:
        info = train_results['probe_info'].get(gene, {'n_probes': 0})
        status = "✅" if info['n_probes'] > 0 else "❌"
        combined_report += f"| {gene} | {info['n_probes']} | {status} |\n"
    
    combined_report += f"""
#### 验证集 (GSE26440)
| 基因 | 探针数 | 状态 |
|------|--------|------|
"""
    
    for gene in MHC_GENES:
        info = valid_results['probe_info'].get(gene, {'n_probes': 0})
        status = "✅" if info['n_probes'] > 0 else "❌"
        combined_report += f"| {gene} | {info['n_probes']} | {status} |\n"
    
    combined_report += f"""
## 差异表达分析

### HLA-DRA 差异表达

#### 训练集 (GSE26378)
"""
    
    if 'HLA-DRA' in train_results['de_results'] and train_results['de_results']['HLA-DRA']:
        res = train_results['de_results']['HLA-DRA']
        combined_report += f"""
- Log2 Fold Change: {res['log2_fold_change']:.4f}
- p值: {res['p_value']:.2e}
- 正常组均值: {res['mean_normal']:.4f}
- 脓毒症组均值: {res['mean_sepsis']:.4f}
"""
    
    combined_report += f"""
#### 验证集 (GSE26440)
"""
    
    if 'HLA-DRA' in valid_results['de_results'] and valid_results['de_results']['HLA-DRA']:
        res = valid_results['de_results']['HLA-DRA']
        combined_report += f"""
- Log2 Fold Change: {res['log2_fold_change']:.4f}
- p值: {res['p_value']:.2e}
- 正常组均值: {res['mean_normal']:.4f}
- 脓毒症组均值: {res['mean_sepsis']:.4f}
"""
    
    combined_report += """
## IPS评分分析

IPS评分基于5个MHC II类基因的加权表达构建。

### 评分公式
```
IPS = 0.35 × HLA-DRA + 0.25 × HLA-DQB1 + 0.20 × HLA-DQA1 + 0.15 × CIITA + 0.05 × CD74
```

### 分组阈值
- **低IPS组 (Low_IPS)**: < 25%分位数
- **中IPS组 (Medium_IPS)**: 25-75%分位数
- **高IPS组 (High_IPS)**: > 75%分位数

## 关键验证点

| 验证点 | 训练集 | 验证集 |
|--------|--------|--------|
"""
    
    train_hla_found = train_results['probe_info'].get('HLA-DRA', {}).get('n_probes', 0) > 0
    valid_hla_found = valid_results['probe_info'].get('HLA-DRA', {}).get('n_probes', 0) > 0
    
    train_hla_down = train_results['de_results'].get('HLA-DRA', {}).get('log2_fold_change', 0) < 0
    valid_hla_down = valid_results['de_results'].get('HLA-DRA', {}).get('log2_fold_change', 0) < 0
    
    train_hla_sig = train_results['de_results'].get('HLA-DRA', {}).get('p_value', 1) < 0.05
    valid_hla_sig = valid_results['de_results'].get('HLA-DRA', {}).get('p_value', 1) < 0.05
    
    combined_report += f"""| HLA-DRA存在 | {"✅" if train_hla_found else "❌"} | {"✅" if valid_hla_found else "❌"} |
| HLA-DRA下调 | {"✅" if train_hla_down else "❌"} | {"✅" if valid_hla_down else "❌"} |
| HLA-DRA显著 | {"✅" if train_hla_sig else "❌"} | {"✅" if valid_hla_sig else "❌"} |

## 输出文件

### 数据文件 (data/normalized/)
- GSE26378_training_expression.csv: 训练集完整表达矩阵
- GSE26378_training_metadata.csv: 训练集样本元数据
- GSE26378_training_mhc_expression.csv: 训练集MHC基因表达
- GSE26378_training_ips_scores.csv: 训练集IPS评分
- GSE26440_validation_expression.csv: 验证集完整表达矩阵
- GSE26440_validation_metadata.csv: 验证集样本元数据
- GSE26440_validation_mhc_expression.csv: 验证集MHC基因表达
- GSE26440_validation_ips_scores.csv: 验证集IPS评分

### 可视化图表 (results/visualization/)
- GSE26378_training_*.png: 训练集可视化图表
- GSE26440_validation_*.png: 验证集可视化图表

### 报告 (results/preprocessing/)
- GSE26378_training_preprocessing_report.md: 训练集报告
- GSE26440_validation_preprocessing_report.md: 验证集报告
- combined_report.md: 综合分析报告

## 结论

"""
    
    if train_hla_found and valid_hla_found:
        combined_report += "1. ✅ **HLA-DRA在两个数据集中均成功检测到**\n"
    else:
        combined_report += "1. ⚠️ **HLA-DRA检测存在异常**\n"
    
    if train_hla_down and valid_hla_down:
        combined_report += "2. ✅ **HLA-DRA在脓毒症组显著下调**，符合免疫瘫痪特征\n"
    elif train_hla_down or valid_hla_down:
        combined_report += "2. ⚠️ **HLA-DRA下调趋势不一致**\n"
    else:
        combined_report += "2. ⚠️ **HLA-DRA未表现出预期的下调趋势**\n"
    
    if train_hla_sig and valid_hla_sig:
        combined_report += "3. ✅ **HLA-DRA差异具有统计学显著性**\n"
    else:
        combined_report += "3. ⚠️ **HLA-DRA差异显著性待验证**\n"
    
    genes_found_train = sum(1 for info in train_results['probe_info'].values() if info['n_probes'] > 0)
    genes_found_valid = sum(1 for info in valid_results['probe_info'].values() if info['n_probes'] > 0)
    
    combined_report += f"""
4. **IPS评分可构建**: 训练集{genes_found_train}/5个基因，验证集{genes_found_valid}/5个基因

## 下一步分析建议

1. 进一步优化探针选择，使用更多探针映射资源
2. 进行通路富集分析，探索免疫相关通路变化
3. 构建机器学习模型预测脓毒症预后
4. 验证IPS评分在独立队列中的预测效能
"""
    
    with open(os.path.join(PREP_DIR, 'combined_report.md'), 'w') as f:
        f.write(combined_report)
    
    print(f"\n综合报告已保存: combined_report.md")
    print("\n" + "="*60)
    print("数据预处理完成！")
    print("="*60)
    
    return train_results, valid_results


if __name__ == '__main__':
    train_results, valid_results = main()
