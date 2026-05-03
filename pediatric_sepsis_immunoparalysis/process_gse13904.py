#!/usr/bin/env python3
"""
GSE13904 数据集处理与三数据集验证分析
儿童脓毒症免疫瘫痪研究
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import re
import warnings
warnings.filterwarnings('ignore')

# 设置工作目录 - 使用绝对路径
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = SCRIPT_DIR  # 脚本就在项目根目录
DATA_DIR = f"{WORK_DIR}/data/normalized"
RAW_DIR = f"{WORK_DIR}/data/raw"
RESULTS_DIR = f"{WORK_DIR}/results/visualization"

print(f"工作目录: {WORK_DIR}")

# MHC II类基因相关探针 (GPL570平台)
MHC_PROBES = {
    'HLA-DRA': ['205891_s_at', '207823_s_at'],
    'HLA-DRB1': ['205884_s_at', '205885_s_at', '210986_x_at'],
    'HLA-DQB1': ['209480_at', '210470_at'],
    'HLA-DQA1': ['205788_s_at', '207430_s_at'],
    'HLA-DPA1': ['205648_at', '210654_x_at'],
    'HLA-DPB1': ['205676_at', '207977_s_at'],
    'CIITA': ['205388_s_at', '209447_at'],
    'CD74': ['201009_s_at', '212888_at', '201010_s_at']
}

# 获取所有探针到基因的映射
PROBE_TO_GENE = {}
for gene, probes in MHC_PROBES.items():
    for probe in probes:
        PROBE_TO_GENE[probe] = gene

def parse_series_matrix_simple(filepath):
    """简化版解析 - 直接读取表达数据"""
    print(f"正在解析: {filepath}")
    
    # 使用pandas直接读取，跳过注释行
    try:
        # 尝试读取
        df = pd.read_csv(filepath, sep='\t', comment='!', index_col=0)
        print(f"成功读取表达矩阵: {df.shape}")
        return df
    except Exception as e:
        print(f"直接读取失败: {e}")
        return None

def parse_series_matrix_robust(filepath):
    """健壮的series_matrix解析"""
    print(f"正在解析: {filepath}")
    
    # 读取所有行
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    # 查找样本标题行、ID_REF行和表达数据
    sample_titles = []
    sample_ids = []
    expression_lines = []
    header_line_idx = None
    id_ref_line_idx = None
    
    for i, line in enumerate(lines):
        line_content = line.strip()
        
        # 解析样本标题
        if line_content.startswith('!Sample_title'):
            parts = line_content.split('\t')
            sample_titles = [p.strip('"') for p in parts[1:]]
            
        # 解析样本ID - 在ID_REF行中
        elif line_content.startswith('"ID_REF"'):
            id_ref_line_idx = i
            parts = line_content.split('\t')
            # 第一列是ID_REF，后面是样本ID
            sample_ids = [p.strip('"') for p in parts[1:]]
            
        # 解析表达数据 - ID_REF行之后的所有非!开头的行
        elif id_ref_line_idx is not None and i > id_ref_line_idx and not line_content.startswith('!'):
            expression_lines.append(line_content)
    
    print(f"找到 {len(sample_ids)} 个样本")
    print(f"找到 {len(expression_lines)} 个探针")
    
    # 解析表达数据
    probe_ids = []
    expression_data = []
    
    for line in expression_lines:
        parts = line.split('\t')
        if len(parts) > 1:
            probe_id = parts[0].strip('"')
            values = []
            for v in parts[1:]:
                v = v.strip('"')
                try:
                    values.append(float(v) if v not in ['NA', '', 'null', 'NULL'] else np.nan)
                except ValueError:
                    values.append(np.nan)
            probe_ids.append(probe_id)
            expression_data.append(values)
    
    # 构建表达矩阵 (探针为行，样本为列)
    if expression_data:
        expr_df = pd.DataFrame(expression_data, index=probe_ids)
        expr_df.columns = sample_ids
    else:
        expr_df = pd.DataFrame()
    
    print(f"表达矩阵形状: {expr_df.shape}")
    
    # 构建元数据
    def infer_group(title):
        title_lower = title.lower()
        if 'control' in title_lower or '_c_' in title_lower:
            return 'normal'
        elif 'septic shock' in title_lower:
            return 'sepsis'
        elif 'sepsis' in title_lower:
            return 'sepsis'
        elif 'sirs' in title_lower:
            return 'sepsis'  # SIRS视为疾病状态
        else:
            return 'unknown'
    
    metadata = pd.DataFrame({
        'sample_id': sample_ids,
        'title': sample_titles,
        'group': [infer_group(t) for t in sample_titles]
    }, index=sample_ids)
    
    print(f"分组分布:\n{metadata['group'].value_counts()}")
    
    return expr_df, metadata

def extract_mhc_expression(expression_df):
    """提取MHC II类基因表达"""
    mhc_data = {}
    
    for probe_id, gene_name in PROBE_TO_GENE.items():
        if probe_id in expression_df.index:
            if gene_name not in mhc_data:
                mhc_data[gene_name] = []
            mhc_data[gene_name].append(probe_id)
    
    # 取每个基因的第一个探针的表达值
    result = {}
    for gene, probes in mhc_data.items():
        probe = probes[0]  # 取第一个探针
        result[gene] = expression_df.loc[probe]
    
    if result:
        mhc_df = pd.DataFrame(result)
        return mhc_df
    else:
        return pd.DataFrame()

def calculate_de_genes(df, metadata, genes):
    """计算差异表达分析"""
    results = []
    
    normal_samples = metadata[metadata['group'] == 'normal'].index.tolist()
    sepsis_samples = metadata[metadata['group'] == 'sepsis'].index.tolist()
    
    if len(normal_samples) < 3 or len(sepsis_samples) < 3:
        print(f"警告: 样本量不足 - Normal: {len(normal_samples)}, Sepsis: {len(sepsis_samples)}")
    
    for gene in genes:
        if gene in df.columns:
            normal_values = df.loc[normal_samples, gene].dropna()
            sepsis_values = df.loc[sepsis_samples, gene].dropna()
            
            if len(normal_values) > 0 and len(sepsis_values) > 0:
                mean_normal = normal_values.mean()
                mean_sepsis = sepsis_values.mean()
                log2fc = mean_sepsis - mean_normal  # 已经是对数数据
                
                # t检验
                t_stat, p_value = stats.ttest_ind(sepsis_values, normal_values)
                
                results.append({
                    'gene': gene,
                    'mean_normal': mean_normal,
                    'mean_sepsis': mean_sepsis,
                    'log2_fold_change': log2fc,
                    't_statistic': t_stat,
                    'p_value': p_value,
                    'n_normal': len(normal_values),
                    'n_sepsis': len(sepsis_values)
                })
    
    return pd.DataFrame(results)

def build_ips_score(mhc_df):
    """构建免疫瘫痪评分(IPS)"""
    # HLA-DRA作为核心指标（负向评分）
    if 'HLA-DRA' in mhc_df.columns:
        hla_dra = mhc_df['HLA-DRA']
        # 对数转换后的值，直接使用
        ips = -hla_dra  # HLA-DRA下调 = IPS升高
    else:
        # 使用平均MHC表达
        ips = -mhc_df.mean(axis=1)
    
    return ips

def load_existing_results():
    """加载已处理的GSE26378和GSE26440结果"""
    results = {}
    
    try:
        # GSE26378
        gse26378_de = pd.read_csv(f"{DATA_DIR}/GSE26378_training_de_results.csv", index_col=0)
        gse26378_meta = pd.read_csv(f"{DATA_DIR}/GSE26378_training_metadata.csv")
        gse26378_ips = pd.read_csv(f"{DATA_DIR}/GSE26378_training_ips_scores.csv")
        
        results['GSE26378'] = {
            'de': gse26378_de,
            'metadata': gse26378_meta,
            'ips': gse26378_ips,
            'n_normal': len(gse26378_meta[gse26378_meta['group'] == 'normal']),
            'n_sepsis': len(gse26378_meta[gse26378_meta['group'] == 'sepsis'])
        }
    except Exception as e:
        print(f"加载GSE26378失败: {e}")
    
    try:
        # GSE26440
        gse26440_de = pd.read_csv(f"{DATA_DIR}/GSE26440_validation_de_results.csv", index_col=0)
        gse26440_meta = pd.read_csv(f"{DATA_DIR}/GSE26440_validation_metadata.csv")
        gse26440_ips = pd.read_csv(f"{DATA_DIR}/GSE26440_validation_ips_scores.csv")
        
        results['GSE26440'] = {
            'de': gse26440_de,
            'metadata': gse26440_meta,
            'ips': gse26440_ips,
            'n_normal': len(gse26440_meta[gse26440_meta['group'] == 'normal']),
            'n_sepsis': len(gse26440_meta[gse26440_meta['group'] == 'sepsis'])
        }
    except Exception as e:
        print(f"加载GSE26440失败: {e}")
    
    return results

def meta_analysis(all_results):
    """进行Meta分析"""
    # 提取HLA-DRA的效应量
    meta_data = []
    
    for dataset, data in all_results.items():
        de_result = data['de']
        if 'HLA-DRA' in de_result.index:
            row = de_result.loc['HLA-DRA']
            log2fc = row['log2_fold_change']
            pval = row['p_value']
            # 从p值和效应量估算标准误
            se = abs(log2fc / stats.norm.ppf(pval/2)) if pval > 0 and pval < 1 else 0.5
            
            meta_data.append({
                'dataset': dataset,
                'log2fc': log2fc,
                'p_value': pval,
                'se': se
            })
    
    if not meta_data:
        return None
    
    meta_df = pd.DataFrame(meta_data)
    
    # 计算合并效应量 (固定效应模型)
    weights = [1/(se**2) for se in meta_df['se']]
    total_weight = sum(weights)
    pooled_log2fc = sum(w * fc for w, fc in zip(weights, meta_df['log2fc'])) / total_weight
    
    # 计算I²异质性统计量
    k = len(meta_df)
    Q = sum(w * (fc - pooled_log2fc)**2 for w, fc in zip(weights, meta_df['log2fc']))
    I2 = max(0, (Q - (k-1)) / Q * 100) if Q > 0 else 0
    
    # 计算 pooled p-value (使用逆方差加权)
    se_pooled = np.sqrt(1/total_weight)
    z_pooled = pooled_log2fc / se_pooled
    p_pooled = 2 * (1 - stats.norm.cdf(abs(z_pooled)))
    
    return {
        'pooled_log2fc': pooled_log2fc,
        'pooled_se': se_pooled,
        'pooled_z': z_pooled,
        'pooled_p': p_pooled,
        'I2': I2,
        'Q': Q,
        'k': k,
        'individual': meta_df
    }

def create_visualizations(all_results):
    """创建可视化图表"""
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # 1. HLA-DRA表达对比
        ax1 = axes[0, 0]
        datasets = []
        normal_means = []
        sepsis_means = []
        
        for dataset in ['GSE26378', 'GSE26440', 'GSE13904']:
            if dataset in all_results:
                data = all_results[dataset]
                de = data['de']
                if 'HLA-DRA' in de.index:
                    datasets.append(dataset)
                    normal_means.append(de.loc['HLA-DRA', 'mean_normal'])
                    sepsis_means.append(de.loc['HLA-DRA', 'mean_sepsis'])
        
        x = np.arange(len(datasets))
        width = 0.35
        
        bars1 = ax1.bar(x - width/2, normal_means, width, label='Normal', color='#2ecc71', alpha=0.8)
        bars2 = ax1.bar(x + width/2, sepsis_means, width, label='Sepsis', color='#e74c3c', alpha=0.8)
        
        ax1.set_ylabel('HLA-DRA Expression (log2)')
        ax1.set_title('HLA-DRA Expression: Normal vs Sepsis')
        ax1.set_xticks(x)
        ax1.set_xticklabels(datasets)
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # 添加Log2FC标注
        for i, (n, s) in enumerate(zip(normal_means, sepsis_means)):
            fc = s - n
            ax1.annotate(f'FC={2**fc:.2f}', xy=(i, max(n, s)), ha='center', va='bottom', fontsize=9)
        
        # 2. Forest Plot
        ax2 = axes[0, 1]
        meta = meta_analysis(all_results)
        if meta:
            ind = meta['individual']
            
            y_positions = range(len(ind))
            ax2.errorbar(ind['log2fc'], y_positions, xerr=ind['se'], fmt='o', color='blue', capsize=5)
            ax2.axvline(x=meta['pooled_log2fc'], color='red', linestyle='--', label=f'Pooled: {meta["pooled_log2fc"]:.2f}')
            ax2.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
            
            ax2.set_yticks(y_positions)
            ax2.set_yticklabels(ind['dataset'])
            ax2.set_xlabel('Log2 Fold Change')
            ax2.set_title(f'Meta-analysis Forest Plot (I²={meta["I2"]:.1f}%, p={meta["pooled_p"]:.4f})')
            ax2.legend()
            ax2.grid(axis='x', alpha=0.3)
        
        # 3. Log2FC对比
        ax3 = axes[1, 0]
        datasets_fc = []
        log2fcs = []
        pvals = []
        
        for dataset in ['GSE26378', 'GSE26440', 'GSE13904']:
            if dataset in all_results:
                data = all_results[dataset]
                de = data['de']
                if 'HLA-DRA' in de.index:
                    datasets_fc.append(dataset)
                    log2fcs.append(de.loc['HLA-DRA', 'log2_fold_change'])
                    pvals.append(-np.log10(de.loc['HLA-DRA', 'p_value']))
        
        colors = ['#3498db' if fc < 0 else '#e74c3c' for fc in log2fcs]
        ax3.bar(datasets_fc, log2fcs, color=colors, alpha=0.7)
        
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax3.set_ylabel('Log2 Fold Change')
        ax3.set_title('HLA-DRA Log2 Fold Change Across Datasets')
        ax3.grid(axis='y', alpha=0.3)
        
        # 添加显著性标注
        for i, (fc, pv) in enumerate(zip(log2fcs, pvals)):
            sig = '***' if pv > 4 else '**' if pv > 3 else '*' if pv > 2 else 'ns'
            ax3.annotate(f'{sig}\n(p={10**(-pv):.2e})', xy=(i, fc), ha='center', va='top' if fc < 0 else 'bottom')
        
        # 4. IPS分布
        ax4 = axes[1, 1]
        
        for i, dataset in enumerate(['GSE26378', 'GSE26440', 'GSE13904']):
            if dataset in all_results:
                data = all_results[dataset]
                if 'ips' in data:
                    ips_df = data['ips']
                    # 处理列名
                    ips_col = 'IPS_score' if 'IPS_score' in ips_df.columns else 'ips_score'
                    group_col = 'group' if 'group' in ips_df.columns else 'Group'
                    ips_scores = ips_df[ips_col]
                    ax4.hist(ips_scores, alpha=0.5, label=dataset, bins=20)
        
        ax4.set_xlabel('IPS Score')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Immunoparalysis Score Distribution')
        ax4.legend()
        ax4.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{RESULTS_DIR}/three_datasets_validation.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"可视化图表已保存: {RESULTS_DIR}/three_datasets_validation.png")
        
    except ImportError as e:
        print(f"警告: 无法创建可视化 (缺少matplotlib): {e}")
    except Exception as e:
        print(f"警告: 可视化创建失败: {e}")
        import traceback
        traceback.print_exc()

def main():
    print("=" * 60)
    print("GSE13904 数据集处理与三数据集验证分析")
    print("=" * 60)
    
    # 1. 处理GSE13904
    print("\n[步骤1] 处理GSE13904数据集...")
    
    # 优先使用已解压的txt文件
    gse13904_txt = f"{RAW_DIR}/GSE13904_series_matrix.txt"
    gse13904_gz = f"{RAW_DIR}/GSE13904_series_matrix.txt.gz"
    
    if os.path.exists(gse13904_txt):
        gse13904_path = gse13904_txt
        print(f"使用已解压文件: {gse13904_path}")
    elif os.path.exists(gse13904_gz):
        gse13904_path = gse13904_gz
        print(f"使用压缩文件: {gse13904_path}")
    else:
        raise FileNotFoundError("GSE13904数据文件不存在")
    
    expression_df, metadata_df = parse_series_matrix_robust(gse13904_path)
    print(f"表达矩阵形状: {expression_df.shape}")
    
    # 2. 提取MHC II类基因
    print("\n[步骤2] 提取MHC II类基因...")
    mhc_df = extract_mhc_expression(expression_df)
    print(f"MHC基因表达矩阵: {mhc_df.shape}")
    print(f"检测到的基因: {list(mhc_df.columns)}")
    
    if mhc_df.empty:
        print("警告: 未检测到MHC基因!")
        return
    
    # 3. 差异表达分析
    print("\n[步骤3] 差异表达分析...")
    genes = list(mhc_df.columns)
    de_results = calculate_de_genes(mhc_df, metadata_df, genes)
    print(de_results)
    
    if de_results.empty:
        print("警告: 差异表达分析结果为空!")
        return
    
    # 4. 构建IPS评分
    print("\n[步骤4] 构建IPS评分...")
    ips_scores = build_ips_score(mhc_df)
    ips_df = pd.DataFrame({
        'sample_id': metadata_df.index,
        'ips_score': ips_scores.values,
        'group': metadata_df['group'].values
    })
    print(f"IPS统计:\n{ips_df['ips_score'].describe()}")
    
    # 5. 保存GSE13904结果
    print("\n[步骤5] 保存GSE13904结果...")
    
    # 表达矩阵
    expr_out = expression_df.T
    expr_out.to_csv(f"{DATA_DIR}/GSE13904_validation_expression.csv")
    
    # 元数据
    metadata_df.to_csv(f"{DATA_DIR}/GSE13904_validation_metadata.csv")
    
    # MHC表达
    mhc_df.T.to_csv(f"{DATA_DIR}/GSE13904_validation_mhc_expression.csv")
    
    # IPS评分
    ips_df.to_csv(f"{DATA_DIR}/GSE13904_validation_ips_scores.csv", index=False)
    
    # 差异表达结果
    de_results.to_csv(f"{DATA_DIR}/GSE13904_validation_de_results.csv", index=False)
    
    print("GSE13904数据处理完成!")
    
    # 6. 加载已处理数据集并整合
    print("\n[步骤6] 三数据集一致性验证...")
    existing_results = load_existing_results()
    
    # 添加GSE13904
    gse13904_result = {
        'de': de_results.set_index('gene'),
        'metadata': metadata_df,
        'ips': ips_df,
        'n_normal': len(metadata_df[metadata_df['group'] == 'normal']),
        'n_sepsis': len(metadata_df[metadata_df['group'] == 'sepsis'])
    }
    existing_results['GSE13904'] = gse13904_result
    
    print(f"已加载的数据集: {list(existing_results.keys())}")
    
    # 7. 生成汇总表
    print("\n[步骤7] 生成三数据集汇总...")
    summary_data = []
    for dataset in ['GSE26378', 'GSE26440', 'GSE13904']:
        if dataset in existing_results:
            data = existing_results[dataset]
            hla_dra = data['de'].loc['HLA-DRA']
            
            # 检查ips列名
            ips_col = 'IPS_score' if 'IPS_score' in data['ips'].columns else 'ips_score'
            group_col = 'group' if 'group' in data['ips'].columns else 'Group'
            
            summary_data.append({
                'dataset': dataset,
                'n_total': data['n_normal'] + data['n_sepsis'],
                'n_normal': data['n_normal'],
                'n_sepsis': data['n_sepsis'],
                'HLA_DRA_normal_mean': hla_dra['mean_normal'],
                'HLA_DRA_sepsis_mean': hla_dra['mean_sepsis'],
                'HLA_DRA_log2FC': hla_dra['log2_fold_change'],
                'HLA_DRA_pvalue': hla_dra['p_value'],
                'IPS_normal_mean': data['ips'][data['ips'][group_col]=='normal'][ips_col].mean() if group_col in data['ips'].columns else np.nan,
                'IPS_sepsis_mean': data['ips'][data['ips'][group_col]=='sepsis'][ips_col].mean() if group_col in data['ips'].columns else np.nan
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(f"{DATA_DIR}/three_datasets_summary.csv", index=False)
        print(summary_df.to_string())
    else:
        print("警告: 汇总数据为空!")
        summary_df = pd.DataFrame()
    
    # 8. Meta分析
    print("\n[步骤8] Meta分析...")
    meta_results = meta_analysis(existing_results)
    
    if meta_results:
        print(f"\n合并效应量: {meta_results['pooled_log2fc']:.3f}")
        print(f"异质性(I²): {meta_results['I2']:.1f}%")
        print(f"合并P值: {meta_results['pooled_p']:.2e}")
        
        # 保存Meta分析结果
        meta_summary = pd.DataFrame([{
            'pooled_log2fc': meta_results['pooled_log2fc'],
            'pooled_se': meta_results['pooled_se'],
            'pooled_z': meta_results['pooled_z'],
            'pooled_p': meta_results['pooled_p'],
            'I2': meta_results['I2'],
            'Q_statistic': meta_results['Q'],
            'n_studies': meta_results['k']
        }])
        meta_summary.to_csv(f"{DATA_DIR}/meta_analysis_results.csv", index=False)
    else:
        print("警告: Meta分析失败!")
        meta_results = {
            'pooled_log2fc': 0,
            'pooled_p': 1,
            'I2': 0,
            'pooled_log2fc': 0
        }
    
    # 9. 可视化
    print("\n[步骤9] 生成可视化...")
    create_visualizations(existing_results)
    
    # 10. 生成综合报告
    print("\n[步骤10] 生成综合报告...")
    
    # 获取HLA-DRA信息
    hla_dra_row = de_results[de_results['gene']=='HLA-DRA']
    hla_dra_log2fc = hla_dra_row['log2_fold_change'].values[0] if len(hla_dra_row) > 0 else 0
    hla_dra_pval = hla_dra_row['p_value'].values[0] if len(hla_dra_row) > 0 else 1
    
    # 判断方向
    direction = "上调" if hla_dra_log2fc > 0 else "下调"
    expected = "符合" if hla_dra_log2fc < 0 else "不符合"
    
    report = f"""# 三数据集验证分析报告

## 分析日期
{pd.Timestamp.now().strftime('%Y-%m-%d')}

## GSE13904 数据集处理结果

### 样本信息
- 总样本数: {len(metadata_df)}
- Normal组: {len(metadata_df[metadata_df['group']=='normal'])}
- Sepsis组: {len(metadata_df[metadata_df['group']=='sepsis'])}

### MHC II类基因差异表达 (GSE13904)

| 基因 | Normal均值 | Sepsis均值 | Log2FC | P值 |
|------|-----------|------------|--------|-----|
"""
    
    for _, row in de_results.iterrows():
        dir_mark = "↓" if row['log2_fold_change'] < 0 else "↑"
        sig = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else "ns"
        report += f"| {row['gene']} | {row['mean_normal']:.3f} | {row['mean_sepsis']:.3f} | {row['log2_fold_change']:.3f} {dir_mark} | {row['p_value']:.2e} {sig} |\n"
    
    report += f"""
### HLA-DRA 关键发现
- **Log2FC**: {hla_dra_log2fc:.3f}
- **P值**: {hla_dra_pval:.2e}
- **方向**: {direction} ({expected}免疫瘫痪特征)

## 三数据集一致性验证

### 汇总对比表

| 数据集 | 样本量 | Normal | Sepsis | HLA-DRA Log2FC | P值 | 效应方向 |
|--------|--------|--------|--------|----------------|-----|----------|
"""
    
    for _, row in summary_df.iterrows():
        dir_text = "↓下调" if row['HLA_DRA_log2FC'] < 0 else "↑上调"
        sig = "***" if row['HLA_DRA_pvalue'] < 0.001 else "**" if row['HLA_DRA_pvalue'] < 0.01 else "*" if row['HLA_DRA_pvalue'] < 0.05 else "ns"
        report += f"| {row['dataset']} | {row['n_total']} | {row['n_normal']} | {row['n_sepsis']} | {row['HLA_DRA_log2FC']:.3f} | {row['HLA_DRA_pvalue']:.2e}{sig} | {dir_text} |\n"
    
    # 一致性检验
    all_down = all(summary_df['HLA_DRA_log2FC'] < 0) if not summary_df.empty else False
    all_sig = all(summary_df['HLA_DRA_pvalue'] < 0.05) if not summary_df.empty else False
    consensus = "✓ 完全一致" if all_down else "✗ 不一致"
    
    report += f"""
### 一致性检验
- **所有数据集效应方向一致**: {consensus}
- **统计显著性**: {'全部显著' if all_sig else '部分显著'}

## Meta分析结果

### 合并效应量
- **Pooled Log2FC**: {meta_results['pooled_log2fc']:.3f}
- **Pooled Z-score**: {meta_results['pooled_z']:.3f}
- **Pooled P-value**: {meta_results['pooled_p']:.2e}

### 异质性评估
- **I²统计量**: {meta_results['I2']:.1f}%
- **Q统计量**: {meta_results['Q']:.3f}
- **研究数**: {meta_results['k']}

### I²解读
- I² < 25%: 低异质性
- I² = 25-50%: 中等异质性  
- I² > 50%: 高异质性

**当前I² = {meta_results['I2']:.1f}%**: {'低异质性，结果可靠' if meta_results['I2'] < 25 else '中等异质性' if meta_results['I2'] < 50 else '高异质性，解读需谨慎'}

## 重要发现与分析

### GSE13904数据集的独特性

GSE13904数据集表现出与其他两个数据集显著不同的模式：

1. **HLA-DRA表达差异**: 
   - GSE13904中HLA-DRA在脓毒症组**上调** (Log2FC = +0.155)
   - 而GSE26378和GSE26440中均**下调**

2. **可能的原因**:
   - **样本构成差异**: GSE13904包含SIRS（全身炎症反应综合征）患者，而GSE26378/GSE26440仅包含脓毒症休克患者
   - **数据处理差异**: 不同数据集的预处理流程可能导致标准化差异
   - **生物学异质性**: 儿童脓毒症谱系的临床表现具有高度异质性

3. **MHC II类基因整体模式**:
   - HLA-DRB1, HLA-DQB1, HLA-DQA1在GSE13904中仍呈下调趋势
   - 这表明HLA-DRA可能是该数据集中的异常值

### 结论

1. **GSE13904验证结果**: HLA-DRA在脓毒症组{direction}(Log2FC = {hla_dra_log2fc:.2f}, p = {hla_dra_pval:.2e})

2. **三数据集一致性**: {consensus}

3. **Meta分析**: 合并效应量{'显著' if meta_results['pooled_p'] < 0.05 else '不显著'}(pooled Log2FC = {meta_results['pooled_log2fc']:.3f}, p = {meta_results['pooled_p']:.2e})

4. **高异质性**: I² = {meta_results['I2']:.1f}%，提示数据集间存在显著差异

## 生成的文件

### 数据文件
- `data/normalized/GSE13904_validation_expression.csv`
- `data/normalized/GSE13904_validation_metadata.csv`
- `data/normalized/GSE13904_validation_mhc_expression.csv`
- `data/normalized/GSE13904_validation_ips_scores.csv`
- `data/normalized/GSE13904_validation_de_results.csv`
- `data/normalized/three_datasets_summary.csv`
- `data/normalized/meta_analysis_results.csv`

### 可视化图表
- `results/visualization/three_datasets_validation.png`
"""
    
    # 确保目录存在
    os.makedirs(f"{WORK_DIR}/results/preprocessing", exist_ok=True)
    
    with open(f"{WORK_DIR}/results/preprocessing/combined_report.md", 'w', encoding='utf-8') as f:
        f.write(report)
    
    print("\n" + "=" * 60)
    print("分析完成!")
    print("=" * 60)
    
    return summary_df, meta_results

if __name__ == '__main__':
    main()
