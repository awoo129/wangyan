#!/usr/bin/env python3
"""
路线C：外部验证数据集分析
分析GSE66099和GSE13904中的CIITA和MHC II基因表达
"""

import re
import gzip
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# MHC II相关基因
MHC_II_GENES = {
    'CIITA': 'class II transactivator',
    'HLA-DRA': 'major histocompatibility complex class II DR alpha',
    'HLA-DPA1': 'HLA class II DP alpha 1',
    'HLA-DPB1': 'HLA class II DP beta 1',
    'HLA-DQA1': 'HLA class II DQ alpha 1',
    'HLA-DQB1': 'HLA class II DQ beta 1',
    'HLA-DMA': 'HLA class II DM alpha',
    'HLA-DMB': 'HLA class II DM beta',
    'HLA-DOB': 'HLA class II DO beta',
    'CD74': 'invariant chain (Ii)'
}

# GPL570探针注释（基于已知信息）
# CIITA: 208981_at, 208983_s_at, 214523_at
# HLA-DRA: 212828_s_at, 214486_at
PROBE_MAP = {
    'CIITA': ['208981_at', '208983_s_at', '214523_at'],
    'HLA-DRA': ['212828_s_at', '214486_at', '201137_s_at'],
    'HLA-DPA1': ['210044_s_at', '209823_at'],
    'HLA-DPB1': ['204790_s_at', '210331_x_at'],
    'HLA-DQA1': ['209480_x_at', '211990_x_at'],
    'HLA-DQB1': ['209480_x_at', '210331_x_at', '209482_at'],
    'HLA-DMA': ['202411_at', '208691_s_at'],
    'HLA-DMB': ['202409_s_at', '203798_at'],
    'HLA-DOB': ['206239_at', '214038_x_at'],
    'CD74': ['200664_at', '212068_at', '201009_s_at']
}

def parse_series_matrix(filepath, is_gz=False):
    """解析series_matrix文件"""
    print(f"\n解析文件: {filepath}")
    
    data_lines = []
    sample_info = {}
    
    if is_gz:
        opener = gzip.open
    else:
        opener = open
    
    with opener(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!') and not line.startswith('!Sample_table_begin'):
                # 解析样本元信息
                match = re.match(r'(!Sample_\w+)\t(.+)', line)
                if match:
                    key = match.group(1).replace('!', '')
                    values = [v.strip('"') for v in match.group(2).split('\t')]
                    sample_info[key] = values
            elif line.startswith('!Sample_table_begin'):
                break
    
    return sample_info

def extract_expression(filepath, is_gz=False):
    """提取表达数据"""
    print(f"\n提取表达数据: {filepath}")
    
    probe_ids = []
    expression_data = []
    
    if is_gz:
        opener = gzip.open
    else:
        opener = open
    
    reading_data = False
    with opener(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!Sample_table_begin'):
                reading_data = True
                continue
            elif line.startswith('!') and reading_data:
                break
            elif reading_data:
                parts = line.split('\t')
                if parts:
                    probe_id = parts[0].strip('"')
                    try:
                        expr_values = [float(v.strip('"')) for v in parts[1:]]
                        probe_ids.append(probe_id)
                        expression_data.append(expr_values)
                    except:
                        continue
    
    print(f"  提取到 {len(probe_ids)} 个探针, {len(expression_data[0]) if expression_data else 0} 个样本")
    return probe_ids, expression_data

def extract_gse13904_samples(sample_titles):
    """从GSE13904样本标题中提取分组信息"""
    samples = []
    for title in sample_titles:
        parts = title.split('_')
        if len(parts) >= 3:
            sample_id = parts[0]
            group = parts[1]
            timepoint = parts[2] if len(parts) > 2 else 'unknown'
            
            # 标准化分组名称
            if 'Control' in group or 'control' in group:
                group = 'Control'
            elif 'SIRS' in group or 'sirs' in group:
                if 'resolved' in title:
                    group = 'SIRS_resolved'
                else:
                    group = 'SIRS'
            elif 'Sepsis' in group or 'sepsis' in group:
                group = 'Sepsis'
            elif 'Septic' in group or 'septic' in group:
                group = 'SepticShock'
            
            samples.append({
                'title': title,
                'group': group,
                'timepoint': timepoint
            })
    return samples

def analyze_gene_expression(probe_ids, expression_data, gene_probes):
    """分析基因表达"""
    results = {}
    for gene, probes in gene_probes.items():
        gene_expr = []
        for probe in probes:
            if probe in probe_ids:
                idx = probe_ids.index(probe)
                gene_expr.append(np.array(expression_data[idx]))
        if gene_expr:
            results[gene] = np.mean(gene_expr, axis=0)
    return results

def compare_groups(expression, sample_groups, group1, group2):
    """比较两组表达差异"""
    idx1 = [i for i, g in enumerate(sample_groups) if g == group1]
    idx2 = [i for i, g in enumerate(sample_groups) if g == group2]
    
    if len(idx1) == 0 or len(idx2) == 0:
        return None, None, None
    
    expr1 = expression[idx1]
    expr2 = expression[idx2]
    
    # t检验
    t_stat, p_value = stats.ttest_ind(expr1, expr2)
    
    # 效应量 (Cohen's d)
    pooled_std = np.sqrt((np.std(expr1)**2 + np.std(expr2)**2) / 2)
    cohens_d = (np.mean(expr1) - np.mean(expr2)) / pooled_std if pooled_std > 0 else 0
    
    return np.mean(expr1), np.mean(expr2), p_value, cohens_d

def main():
    print("="*80)
    print("路线C: 外部验证数据集分析")
    print("="*80)
    
    # ============ GSE13904分析 ============
    print("\n" + "="*80)
    print("GSE13904 数据集分析")
    print("="*80)
    
    gse13904_path = "长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt"
    sample_info = parse_series_matrix(gse13904_path, is_gz=False)
    probe_ids, expression_data = extract_expression(gse13904_path, is_gz=False)
    
    # 提取样本分组
    sample_titles = sample_info.get('title', [])
    samples = extract_gse13904_samples(sample_titles)
    sample_groups = [s['group'] for s in samples]
    
    # 统计各组样本数
    print("\n样本分组统计:")
    group_counts = {}
    for g in sample_groups:
        group_counts[g] = group_counts.get(g, 0) + 1
    for g, count in sorted(group_counts.items()):
        print(f"  {g}: {count}")
    
    # 分析MHC II基因表达
    print("\nMHC II基因表达分析 (GSE13904 - 儿童数据):")
    gene_expr = analyze_gene_expression(probe_ids, expression_data, PROBE_MAP)
    
    # 创建表达矩阵
    expr_df = pd.DataFrame(gene_expr).T
    expr_df.columns = sample_titles
    expr_df['group'] = sample_groups
    
    results_13904 = {}
    for gene in gene_expr:
        results_13904[gene] = {}
        print(f"\n{gene} ({MHC_II_GENES.get(gene, '')}):")
        for group in ['Control', 'SIRS', 'Sepsis', 'SepticShock']:
            group_expr = [gene_expr[gene][i] for i in range(len(sample_groups)) if sample_groups[i] == group]
            if group_expr:
                results_13904[gene][group] = {
                    'mean': np.mean(group_expr),
                    'std': np.std(group_expr),
                    'n': len(group_expr)
                }
                print(f"  {group}: {np.mean(group_expr):.3f} ± {np.std(group_expr):.3f} (n={len(group_expr)})")
        
        # Sepsis vs Control 比较
        if 'Control' in results_13904[gene] and 'Sepsis' in results_13904[gene]:
            _, _, p_val, d = compare_groups(gene_expr[gene], sample_groups, 'Sepsis', 'Control')
            print(f"  Sepsis vs Control: p={p_val:.4f}, Cohen's d={d:.3f}")
    
    # ============ GSE66099分析 ============
    print("\n" + "="*80)
    print("GSE66099 数据集分析")
    print("="*80)
    
    gse66099_path = "长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE66099_series_matrix.txt.gz"
    sample_info_66099 = parse_series_matrix(gse66099_path, is_gz=True)
    probe_ids_66099, expression_data_66099 = extract_expression(gse66099_path, is_gz=True)
    
    # GSE66099样本标题格式不同，需要解析
    sample_titles_66099 = sample_info_66099.get('title', [])
    
    # 从样本标题中推断分组（GSE66099包含多个原始研究的合并数据）
    # 根据PMID 25972003和27635771，该数据集包含SIRS、sepsis、septic shock和健康对照
    print("\n分析GSE66099样本...")
    print(f"总样本数: {len(sample_titles_66099)}")
    
    # GSE66099的样本标题格式: "XX-XXXX" 格式，没有明确的分组信息
    # 需要从其他注释字段中提取或搜索外部信息
    print("\n样本标题示例:")
    for i, title in enumerate(sample_titles_66099[:10]):
        print(f"  {title}")
    
    # 分析MHC II基因表达
    print("\nMHC II基因表达分析 (GSE66099):")
    gene_expr_66099 = analyze_gene_expression(probe_ids_66099, expression_data_66099, PROBE_MAP)
    
    results_66099 = {}
    for gene in gene_expr_66099:
        results_66099[gene] = {
            'mean': np.mean(gene_expr_66099[gene]),
            'std': np.std(gene_expr_66099[gene]),
            'median': np.median(gene_expr_66099[gene]),
            'range': (np.min(gene_expr_66099[gene]), np.max(gene_expr_66099[gene]))
        }
        print(f"  {gene}: mean={results_66099[gene]['mean']:.3f}, std={results_66099[gene]['std']:.3f}")
    
    # ============ 比较分析 ============
    print("\n" + "="*80)
    print("年龄异质性验证")
    print("="*80)
    
    # 检查是否有年龄信息
    age_info = sample_info_66099.get('characteristics_ch1', [])
    print(f"\nGSE66099样本特征信息存在: {'是' if age_info else '否'}")
    if age_info and len(age_info) > 0:
        print("特征信息示例:")
        for i, char in enumerate(age_info[:5]):
            print(f"  {char}")
    
    # 分析结论
    print("\n" + "="*80)
    print("分析结论")
    print("="*80)
    
    print("""
    1. GSE13904 数据集评估:
       - 样本组成: 儿童对照、SIRS、脓毒症、脓毒性休克
       - 可用于分析CIITA等MHC II基因在不同疾病严重程度下的表达
       - 适合与路线A(儿童IPS)结果对比
    
    2. GSE66099 数据集评估:
       - 这是6个儿童SIRS/脓毒症研究的合并数据集
       - 包含276个样本(Day 1)
       - 问题: 原始series_matrix文件中缺少明确的分组标注
       - 需要额外元数据或原始CEL文件才能确定样本分组
    
    3. 年龄异质性验证状态:
       - 路线A: 儿童脓毒症数据已分析
       - 路线B: 成人脓毒症数据已分析
       - 当前两个外部验证数据集都是儿童数据
       - 需要成人外部数据集进行年龄分层验证
    """)

if __name__ == "__main__":
    main()
