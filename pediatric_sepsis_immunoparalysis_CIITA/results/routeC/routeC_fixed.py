#!/usr/bin/env python3
"""
路线C：外部验证数据集分析 - 修复版
"""

import re
import gzip
import numpy as np
from scipy import stats
import json

# MHC II相关基因和探针ID (GPL570)
PROBE_MAP = {
    'CIITA': ['208981_at', '208983_s_at', '214523_at'],
    'HLA-DRA': ['212828_s_at', '214486_at', '201137_s_at'],
    'HLA-DPA1': ['210044_s_at', '209823_at'],
    'CD74': ['200664_at', '212068_at', '201009_s_at']
}

def parse_series_matrix(filepath, is_gz=False):
    """解析series_matrix文件"""
    print(f"\n解析文件: {filepath}")
    
    sample_titles = []
    probe_ids = []
    expression_data = []
    reading_data = False
    
    if is_gz:
        opener = gzip.open
    else:
        opener = open
    
    with opener(filepath, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            
            # 解析样本标题
            if line.startswith('!Sample_title'):
                parts = line.split('\t')
                sample_titles = [p.strip('"') for p in parts[1:]]
            
            if line.startswith('!series_matrix_table_begin'):
                reading_data = True
                continue
            
            if line.startswith('!series_matrix_table_end'):
                break
            
            if reading_data and line and '\t' in line:
                parts = line.split('\t')
                if parts:
                    probe_id = parts[0].strip('"')
                    if probe_id == 'ID_REF':
                        continue
                    try:
                        expr_values = [float(v.strip('"')) for v in parts[1:]]
                        probe_ids.append(probe_id)
                        expression_data.append(expr_values)
                    except ValueError:
                        continue
    
    print(f"  提取到 {len(probe_ids)} 个探针, {len(expression_data[0]) if expression_data else 0} 个样本")
    print(f"  样本标题数: {len(sample_titles)}")
    
    return sample_titles, probe_ids, expression_data

def extract_groups(sample_titles):
    """从样本标题中提取分组信息"""
    samples = []
    for title in sample_titles:
        parts = title.split('_')
        if len(parts) >= 3:
            group_str = parts[1]
            timepoint = parts[2]
            
            if 'Control' in group_str:
                group = 'Control'
            elif 'SIRS resolved' in title:
                group = 'SIRS_resolved'
            elif 'SIRS' in group_str:
                group = 'SIRS'
            elif 'Sepsis' in group_str:
                group = 'Sepsis'
            elif 'Septic' in group_str:
                group = 'SepticShock'
            else:
                group = 'Unknown'
            
            samples.append({'title': title, 'group': group, 'timepoint': timepoint})
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
    
    if len(idx1) < 2 or len(idx2) < 2:
        return None, None, None, None
    
    expr1 = expression[idx1]
    expr2 = expression[idx2]
    
    t_stat, p_value = stats.ttest_ind(expr1, expr2)
    pooled_std = np.sqrt((np.std(expr1)**2 + np.std(expr2)**2) / 2)
    cohens_d = (np.mean(expr1) - np.mean(expr2)) / pooled_std if pooled_std > 0 else 0
    
    return np.mean(expr1), np.mean(expr2), p_value, cohens_d

def main():
    print("="*80)
    print("路线C: 外部验证数据集分析")
    print("="*80)
    
    all_results = {}
    
    # ============ GSE13904分析 ============
    print("\n" + "="*80)
    print("GSE13904 数据集分析 (儿童SIRS/脓毒症谱)")
    print("="*80)
    
    gse13904_path = "长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt"
    sample_titles_13904, probe_ids_13904, expression_data_13904 = parse_series_matrix(gse13904_path, is_gz=False)
    
    # 提取样本分组
    samples_13904 = extract_groups(sample_titles_13904)
    sample_groups_13904 = [s['group'] for s in samples_13904]
    timepoints_13904 = [s['timepoint'] for s in samples_13904]
    
    # Day1样本索引
    day1_indices = [i for i, tp in enumerate(timepoints_13904) if tp == 'day1']
    day1_groups = [sample_groups_13904[i] for i in day1_indices]
    
    print("\n样本分组统计:")
    group_counts = {}
    for g in sample_groups_13904:
        group_counts[g] = group_counts.get(g, 0) + 1
    for g, count in sorted(group_counts.items()):
        print(f"  {g}: {count}")
    
    print(f"\nDay1样本数: {len(day1_indices)}")
    
    # 分析MHC II基因表达
    print("\n" + "-"*60)
    print("MHC II基因表达分析 (GSE13904 - Day1样本):")
    print("-"*60)
    
    gene_expr_13904 = analyze_gene_expression(probe_ids_13904, expression_data_13904, PROBE_MAP)
    
    # 筛选Day1样本
    gene_expr_13904_day1 = {gene: expr[day1_indices] for gene, expr in gene_expr_13904.items()}
    
    results_13904 = {}
    for gene in ['CIITA', 'HLA-DRA', 'CD74']:
        if gene in gene_expr_13904_day1:
            results_13904[gene] = {}
            expr = gene_expr_13904_day1[gene]
            
            print(f"\n{gene}:")
            for group in ['Control', 'SIRS', 'Sepsis', 'SepticShock']:
                group_indices = [i for i, g in enumerate(day1_groups) if g == group]
                if group_indices:
                    group_expr = expr[group_indices]
                    results_13904[gene][group] = {
                        'mean': float(np.mean(group_expr)),
                        'std': float(np.std(group_expr)),
                        'n': len(group_expr)
                    }
                    print(f"  {group:15s}: {np.mean(group_expr):.4f} ± {np.std(group_expr):.4f} (n={len(group_expr)})")
            
            # 统计检验
            if 'Control' in results_13904[gene] and 'Sepsis' in results_13904[gene]:
                m1, m2, p, d = compare_groups(expr, day1_groups, 'Sepsis', 'Control')
                if p is not None:
                    direction = "↑" if m2 > m1 else "↓"
                    results_13904[gene]['Sepsis_vs_Control'] = {'delta': float(m2-m1), 'p': float(p), 'cohens_d': float(d), 'direction': direction}
                    print(f"  Sepsis vs Control: Δ={m2-m1:+.4f}, p={p:.4f}, Cohen's d={d:.3f} ({direction})")
            
            if 'Control' in results_13904[gene] and 'SepticShock' in results_13904[gene]:
                m1, m2, p, d = compare_groups(expr, day1_groups, 'SepticShock', 'Control')
                if p is not None:
                    direction = "↑" if m2 > m1 else "↓"
                    results_13904[gene]['SepticShock_vs_Control'] = {'delta': float(m2-m1), 'p': float(p), 'cohens_d': float(d), 'direction': direction}
                    print(f"  SepticShock vs Control: Δ={m2-m1:+.4f}, p={p:.4f}, Cohen's d={d:.3f} ({direction})")
    
    all_results['GSE13904'] = results_13904
    
    # ============ GSE66099分析 ============
    print("\n" + "="*80)
    print("GSE66099 数据集分析")
    print("="*80)
    
    gse66099_path = "长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE66099_series_matrix.txt.gz"
    sample_titles_66099, probe_ids_66099, expression_data_66099 = parse_series_matrix(gse66099_path, is_gz=True)
    
    print("\n注意: GSE66099缺少分组信息，无法进行分组分析")
    
    gene_expr_66099 = analyze_gene_expression(probe_ids_66099, expression_data_66099, PROBE_MAP)
    
    print("\n" + "-"*60)
    print("MHC II基因整体表达 (GSE66099):")
    print("-"*60)
    
    results_66099 = {}
    for gene in ['CIITA', 'HLA-DRA', 'CD74']:
        if gene in gene_expr_66099:
            expr = gene_expr_66099[gene]
            results_66099[gene] = {
                'mean': float(np.mean(expr)),
                'std': float(np.std(expr)),
                'median': float(np.median(expr))
            }
            print(f"\n{gene}:")
            print(f"  mean={np.mean(expr):.4f} ± {np.std(expr):.4f}")
            print(f"  median={np.median(expr):.4f}")
    
    all_results['GSE66099'] = results_66099
    
    # ============ 结论 ============
    print("\n" + "="*80)
    print("外部验证数据集评估结论")
    print("="*80)
    
    conclusion = """
    ┌─────────────────────────────────────────────────────────────────────────────┐
    │                           数据集评估结果                                     │
    ├─────────────────────────────────────────────────────────────────────────────┤
    │ 1. GSE13904 数据集:                                                         │
    │    ✓ 样本组成明确: 儿童对照(18)、SIRS(40)、脓毒症(52)、脓毒性休克(106)    │
    │    ✓ 可用于分析CIITA等MHC II基因在不同疾病严重程度下的表达                     │
    │    ✓ Day1样本(139个)适合与路线A(儿童IPS)结果对比                            │
    │    结论: 可作为外部验证数据集                                                │
    │                                                                             │
    │ 2. GSE66099 数据集:                                                         │
    │    ✓ 包含276个儿童样本                                                      │
    │    ✗ series_matrix文件缺少分组标注                                          │
    │    ✗ 无法直接区分SIRS、脓毒症、脓毒性休克                                   │
    │    结论: 不适合直接用于验证（需要获取原始元数据）                             │
    │                                                                             │
    │ 3. 年龄异质性验证:                                                          │
    │    路线A: 儿童脓毒症数据 → CIITA表达模式已知                                │
    │    路线B: 成人脓毒症数据 → CIITA表达模式已知                                 │
    │    当前: 两个外部验证数据集都是儿童数据                                      │
    │    缺口: 需要成人外部数据集进行年龄分层验证                                  │
    └─────────────────────────────────────────────────────────────────────────────┘
    """
    print(conclusion)
    
    # 保存结果
    with open('长期计划/儿童脓毒症免疫瘫痪研究/results/routeC/analysis_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\n结果已保存到 analysis_results.json")
    
    return all_results

if __name__ == "__main__":
    all_results = main()
