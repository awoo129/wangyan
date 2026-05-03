#!/usr/bin/env python3
"""
儿童脓毒症免疫瘫痪研究 - GEO数据下载与预处理
Author: Bioinformatic Analysis Pipeline
"""

import os
import sys
import gzip
import logging
import warnings
from pathlib import Path
from datetime import datetime

import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

warnings.filterwarnings('ignore')

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

# ========================
# 配置
# ========================
MAIN_DIR = Path("./长期计划/儿童脓毒症免疫瘫痪研究")
DATA_DIR = MAIN_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
NORM_DIR = DATA_DIR / "normalized"
META_DIR = DATA_DIR / "metadata"
QC_DIR = MAIN_DIR / "results" / "QC"
PREP_DIR = MAIN_DIR / "results" / "preprocessing"
SCRIPT_DIR = MAIN_DIR / "scripts"

# 目标数据集
GSE_IDS = ["GSE26378", "GSE26440", "GSE13904"]

# MHC II类基因列表
MHC_GENES = [
    "HLA-DRA", "HLA-DQB1", "HLA-DQA1", "HLA-DRB1", 
    "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB",
    "CIITA", "CD74", "HLA-DOB", "HLA-DPA2", "HLA-DPB2"
]

# ========================
# 工具函数
# ========================
def create_directories():
    """创建必要的目录结构"""
    for dir_path in [RAW_DIR, NORM_DIR, META_DIR, QC_DIR, PREP_DIR]:
        dir_path.mkdir(parents=True, exist_ok=True)
        log.info(f"创建目录: {dir_path}")

def download_file(url, dest_path, chunk_size=8192):
    """下载文件并显示进度"""
    try:
        response = requests.get(url, stream=True, timeout=120)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(dest_path, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    
        return True
    except Exception as e:
        log.error(f"下载失败 {url}: {e}")
        return False

def parse_series_matrix(matrix_file):
    """解析GEO series matrix文件"""
    log.info(f"解析文件: {matrix_file}")
    
    sample_ids = []
    expression_data = []
    metadata = {}
    current_section = None
    
    with gzip.open(matrix_file, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            # 获取样本ID
            if line.startswith('!Sample_geo_accession'):
                sample_ids = [x.strip('"') for x in line.split('\t')[1:]]
                
            # 获取样本特征
            elif line.startswith('!Sample_') and not line.startswith('!Sample_table_begin'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    key = parts[0].replace('!Sample_', '').replace('_', ' ').title()
                    values = [x.strip('"') for x in parts[1:]]
                    metadata[key] = values
                    
            # 表达数据
            elif not line.startswith('!') and not line.startswith('#'):
                if line and sample_ids:
                    parts = line.split('\t')
                    if len(parts) >= len(sample_ids):
                        try:
                            gene_id = parts[0]
                            values = [float(x) if x != 'NA' else np.nan for x in parts[1:len(sample_ids)+1]]
                            expression_data.append([gene_id] + values)
                        except ValueError:
                            continue
                            
            elif line.startswith('!Sample_table_end'):
                break
    
    # 创建DataFrame
    if expression_data:
        df = pd.DataFrame(expression_data)
        df.columns = ['ID'] + sample_ids
        df = df.set_index('ID')
        df = df.apply(pd.to_numeric, errors='coerce')
        
        sample_metadata = pd.DataFrame(metadata, index=sample_ids).T if metadata else None
        
        return df, sample_metadata
    
    return None, None

def get_gse_download_url(gse_id):
    """获取GEO数据集的下载URL"""
    base_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file"
    return base_url

# ========================
# 数据下载
# ========================
def download_gse_data(gse_id):
    """下载GEO数据集"""
    log.info(f"========== 开始下载数据集: {gse_id} ==========")
    
    matrix_file = RAW_DIR / f"{gse_id}_series_matrix.txt.gz"
    
    if matrix_file.exists():
        log.info(f"{gse_id} 已存在，跳过下载")
        return matrix_file
    
    # 尝试直接下载series matrix
    url = get_gse_download_url(gse_id)
    log.info(f"下载URL: {url}")
    
    if download_file(url, matrix_file):
        log.info(f"{gse_id} 下载成功: {matrix_file}")
        return matrix_file
    else:
        log.error(f"{gse_id} 下载失败")
        return None

# ========================
# 数据预处理
# ========================
def preprocess_expression(expr_df):
    """预处理表达矩阵"""
    log.info("开始数据预处理...")
    
    # 移除缺失值过多的基因
    missing_ratio = expr_df.isnull().sum(axis=1) / expr_df.shape[1]
    valid_genes = missing_ratio[missing_ratio < 0.5].index
    expr_filtered = expr_df.loc[valid_genes]
    
    # 用中位数填充剩余缺失值
    expr_filled = expr_filtered.fillna(expr_filtered.median(axis=1))
    
    # 移除低表达基因（所有样本表达量都低于中位数的基因）
    median_expr = expr_filled.median(axis=1)
    high_expr_genes = median_expr[median_expr > np.percentile(median_expr, 10)].index
    expr_cleaned = expr_filled.loc[high_expr_genes]
    
    log.info(f"预处理后: {expr_cleaned.shape[0]} 基因, {expr_cleaned.shape[1]} 样本")
    
    return expr_cleaned

def map_probes_to_genes(expr_df, annotation_df=None):
    """探针到基因符号映射（如果注释可用）"""
    # 如果ID已经是基因符号，直接返回
    if expr_df.index[0].upper() in [g.upper() for g in MHC_GENES] or \
       not expr_df.index[0].startswith('GPL'):
        return expr_df
    
    # 尝试提取基因符号（有些ID格式: ABC123_at|GENE_SYMBOL）
    def extract_gene_symbol(probe_id):
        if '|' in str(probe_id):
            parts = str(probe_id).split('|')
            for part in parts:
                part = part.strip()
                if part.isalpha() or (len(part) <= 10 and '_' not in part and '-' not in part):
                    return part
        return str(probe_id)
    
    gene_symbols = [extract_gene_symbol(idx) for idx in expr_df.index]
    
    # 聚合重复基因
    expr_df_copy = expr_df.copy()
    expr_df_copy['gene_symbol'] = gene_symbols
    
    expr_aggregated = expr_df_copy.groupby('gene_symbol').mean()
    
    log.info(f"映射后: {expr_aggregated.shape[0]} 唯一基因")
    
    return expr_aggregated

def extract_mhc_genes(expr_df):
    """提取MHC II类基因"""
    log.info("提取MHC II类基因...")
    
    # 标准化基因符号
    expr_df.index = expr_df.index.str.upper()
    
    found_genes = []
    not_found_genes = []
    
    for gene in MHC_GENES:
        gene_upper = gene.upper()
        if gene_upper in expr_df.index:
            found_genes.append(gene)
        elif gene.upper() in ['CIITA', 'CD74']:  # 这些不是GPL570探针
            not_found_genes.append(gene)
        else:
            # 尝试模糊匹配
            matches = [idx for idx in expr_df.index if gene.upper() in idx.upper()]
            if matches:
                found_genes.append(gene)
            else:
                not_found_genes.append(gene)
    
    if found_genes:
        mhc_expr = expr_df.loc[[g.upper() for g in found_genes]]
        log.info(f"成功提取 {len(found_genes)} 个MHC II类基因")
    else:
        mhc_expr = None
        log.warning("未找到任何MHC II类基因!")
    
    return {
        'expression': mhc_expr,
        'found_genes': found_genes,
        'not_found_genes': not_found_genes
    }

# ========================
# 数据质控
# ========================
def perform_qc(expr_df, metadata, gse_id, sample_info=None):
    """执行数据质控"""
    log.info(f"========== 质控分析: {gse_id} ==========")
    
    qc_gse_dir = QC_DIR / gse_id
    qc_gse_dir.mkdir(exist_ok=True)
    
    results = {}
    
    # 1. 箱线图
    plt.figure(figsize=(12, 6))
    expr_melted = expr_df.iloc[:, :min(50, expr_df.shape[1])].T  # 限制50个样本
    plt.boxplot(expr_melted.values, showfliers=False)
    plt.title(f'{gse_id} - 样本表达分布 (箱线图)')
    plt.xlabel('Samples')
    plt.ylabel('Expression')
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(qc_gse_dir / f'{gse_id}_boxplot.png', dpi=150)
    plt.close()
    results['boxplot'] = str(qc_gse_dir / f'{gse_id}_boxplot.png')
    
    # 2. 密度图
    plt.figure(figsize=(10, 6))
    for i, col in enumerate(expr_df.columns[:20]):  # 限制20个样本
        density = stats.gaussian_kde(expr_df[col].dropna())
        xs = np.linspace(expr_df[col].min(), expr_df[col].max(), 200)
        plt.plot(xs, density(xs), alpha=0.7, label=f'S{i+1}')
    plt.title(f'{gse_id} - 表达密度分布')
    plt.xlabel('Expression')
    plt.ylabel('Density')
    plt.legend().set_visible(False)
    plt.tight_layout()
    plt.savefig(qc_gse_dir / f'{gse_id}_density.png', dpi=150)
    plt.close()
    results['density'] = str(qc_gse_dir / f'{gse_id}_density.png')
    
    # 3. PCA分析
    try:
        scaler = StandardScaler()
        expr_scaled = scaler.fit_transform(expr_df.T)
        
        pca = PCA(n_components=min(10, expr_df.shape[1]-1))
        pca_result = pca.fit_transform(expr_scaled)
        
        plt.figure(figsize=(10, 8))
        plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)
        for i, txt in enumerate(expr_df.columns):
            plt.annotate(txt[:10], (pca_result[i, 0], pca_result[i, 1]), fontsize=8)
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        plt.title(f'{gse_id} - PCA分析')
        plt.tight_layout()
        plt.savefig(qc_gse_dir / f'{gse_id}_pca.png', dpi=150)
        plt.close()
        
        results['pca'] = {
            'explained_variance': pca.explained_variance_ratio_[:5].tolist(),
            'plot': str(qc_gse_dir / f'{gse_id}_pca.png')
        }
    except Exception as e:
        log.warning(f"PCA分析失败: {e}")
        results['pca'] = None
    
    # 4. 层次聚类
    try:
        # 计算距离矩阵
        from scipy.spatial.distance import pdist, squareform
        dist_matrix = pdist(expr_df.T, metric='euclidean')
        linkage_matrix = linkage(dist_matrix, method='ward')
        
        plt.figure(figsize=(12, 6))
        dendrogram(linkage_matrix, labels=expr_df.columns.tolist(), leaf_rotation=90)
        plt.title(f'{gse_id} - 样本聚类树')
        plt.xlabel('Samples')
        plt.ylabel('Distance')
        plt.tight_layout()
        plt.savefig(qc_gse_dir / f'{gse_id}_dendrogram.png', dpi=150)
        plt.close()
        results['dendrogram'] = str(qc_gse_dir / f'{gse_id}_dendrogram.png')
    except Exception as e:
        log.warning(f"聚类分析失败: {e}")
        results['dendrogram'] = None
    
    # 5. 质控统计
    qc_stats = {
        'total_samples': int(expr_df.shape[1]),
        'total_genes': int(expr_df.shape[0]),
        'mean_expression': float(expr_df.values.mean()),
        'median_expression': float(np.median(expr_df.values)),
        'std_expression': float(expr_df.values.std()),
        'min_expression': float(expr_df.values.min()),
        'max_expression': float(expr_df.values.max()),
        'missing_rate': float(expr_df.isnull().sum().sum() / expr_df.size)
    }
    
    results['stats'] = qc_stats
    
    # 保存统计结果
    stats_df = pd.DataFrame([qc_stats])
    stats_df.to_csv(qc_gse_dir / f'{gse_id}_qc_stats.csv', index=False)
    
    return results

# ========================
# MHC基因分析
# ========================
def analyze_mhc_genes(mhc_results, gse_ids):
    """分析MHC基因检测结果"""
    log.info("========== MHC基因分析 ==========")
    
    # 创建汇总表
    summary_data = []
    
    for gse_id in gse_ids:
        if gse_id in mhc_results:
            result = mhc_results[gse_id]
            summary_data.append({
                'GSE_ID': gse_id,
                'Found_Genes': len(result['found_genes']),
                'Found_List': ', '.join(result['found_genes']),
                'Not_Found_List': ', '.join(result['not_found_genes']) if result['not_found_genes'] else 'N/A'
            })
        else:
            summary_data.append({
                'GSE_ID': gse_id,
                'Found_Genes': 0,
                'Found_List': 'N/A',
                'Not_Found_List': 'N/A'
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(PREP_DIR / 'mhc_gene_detection.csv', index=False)
    
    # HLA-DRA状态
    hla_dra_status = {}
    for gse_id in gse_ids:
        if gse_id in mhc_results:
            hla_dra_status[gse_id] = 'HLA-DRA' in mhc_results[gse_id]['found_genes']
        else:
            hla_dra_status[gse_id] = False
    
    log.info("\nHLA-DRA检测结果:")
    for gse_id, present in hla_dra_status.items():
        status = "✓ 存在" if present else "✗ 不存在"
        log.info(f"  {gse_id}: {status}")
    
    return summary_df, hla_dra_status

# ========================
# 生成报告
# ========================
def generate_report(all_data, mhc_results, qc_results, hla_dra_status, gse_ids):
    """生成分析报告"""
    log.info("========== 生成分析报告 ==========")
    
    # 数据集汇总
    summary_data = []
    for gse_id in gse_ids:
        if gse_id in all_data:
            expr = all_data[gse_id]['expression']
            meta = all_data[gse_id]['metadata']
            
            normal_count = (meta['group'] == 'Normal').sum() if 'group' in meta.columns else 0
            sepsis_count = (meta['group'] == 'Sepsis').sum() if 'group' in meta.columns else 0
            
            summary_data.append({
                'GSE_ID': gse_id,
                'Platform': all_data[gse_id].get('platform', 'GPL570'),
                'Total_Samples': expr.shape[1],
                'Normal_Samples': int(normal_count),
                'Sepsis_Samples': int(sepsis_count),
                'Total_Genes': expr.shape[0],
                'MHC_Genes_Found': len(mhc_results[gse_id]['found_genes']) if gse_id in mhc_results else 0,
                'HLA_DRA_Present': 'Yes' if hla_dra_status.get(gse_id, False) else 'No'
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(PREP_DIR / 'dataset_summary.csv', index=False)
    
    # 生成Markdown报告
    report = f"""# 儿童脓毒症免疫瘫痪研究 - GEO数据预处理报告

## 1. 研究概述

本分析旨在筛选儿童脓毒症免疫瘫痪的分子标志物，采用GEO公共数据库中的芯片数据进行生物信息学分析。

### 目标基因
关键基因：HLA-DRA, HLA-DQB1, HLA-DQA1, CIITA, CD74
扩展基因：HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DMA, HLA-DMB

## 2. 数据集汇总

| 数据集 | 平台 | 总样本 | 正常对照 | 脓毒症 | 基因数 | MHC基因 | HLA-DRA |
|--------|------|--------|----------|--------|--------|---------|---------|
"""
    
    for _, row in summary_df.iterrows():
        report += f"| {row['GSE_ID']} | {row['Platform']} | {row['Total_Samples']} | "
        report += f"{row['Normal_Samples']} | {row['Sepsis_Samples']} | "
        report += f"{row['Total_Genes']} | {row['MHC_Genes_Found']} | {row['HLA_DRA_Present']} |\n"
    
    report += """
## 3. MHC II类基因检测结果

### 检测到的基因
"""
    
    for gse_id in gse_ids:
        if gse_id in mhc_results and mhc_results[gse_id]['found_genes']:
            report += f"\n**{gse_id}**:\n"
            report += f"- 找到: {', '.join(mhc_results[gse_id]['found_genes'])}\n"
            if mhc_results[gse_id]['not_found_genes']:
                report += f"- 未找到: {', '.join(mhc_results[gse_id]['not_found_genes'])}\n"
    
    report += """
## 4. 质控结果

### 4.1 质控指标
"""
    
    for gse_id in gse_ids:
        if gse_id in qc_results and qc_results[gse_id]:
            stats = qc_results[gse_id].get('stats', {})
            report += f"\n**{gse_id}**:\n"
            report += f"- 样本数: {stats.get('total_samples', 'N/A')}\n"
            report += f"- 基因数: {stats.get('total_genes', 'N/A')}\n"
            report += f"- 平均表达: {stats.get('mean_expression', 0):.2f}\n"
            report += f"- 中位表达: {stats.get('median_expression', 0):.2f}\n"
    
    report += """
### 4.2 质控图表

质控图表保存于: `results/QC/{GSE_ID}/`

- Boxplot: 样本表达分布箱线图
- Density: 表达密度分布图
- PCA: 主成分分析图
- Dendrogram: 层次聚类树

## 5. HLA-DRA验证结果

| 数据集 | HLA-DRA状态 |
|--------|-------------|
"""
    
    for gse_id, present in hla_dra_status.items():
        status = "✓ 存在" if present else "✗ 不存在"
        report += f"| {gse_id} | {status} |\n"
    
    # 检查是否所有数据集都包含HLA-DRA
    all_present = all(hla_dra_status.values())
    if all_present:
        report += """
## 6. 结论与建议

✓ **HLA-DRA在所有GPL570数据集中均存在**，可作为免疫瘫痪评分的核心指标。

### 下游分析建议

1. **差异表达分析**：比较正常对照组 vs 脓毒症组的HLA-DRA表达差异
2. **免疫评分构建**：基于HLA-DRA等MHC II类基因表达构建免疫瘫痪评分(IPS)
3. **预后分析**：结合临床预后信息评估IPS的预测价值
4. **独立验证**：使用GSE145227（lncRNA数据）进行交叉验证

## 7. 文件清单

### 原始数据
- `data/raw/{GSE_ID}_series_matrix.txt.gz`

### 处理后数据
- `data/normalized/{GSE_ID}_expression.csv` - 标准化表达矩阵
- `data/normalized/{GSE_ID}_mhc_expression.csv` - MHC基因表达矩阵
- `data/metadata/{GSE_ID}_metadata.csv` - 样本元数据

### 质控结果
- `results/QC/{GSE_ID}/` - 各数据集质控图表
- `results/QC/{GSE_ID}_qc_stats.csv` - 质控统计指标

### 汇总文件
- `results/preprocessing/dataset_summary.csv` - 数据集汇总表
- `results/preprocessing/mhc_gene_detection.csv` - MHC基因检测结果

---
*报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""
    
    # 保存报告
    report_file = PREP_DIR / 'preprocessing_report.md'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    log.info(f"报告已保存: {report_file}")
    
    return summary_df, report

# ========================
# 主程序
# ========================
def main():
    log.info("=" * 60)
    log.info("儿童脓毒症免疫瘫痪研究 - GEO数据预处理")
    log.info("=" * 60)
    
    # 创建目录
    create_directories()
    
    # 存储结果
    all_data = {}
    mhc_results = {}
    qc_results = {}
    
    # 处理每个数据集
    for gse_id in GSE_IDS:
        log.info(f"\n{'='*60}")
        log.info(f"处理数据集: {gse_id}")
        log.info(f"{'='*60}")
        
        # 1. 下载数据
        matrix_file = download_gse_data(gse_id)
        
        if matrix_file is None or not matrix_file.exists():
            log.error(f"{gse_id} 下载失败，跳过")
            continue
        
        # 2. 解析数据
        expr_df, sample_metadata = parse_series_matrix(matrix_file)
        
        if expr_df is None:
            log.error(f"{gse_id} 解析失败，跳过")
            continue
        
        log.info(f"{gse_id}: {expr_df.shape[0]} 探针, {expr_df.shape[1]} 样本")
        
        # 3. 预处理
        expr_cleaned = preprocess_expression(expr_df)
        expr_mapped = map_probes_to_genes(expr_cleaned)
        
        # 4. 保存处理后数据
        expr_mapped.to_csv(NORM_DIR / f'{gse_id}_expression.csv')
        log.info(f"保存表达矩阵: {NORM_DIR / f'{gse_id}_expression.csv'}")
        
        # 5. 提取样本分组
        if sample_metadata is not None:
            groups = []
            for col in sample_metadata.columns:
                if any(kw in col.lower() for kw in ['disease', 'source', 'status', 'group', 'condition']):
                    for val in sample_metadata[col]:
                        if 'normal' in str(val).lower() or 'healthy' in str(val).lower():
                            groups.append('Normal')
                        elif 'sepsis' in str(val).lower() or 'septic' in str(val).lower():
                            groups.append('Sepsis')
                        else:
                            groups.append('Unknown')
                    break
            
            if len(groups) != len(sample_metadata.columns):
                groups = ['Unknown'] * len(sample_metadata.columns)
        else:
            groups = ['Unknown'] * expr_df.shape[1]
        
        # 创建元数据
        metadata = pd.DataFrame({
            'sample_id': expr_df.columns.tolist(),
            'group': groups,
            'gse_id': gse_id,
            'platform': 'GPL570'
        })
        metadata.to_csv(META_DIR / f'{gse_id}_metadata.csv', index=False)
        
        # 6. 提取MHC基因
        mhc_data = extract_mhc_genes(expr_mapped)
        mhc_results[gse_id] = mhc_data
        
        if mhc_data['expression'] is not None:
            mhc_data['expression'].to_csv(NORM_DIR / f'{gse_id}_mhc_expression.csv')
        
        # 7. 质控分析
        qc_result = perform_qc(expr_mapped, metadata, gse_id)
        qc_results[gse_id] = qc_result
        
        # 保存完整数据
        all_data[gse_id] = {
            'expression': expr_mapped,
            'metadata': metadata,
            'platform': 'GPL570'
        }
        
        # 清理内存
        del expr_df, sample_metadata, expr_cleaned, expr_mapped
        import gc
        gc.collect()
    
    # MHC基因分析
    mhc_summary_df, hla_dra_status = analyze_mhc_genes(mhc_results, GSE_IDS)
    
    # 生成报告
    summary_df, report = generate_report(all_data, mhc_results, qc_results, hla_dra_status, GSE_IDS)
    
    # 保存所有结果
    log.info("\n保存完整结果...")
    
    # 保存为pickle格式
    import pickle
    with open(NORM_DIR / 'all_data.pkl', 'wb') as f:
        pickle.dump(all_data, f)
    
    with open(NORM_DIR / 'mhc_results.pkl', 'wb') as f:
        pickle.dump(mhc_results, f)
    
    with open(QC_DIR / 'qc_results.pkl', 'wb') as f:
        pickle.dump(qc_results, f)
    
    log.info("\n" + "=" * 60)
    log.info("数据预处理完成!")
    log.info("=" * 60)
    
    # 打印汇总
    print("\n" + "=" * 60)
    print("数据集汇总:")
    print("=" * 60)
    print(summary_df.to_string(index=False))
    
    print("\n" + "=" * 60)
    print("MHC基因检测结果:")
    print("=" * 60)
    print(mhc_summary_df.to_string(index=False))
    
    print("\n" + "=" * 60)
    print("HLA-DRA验证结果:")
    print("=" * 60)
    for gse_id, present in hla_dra_status.items():
        status = "✓ 存在" if present else "✗ 不存在"
        print(f"  {gse_id}: {status}")
    
    return all_data, mhc_results, qc_results, summary_df

if __name__ == "__main__":
    results = main()
