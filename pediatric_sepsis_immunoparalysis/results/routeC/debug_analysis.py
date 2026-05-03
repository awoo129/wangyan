#!/usr/bin/env python3
import re
import gzip
import numpy as np

PROBE_MAP = {
    'CIITA': ['208981_at', '208983_s_at', '214523_at'],
    'HLA-DRA': ['212828_s_at', '214486_at', '201137_s_at'],
    'HLA-DPA1': ['210044_s_at', '209823_at'],
    'CD74': ['200664_at', '212068_at', '201009_s_at']
}

# 读取GSE13904
print("读取GSE13904...")
with open("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt", 'r') as f:
    probe_ids = []
    expression_data = []
    sample_info = {}
    reading_data = False
    
    for line in f:
        line = line.rstrip('\n')
        
        if not reading_data and line.startswith('!') and not line.startswith('!series_matrix_table_begin'):
            if '\t' in line:
                parts = line.split('\t', 1)
                if len(parts) >= 2:
                    key_match = re.match(r'!(\w+)', parts[0])
                    if key_match:
                        key = key_match.group(1)
                        values = parts[1].split('\t')
                        sample_info[key] = [v.strip('"') for v in values]
        
        if line.startswith('!series_matrix_table_begin'):
            reading_data = True
            continue
        
        if line.startswith('!series_matrix_table_end'):
            break
        
        if reading_data and line and '\t' in line:
            parts = line.split('\t')
            if parts:
                probe_id = parts[0].strip('"')
                try:
                    expr_values = [float(v.strip('"')) for v in parts[1:]]
                    probe_ids.append(probe_id)
                    expression_data.append(expr_values)
                except ValueError:
                    continue

print(f"探针数: {len(probe_ids)}, 样本数: {len(expression_data[0]) if expression_data else 0}")

# 检查目标探针
for gene, probes in PROBE_MAP.items():
    found_probes = [p for p in probes if p in probe_ids]
    print(f"\n{gene}: 目标探针{probes}, 找到: {found_probes}")

# 提取CIITA表达
ciita_expr = None
for probe in PROBE_MAP['CIITA']:
    if probe in probe_ids:
        idx = probe_ids.index(probe)
        ciita_expr = np.array(expression_data[idx])
        print(f"\nCIITA ({probe}) 找到，索引: {idx}")
        print(f"表达值示例: {ciita_expr[:5]}")
        print(f"CIITA mean: {np.mean(ciita_expr):.4f}")
        break

if ciita_expr is None:
    print("CIITA探针未找到!")
