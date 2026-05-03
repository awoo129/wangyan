#!/usr/bin/env python3
import re

# 读取样本标题
with open("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt", 'r') as f:
    sample_titles = None
    for line in f:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            sample_titles = [p.strip('"') for p in parts[1:]]
            break

print(f"样本标题数: {len(sample_titles)}")

# 解析分组和时间点
day1_count = 0
for title in sample_titles:
    parts = title.split('_')
    if len(parts) >= 3:
        timepoint = parts[2]
        if timepoint == 'day1':
            day1_count += 1

print(f"Day1样本数: {day1_count}")

# 检查问题
for i, title in enumerate(sample_titles[:5]):
    parts = title.split('_')
    print(f"样本{i}: {title}")
    print(f"  分隔后: {parts}")
    print(f"  parts[2]: {parts[2] if len(parts) > 2 else 'N/A'}")
