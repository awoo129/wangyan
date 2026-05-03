#!/usr/bin/env python3
import re

# 直接读取文件调试
with open("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt", 'r') as f:
    sample_titles = None
    reading_data = False
    
    for line in f:
        line = line.rstrip('\n')
        
        if line.startswith('!Sample_title'):
            parts = line.split('\t')
            sample_titles = [p.strip('"') for p in parts[1:]]
            print(f"样本标题数: {len(sample_titles)}")
            print(f"前5个标题: {sample_titles[:5]}")
        
        if line.startswith('!series_matrix_table_begin'):
            reading_data = True
            print("开始读取数据")
        
        if reading_data and line and not line.startswith('!'):
            parts = line.split('\t')
            print(f"第一个探针ID: {parts[0]}")
            print(f"表达值数量: {len(parts)-1}")
            break

# 调试分组提取
if sample_titles:
    print("\n" + "="*50)
    print("分组提取调试:")
    
    samples = []
    for i, title in enumerate(sample_titles):
        parts = title.split('_')
        if len(parts) >= 3:
            timepoint = parts[2]
            print(f"样本{i}: {title} -> timepoint={timepoint}")
            if i < 3:
                samples.append({'title': title, 'timepoint': timepoint})
    
    print("\n前3个样本解析结果:")
    for s in samples:
        print(f"  {s}")
