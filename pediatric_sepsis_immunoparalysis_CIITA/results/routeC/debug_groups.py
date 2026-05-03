#!/usr/bin/env python3
import re

# 读取GSE13904样本标题
with open("长期计划/儿童脓毒症免疫瘫痪研究/data/raw/GSE13904_series_matrix.txt", 'r') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            titles = [p.strip('"') for p in parts[1:]]
            print(f"总样本数: {len(titles)}")
            print("\n前10个样本标题:")
            for t in titles[:10]:
                print(f"  {t}")
            
            # 统计各组
            groups = {'Control': 0, 'SIRS': 0, 'Sepsis': 0, 'SepticShock': 0, 'Other': 0}
            day1_count = 0
            
            for title in titles:
                parts = title.split('_')
                if len(parts) >= 3:
                    group_str = parts[1]
                    timepoint = parts[2]
                    
                    if timepoint == 'day1':
                        day1_count += 1
                    
                    if 'Control' in group_str:
                        groups['Control'] += 1
                    elif 'SIRS' in group_str:
                        groups['SIRS'] += 1
                    elif 'Sepsis' in group_str:
                        groups['Sepsis'] += 1
                    elif 'Septic' in group_str:
                        groups['SepticShock'] += 1
                    else:
                        groups['Other'] += 1
                        print(f"  Other: {title}")
            
            print("\n分组统计:")
            for g, c in groups.items():
                print(f"  {g}: {c}")
            print(f"\nDay1样本数: {day1_count}")
            break
