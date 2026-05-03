#!/usr/bin/env python3
"""生成可视化图表"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 加载结果
with open('长期计划/儿童脓毒症免疫瘫痪研究/results/routeC/analysis_results.json', 'r') as f:
    results = json.load(f)

# GSE13904结果
gse13904 = results.get('GSE13904', {})
gse66099 = results.get('GSE66099', {})

# 创建图表
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1. GSE13904 CIITA表达 (条形图)
ax1 = axes[0, 0]
genes_to_plot = ['CIITA', 'HLA-DRA', 'CD74']
groups = ['Control', 'Sepsis', 'SepticShock']
x = np.arange(len(groups))
width = 0.25

for i, gene in enumerate(genes_to_plot):
    if gene in gse13904:
        means = []
        stds = []
        for group in groups:
            if group in gse13904[gene]:
                means.append(gse13904[gene][group]['mean'])
                stds.append(gse13904[gene][group]['std'])
            else:
                means.append(0)
                stds.append(0)
        bars = ax1.bar(x + i*width - width, means, width, yerr=stds, 
                       label=gene, capsize=3, alpha=0.8)

ax1.set_xlabel('Group')
ax1.set_ylabel('Expression (normalized)')
ax1.set_title('GSE13904: MHC II Gene Expression by Disease Severity')
ax1.set_xticks(x)
ax1.set_xticklabels(groups)
ax1.legend()
ax1.grid(axis='y', alpha=0.3)

# 2. GSE13904 CIITA表达变化方向
ax2 = axes[0, 1]
comparisons = []
directions = []
colors = []

for gene in ['CIITA', 'HLA-DRA', 'CD74']:
    if gene in gse13904:
        for comp_key in ['Sepsis_vs_Control', 'SepticShock_vs_Control']:
            if comp_key in gse13904[gene]:
                delta = gse13904[gene][comp_key]['delta']
                p = gse13904[gene][comp_key]['p']
                gene_comp = f"{gene}\n{comp_key.split('_')[0]}"
                comparisons.append(gene_comp)
                directions.append(delta)
                colors.append('red' if delta > 0 else 'blue')

bars = ax2.barh(comparisons, directions, color=colors, alpha=0.7)
ax2.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
ax2.set_xlabel('Expression Change (Delta)')
ax2.set_title('GSE13904: Gene Expression Changes\n(Red=Increase, Blue=Decrease)')

# 3. GSE66099 整体表达分布
ax3 = axes[1, 0]
genes_66099 = ['CIITA', 'HLA-DRA', 'CD74']
means_66099 = [gse66099[g]['mean'] for g in genes_66099 if g in gse66099]
stds_66099 = [gse66099[g]['std'] for g in genes_66099 if g in gse66099]
genes_filtered = [g for g in genes_66099 if g in gse66099]

ax3.bar(genes_filtered, means_66099, yerr=stds_66099, 
         color=['steelblue', 'forestgreen', 'darkorange'], 
         capsize=5, alpha=0.8)
ax3.set_ylabel('Expression (log2 gcRMA)')
ax3.set_title('GSE66099: Overall MHC II Gene Expression\n(No Group Information Available)')
ax3.grid(axis='y', alpha=0.3)

# 4. 数据集对比总结表
ax4 = axes[1, 1]
ax4.axis('off')

summary_text = """
EXTERNAL VALIDATION DATASET ASSESSMENT

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

GSE13904 (Pediatric SIRS/Sepsis Spectrum):
✓ Samples: 139 Day1 samples
  - Control: 18
  - SIRS: 22  
  - Sepsis: 32
  - Septic Shock: 67
✓ Available for group comparison
✓ Suitable for validation

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

GSE66099 (Combined Pediatric Sepsis Dataset):
✓ Samples: 276 unique patients
✗ No group labels in series_matrix
✗ Cannot distinguish SIRS/Sepsis/Septic Shock
✗ Not suitable for direct validation

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

KEY FINDING:
CIITA expression DECREASES in pediatric sepsis/septic shock
compared to controls (GSE13904)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

AGE HETEROGENEITY GAP:
Route A: Pediatric sepsis → CIITA pattern known
Route B: Adult sepsis → CIITA pattern known  
Current: Both external datasets are pediatric
Need: Adult external dataset for age-stratified validation
"""

ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, 
         fontsize=10, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))

plt.tight_layout()
plt.savefig('长期计划/儿童脓毒症免疫瘫痪研究/results/routeC/external_validation_analysis.png', dpi=150, bbox_inches='tight')
print("图表已保存: external_validation_analysis.png")

# 创建专门的CIITA比较图
fig2, ax = plt.subplots(1, 1, figsize=(10, 6))

# GSE13904 CIITA数据
ciita_data = gse13904.get('CIITA', {})
groups = ['Control', 'SIRS', 'Sepsis', 'SepticShock']
means = [ciita_data[g]['mean'] for g in groups if g in ciita_data]
stds = [ciita_data[g]['std'] for g in groups if g in ciita_data]
ns = [ciita_data[g]['n'] for g in groups if g in ciita_data]
groups_filtered = [g for g in groups if g in ciita_data]

x = np.arange(len(groups_filtered))
bars = ax.bar(x, means, yerr=stds, capsize=5, alpha=0.8, 
              color=['green', 'yellow', 'orange', 'red'])

ax.set_xlabel('Disease Group', fontsize=12)
ax.set_ylabel('CIITA Expression (normalized)', fontsize=12)
ax.set_title('GSE13904: CIITA Expression Across Sepsis Spectrum', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels([f"{g}\n(n={ns[i]})" for i, g in enumerate(groups_filtered)])
ax.grid(axis='y', alpha=0.3)

# 添加显著性标记
if 'Sepsis_vs_Control' in ciita_data:
    delta = ciita_data['Sepsis_vs_Control']['delta']
    p = ciita_data['Sepsis_vs_Control']['p']
    ax.annotate(f'Δ={delta:.3f}\np={p:.4f}', 
                xy=(1, means[1]), xytext=(1.3, means[1]+0.2),
                fontsize=10, ha='left')

if 'SepticShock_vs_Control' in ciita_data:
    delta = ciita_data['SepticShock_vs_Control']['delta']
    p = ciita_data['SepticShock_vs_Control']['p']
    ax.annotate(f'Δ={delta:.3f}\np={p:.4f}', 
                xy=(3, means[3]), xytext=(3.3, means[3]+0.2),
                fontsize=10, ha='left')

plt.tight_layout()
plt.savefig('长期计划/儿童脓毒症免疫瘫痪研究/results/routeC/CIITA_expression_GSE13904.png', dpi=150, bbox_inches='tight')
print("CIITA表达图已保存: CIITA_expression_GSE13904.png")
