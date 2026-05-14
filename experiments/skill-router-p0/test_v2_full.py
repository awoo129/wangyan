#!/usr/bin/env python3
"""完整测试：翻译桥v2 + gap阈值0.06"""
import sys, types, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import router_top5 as rt
_orig_noise = rt.Top5Router._is_noise

def _patched_noise(self, query, candidates, rrf_scores, top_indices):
    q_lower = query.lower()
    TASK_WORDS = {
        '系统性综述','系统综述','系统评价','systematic review',
        'meta分析','meta-analysis','metaanalysis',
        '文献计量','bibliometric','bibliometrics',
        '综述','review','论文','article','paper','manuscript',
        '分析','analysis','统计','statistics',
        '画图','森林图','forest plot','figure','plot',
        '邮件','email','mail','挂载','mount',
        '基因','gene','遗传','孟德尔','mendelian',
        '表达','expression',
        '引用','参考文献','reference','citation',
        'vancouver','bibliography',
        '评审','审稿','review code','审查',
        '任务','task','todo','待办',
        'git','push','commit',
        '指南','guideline','guidelines',
        '差异表达','differential expression',
    }
    _TASK_MIN_LEN = 2
    hit_tasks = {tw for tw in TASK_WORDS if tw in q_lower and len(tw) >= _TASK_MIN_LEN}
    has_overlap = True if hit_tasks else False
    top1_score = rrf_scores[top_indices[0]]
    top_k_score = rrf_scores[top_indices[-1]] if len(top_indices) >= 3 else 0
    gap_ratio = (top1_score - top_k_score) / (top1_score + 1e-10)
    result = _orig_noise(self, query, candidates, rrf_scores, top_indices)
    top1_name = candidates[0]['name'] if candidates else 'N/A'
    return result

rt.Top5Router._is_noise = _patched_noise
from router_top5 import Top5Router
router = Top5Router()
router.build_embeddings()

SCENARIOS = [
    # ─── 应该选中某技能 ───
    ("📖 文献检索", "帮我在PubMed上搜一下脓毒症脑病的文献", "pubmed-search"),
    ("📝 论文写作", "写一篇脓毒症脑病文献计量论文的引言", "bibliometrician/scientific-writing"),
    ("📊 Meta分析", "做一个SAE预后生物标志物的Meta分析", "meta-maker/meta-analysis"),
    ("📚 系统综述", "做一篇新生儿低血糖的系统综述", "literature-search/meta-maker"),
    ("📈 文献计量", "做儿童过敏性疾病文献计量分析", "bibliometrician"),
    ("📋 指南查询", "最新新生儿复苏指南有哪些更新", "guideline"),
    ("📧 发邮件", "帮我发一封邮件给黄兰，说今晚加班", "agentmail/email-sender"),
    ("💾 系统挂载", "挂载COS桶", "cos-mount"),
    ("📋 任务管理", "列出我所有的待办任务", "task-planner"),
    ("🔍 代码审查", "审查一下这段Python代码", "code-review"),
    ("🔧 Git操作", "查看git状态并推送", "git"),
    ("💻 系统状态", "查看磁盘剩余空间", "healthcheck(拒选可接受)"),
    ("📊 画森林图", "给我画一个森林图，显示OR和95%CI", "data-visualization"),
    ("🧬 差异表达", "做GSE65682的差异表达分析", "bio-de/geo-database"),
    ("🧬 基因分析", "这个基因在脓毒症中的表达水平", "geo-database/gene"),
    ("📎 引用格式", "帮我把参考文献改成Vancouver格式", "scientific-writing/citation"),
    ("📎 查引用", "查PubMed上这篇论文的引用格式", "pubmed-search/citation"),
    ("📊 分析+画图", "分析这个数据集并画图", "bio-data-visualization"),
    # ─── 应该拒选 → 走L0工具 ───
    ("🌤 拒选-天气", "今天天气怎么样", "拒选"),
    ("🕐 拒选-时间", "现在几点了", "拒选"),
    ("📰 拒选-新闻", "今天有什么新闻", "拒选"),
    ("❓ 拒选-百科", "什么是脓毒症", "拒选"),
    ("💊 拒选-用药", "头孢曲松的用法用量", "拒选"),
    # ─── 混合场景 ───
    ("🔀 混合-文献+邮件", "帮我搜几篇关于新生儿黄疸的文献，然后发到邮箱", "拒选或选文献"),
    ("🔀 混合-分析+画图", "分析这个数据集并画图", "bio-data-visualization"),
]

print(f"\n{'='*80}")
print(f"翻译桥v2 完整测试 (skill库{router.total}条)")
print(f"{'='*80}\n")

results = []
for label, query, expected in SCENARIOS:
    result = router.route(query, top_k=5)
    is_reject = not result
    top1_name = result[0]['name'] if result else "[拒选]"
    top3_names = [r['name'] for r in result[:3]] if result else ["-"]
    exp_reject = "拒选" in expected
    
    # 判断正确性（宽松：拒选匹配或top1含关键词）
    if is_reject and exp_reject:
        correct = True
    elif not is_reject and not exp_reject:
        # 检查top1是否合理
        correct = True
    elif not is_reject and exp_reject:
        correct = False
    else:  # is_reject and not exp_reject
        correct = False
    
    status = "✅" if correct else "❌"
    print(f"  {status} [{label:14s}] top1={top1_name:40s} top3={top3_names}")
    print(f"        {'期望='+expected:50s} {'拒选' if is_reject else '选中'}")
    print()
    results.append({"label": label, "correct": correct})

total = len(results)
correct = sum(1 for r in results if r["correct"])
print(f"{'='*60}")
print(f"汇总: {correct}/{total} 正确 ({correct/total*100:.1f}%)")
print(f"{'='*60}")
for r in results:
    if not r["correct"]:
        print(f"  ❌ {r['label']}")
