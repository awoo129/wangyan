# 🧠 Skill Router v2 — 智能技能路由系统

> 977个可用技能的智能路由引擎，96%正确率

## 问题

OpenClaw 技能库从几十个膨胀到 977 个后，靠 LLM 自己逐条审视所有 skill 变得不现实——token 成本高、上下文窗口不够、延迟大。

## 方案

**两阶段路由**：Embedding 粗筛 → LLM 精排（本级仅包含粗筛阶段）

```
用户 Query → BM25 + Embedding 双路检索 → RRF 融合 → top5 候选 → _is_noise 拒选 → LLM 精排
```

## ✨ v2 优化（2026-05-14）

### 中英翻译桥 v2

skill 库全英文化（OpenClaw ClawHub 镜像），但用户 query 是中文。all-MiniLM-L6-v2 对中文 embedding 质量差，导致匹配全是模糊近似。

**翻译桥方案**：规则替换中→英关键词，生成纯英文 query，与中文原版按 `0.3zh + 0.7en` 加权融合。

```
"做一篇新生儿低血糖的系统综述"
→ "systematic review meta analysis literature synthesis conduct summary"
→ 正确匹配 literature-search ✅
```

包含 50+ 中→英专业术语映射（系统综述→systematic review、差异表达→differential expression、森林图→forest plot...）

### 拒选防线

| 防线 | 作用 | 触发条件 |
|------|------|---------|
| 0 — 信息查询 | 天气/百科/用药等直接L0工具 | INFO_SEEK 词命中且无强技能词 |
| A — 任务重叠 | query 中有可识别的任务词 | TASK_WORDS 命中（综述/分析/画图...） |
| B — 分数阈值 | 候选与 query 的相关度太低 | RRF top1 < 0.008 或 gap < 0.06 |

双重触发（无任务重叠 + 分数低）才拒选，避免误拦。

## 🎯 测试结果

25个炎哥常用场景测试：

| 场景 | 改前(纯中文) | 改后(v2) |
|------|-------------|----------|
| 📖 PubMed搜文献 | ✅ pubmed-search | ✅ pubmed-search |
| 📝 写文献计量论文 | ❌ humanizer-zh | ✅ bibliometrician |
| 📚 做系统综述 | ❌ patiently-ai | ✅ literature-search |
| 📧 发邮件 | ❌ monday | ✅ agentmail/email-sender |
| 💾 挂载COS | ✅ cos-mount | ✅ cos-mount |
| 📊 画森林图 | ❌ care-coordination | ✅ data-visualization-biomedical |
| 📎 改Vancouver格式 | ❌ humanizer-zh | ✅ scientific-writing |
| 📊 分析数据画图 | ❌ monday | ✅ bio-data-visualization |
| 🧬 差异表达分析 | ❌ | ✅ bio-de-visualization |
| 🌤 天气/百科查询 | ✅ 拒选 | ✅ 拒选 |
| **整体正确率** | **20.8%** | **96%** |

## 🚀 快速体验

```bash
cd experiments/skill-router-p0
python3 -c "
from router_top5 import Top5Router
router = Top5Router()
print(router.route_formatted('做一篇新生儿低血糖的系统综述'))
"
```

## 📁 文件结构

| 文件 | 说明 |
|------|------|
| `router_top5.py` | 核心路由引擎（BM25 + 翻译桥 + Embedding + RRF + 拒选） |
| `schemas.py` | SkillProfile 数据模型 |
| `tool_catalog.json` | L0 工具目录 |
| `tool_first_router.py` | L0 工具优先路由 |
| `reranker.py` | Cross-encoder 精排（预留） |
| `skill_index_v2.json` | 977个 skill 索引 |
| `skill_index_v2_optimized.json` | 优化后的 skill 索引 |
| `skill_desc_optimize.py` | skill 描述优化工具 |
| `test_v2_full.py` | 25场景全覆盖测试 |
| `p0_summary.py` | Project 0 实验总结 |

## 🔮 待办

- [ ] 接入 OpenClaw runtime（当前为独立原型）
- [ ] Cross-encoder 精排整合
- [ ] 多技能编排（orchestration）
- [ ] 用户反馈闭环学习

---

**作者**: 炎哥 🔥 | **仓库**: [awoo129/wangyan](https://github.com/awoo129/wangyan)
