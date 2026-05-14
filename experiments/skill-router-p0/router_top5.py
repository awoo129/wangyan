#!/usr/bin/env python3
"""SkillRouter v2 — Python embedding 粗筛 top5，返回丰富上下文给 LLM 精排"""
import json
import sys
import re
from pathlib import Path
import numpy as np

V2_DIR = Path(__file__).parent
WORKSPACE = V2_DIR.parent

# 尝试加载 embedding 模型
_EMBEDDING_MODEL = None


def _load_embedding_model():
    global _EMBEDDING_MODEL
    if _EMBEDDING_MODEL is not None:
        return _EMBEDDING_MODEL
    try:
        import os
        os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"
        from sentence_transformers import SentenceTransformer
        _EMBEDDING_MODEL = SentenceTransformer("all-MiniLM-L6-v2")
        return _EMBEDDING_MODEL
    except Exception as e:
        print(f"[router_top5] embedding 模型加载失败: {e}", file=sys.stderr)
        return None


class Top5Router:
    """embedding 检索 → top5 候选（含丰富画像）"""

    def __init__(self, index_path: Path = None):
        if index_path is None:
            index_path = V2_DIR / "skill_index_v2.json"

        # 加载索引
        self.data = json.loads(index_path.read_text(encoding="utf-8"))
        self.skills: list[dict] = self.data["skills"]
        self.total = self.data["total"]

        # 加载 embeddings（如果存在）
        emb_path = index_path.with_suffix(".embeddings.npy")
        self.embeddings: np.ndarray | None = None
        if emb_path.exists():
            self.embeddings = np.load(emb_path)

        self._model = None

    @property
    def model(self):
        if self._model is None:
            self._model = _load_embedding_model()
        return self._model

    def build_embeddings(self, force: bool = False):
        """为所有 skill 生成 embeddings（首次运行需要）"""
        emb_path = V2_DIR / "skill_index_v2.embeddings.npy"
        if emb_path.exists() and not force:
            self.embeddings = np.load(emb_path)
            return

        model = self.model
        if model is None:
            print("[router_top5] embedding 模型不可用，无法生成 embeddings", file=sys.stderr)
            return

        # 构建 embedding 文本（仅用 compact 字段避免 OOM；body 留给 LLM 精排阶段）
        texts = []
        for s in self.skills:
            parts = [
                f"名称: {s['name']}",
                f"描述: {s.get('description', '')}",
            ]
            caps = s.get('capabilities', [])
            if caps:
                parts.append(f"能力: {', '.join(caps)}")
            t = '\n'.join(parts)
            texts.append(t)

        print(f"[router_top5] 正在为 {len(texts)} 个 skill 生成 embeddings...", file=sys.stderr)

        # 小批次编码避免 OOM（每次 50 个，编码完立即释放）
        import gc
        batch_size = 50
        all_embeddings = []
        for i in range(0, len(texts), batch_size):
            batch = texts[i:i + batch_size]
            batch_emb = model.encode(batch, show_progress_bar=False)
            all_embeddings.append(batch_emb)
            gc.collect()
            if (i // batch_size) % 5 == 0:
                print(f"  ... {i}/{len(texts)}", file=sys.stderr)

        self.embeddings = np.concatenate(all_embeddings, axis=0)
        np.save(emb_path, self.embeddings)
        print(f"[router_top5] embeddings 已保存到 {emb_path} (shape: {self.embeddings.shape})", file=sys.stderr)

    def _bm25_search(self, query: str, top_k: int = 20) -> list[int]:
        """BM25 关键词检索（name + description + capabilities 作为文档）"""
        import math
        from collections import defaultdict

        # 分词（简单按空格+中文单字切分）
        import re
        def tokenize(text):
            tokens = re.findall(r'[\u4e00-\u9fff]|[a-zA-Z0-9]+', text.lower())
            return tokens

        query_tokens = tokenize(query)
        if not query_tokens:
            return []

        # 预计算文档
        docs = []
        for s in self.skills:
            doc = f"{s['name']} {s.get('description','')} {' '.join(s.get('capabilities',[]))}"
            docs.append(tokenize(doc))

        # IDF
        N = len(docs)
        df = defaultdict(int)
        for doc in docs:
            for t in set(doc):
                df[t] += 1

        k1, b = 1.5, 0.75
        avgdl = sum(len(d) for d in docs) / max(N, 1)

        scores = []
        for i, doc in enumerate(docs):
            dl = len(doc)
            score = 0.0
            tf = defaultdict(int)
            for t in doc:
                tf[t] += 1
            for t in query_tokens:
                if t in tf:
                    idf = math.log((N - df[t] + 0.5) / (df[t] + 0.5) + 1)
                    term_score = idf * (tf[t] * (k1 + 1)) / (tf[t] + k1 * (1 - b + b * dl / avgdl))
                    score += term_score
            scores.append((i, score))

        scores.sort(key=lambda x: -x[1])
        return [idx for idx, s in scores[:top_k] if s > 0]

    # ─── 中英翻译桥 v2 ────────────────────────────────────────────────
    # 思路：识别query中的任务关键词，翻译为纯英文query
    # 用 "; " 分隔不同意图片段，让 embedding 向量方向更精确
    _ZH2EN_MAP = {
        # 通用动作
        '搜索': 'search', '检索': 'search', '查找': 'search',
        '写一篇': 'write paper', '写一': 'write paper',
        '撰写': 'write', '写作': 'draft paper',
        '做一篇': 'conduct', '做一': 'conduct',
        '分析': 'analyze',
        '画图': 'draw visualization', '画一个': 'draw plot',
        '画': 'draw',
        '发邮件': 'send email', '发一封': 'send email',
        '挂载': 'mount',
        '查看': 'check',
        '列出': 'list',
        '审查': 'review code',
        '评审': 'review evaluate',
        '审稿': 'peer review',
        '改成': 'convert to', '改为': 'convert to',
        '查': 'lookup search',
        # 研究类型
        '系统综述': 'systematic review meta analysis literature synthesis',
        '系统评价': 'systematic review evidence evaluation',
        'meta分析': 'meta analysis statistical pooling',
        'Meta分析': 'meta analysis statistical pooling',
        '文献计量': 'bibliometric analysis science mapping',
        '文献检索': 'literature search database',
        '指南': 'clinical guideline',
        '综述': 'literature review summary',
        '论文': 'academic paper manuscript article',
        '写论文': 'write academic paper',
        '引言': 'introduction section',
        # 领域术语
        '文献': 'literature publication',
        '基因': 'gene genetics',
        '差异表达': 'differential expression gene expression',
        '表达水平': 'gene expression level',
        '转录组': 'transcriptome rna sequencing',
        '孟德尔': 'mendelian randomization gwas',
        '森林图': 'forest plot meta analysis',
        'OR': 'odds ratio',
        '95%CI': '95 percent confidence interval',
        '可视化': 'data visualization chart',
        '引用': 'citation reference',
        '参考文献': 'bibliographic reference citation',
        # 系统操作
        '代码': 'source code programming',
        'git': 'git version control',
        '推送': 'git push',
        '状态': 'git status',
        '任务': 'task todo',
        '待办': 'task todo',
        '磁盘': 'disk storage',
        '空间': 'storage space',
        '剩余': 'remaining available',
        # 工具
        '邮件': 'email', '邮箱': 'email',
        '桶': 'bucket cloud storage',
        'COS': 'COS Tencent cloud object storage',
        'PubMed': 'PubMed biomedical literature',
        'PubMed上': 'PubMed database search',
        # 生物信息学
        '表达': 'expression',
        '数据集': 'dataset data',
        '数据': 'data',
        '格式': 'format style',
    }

    def _translate_query(self, query: str) -> str:
        """
        翻译桥 v2：纯英文输出，不再拼接中文原文
        让 all-MiniLM-L6-v2 在纯英文空间中匹配 skill 库
        """
        sorted_map = sorted(
            self._ZH2EN_MAP.items(),
            key=lambda x: -len(x[0])
        )
        en_parts = []
        for zh, en in sorted_map:
            if zh in query:
                en_parts.append(en)
        # 去重 + 整理成句子
        seen = set()
        unique = []
        for w in en_parts:
            for tok in w.split():
                if tok not in seen:
                    seen.add(tok)
                    unique.append(tok)
        if not unique:
            return query  # fallback 到原文
        # 纯英文输出：翻译词重复两次（[翻译词] + [翻译词]）增强语义信号
        return ' '.join(unique) + ' . ' + ' '.join(unique)

    def route(self, query: str, top_k: int = 5) -> list[dict]:
        """
        混合检索（BM25 + embedding）→ 返回 top_k 个最匹配的 skill 完整画像

        返回格式:
        [{name, description, dir, path, capabilities, brief_guide, body_start, relevance}, ...]
        """
        if self.embeddings is None:
            self.build_embeddings()

        model = self.model
        if model is None or self.embeddings is None:
            return self.skills[:top_k]

        # 1. BM25 关键词检索（取 top 20）—— 仍然用中文原文
        bm25_indices = self._bm25_search(query, top_k=20)

        # 2. Embedding 语义检索（双通道加权融合）
        # 中英双通道，各算相似度后加权平均（英文路权重更高）
        # 中文路保留字母数字精确匹配（如 PubMed, COS）
        # 英文路增强语义方向（匹配英文 skill 库）
        q_emb_zh = model.encode([query], show_progress_bar=False)[0]
        en_query = self._translate_query(query)
        q_emb_en = model.encode([en_query], show_progress_bar=False)[0]
        q_norm_zh = q_emb_zh / (np.linalg.norm(q_emb_zh) + 1e-10)
        q_norm_en = q_emb_en / (np.linalg.norm(q_emb_en) + 1e-10)
        emb_norm = self.embeddings / (
            np.linalg.norm(self.embeddings, axis=1, keepdims=True) + 1e-10
        )
        # 加权融合：中文0.3 + 英文0.7
        # 英文路权重更高因为 skill 库全英文化
        emb_scores = (
            0.3 * np.dot(emb_norm, q_norm_zh) +
            0.7 * np.dot(emb_norm, q_norm_en)
        )

        # 3. 加权RRF (Weighted Reciprocal Rank Fusion)
        embed_rank = np.argsort(emb_scores)[::-1]  # 高→低
        bm25_rank_map = {idx: r+1 for r, idx in enumerate(bm25_indices)}  # rank从1开始

        k = 60  # RRF 平滑参数
        BM25_WEIGHT = 0.5
        EMBED_WEIGHT = 2.0

        rrf_scores = np.zeros(len(emb_scores))
        for rank, idx in enumerate(embed_rank):
            rrf_scores[idx] += EMBED_WEIGHT / (k + rank + 1)
        for idx, bm_rank in bm25_rank_map.items():
            rrf_scores[idx] += BM25_WEIGHT / (k + bm_rank)

        # 4. 取 top_k 候选
        top_indices = np.argsort(rrf_scores)[::-1][:top_k]
        top_candidates = [self.skills[int(idx)] for idx in top_indices]

        if self._is_noise(query, top_candidates, rrf_scores, top_indices):
            return []

        results = []
        for idx in top_indices:
            skill = dict(self.skills[int(idx)])
            skill["relevance"] = round(float(rrf_scores[int(idx)]), 4)
            skill["relevance_embed"] = round(float(emb_scores[int(idx)]), 4)
            results.append(skill)

        return results

    def _is_noise(self, query: str, candidates: list[dict],
                  rrf_scores: np.ndarray, top_indices: np.ndarray) -> bool:
        """
        噪声检测：三道防线

        0 — 信息查询意图检测：query只是查/找/搜信息，不要求skill做事
            → 强制拒选，直接进L0工具保底
        A — 任务关键词重叠：query中的任务词与top3候选描述子串匹配
        B — 分数阈值：加权RRF top1太低 或 gap太小

        防线0单独触发即拒绝。防线A+B双重触发才拒绝。
        """
        import sys

        q_lower = query.lower()

        # ================================================================
        # 防线0: 信息查询意图检测 + 垃圾输入检测
        # ================================================================
        INFO_SEEK_WORDS = {'查', '找', '搜', '搜索', '检索', '查找', '查询',
                           'search', 'lookup', 'find',
                           '最新', '最近',
                           '有什么', '有哪些', '告诉我',
                           '是什么', '什么是', '怎么', '如何', '为什么',
                           '哪个', '哪种', '哪个好', '区别', '对比',
                           '天气', '时间', '现在几点',
                           '多少钱', '免费', '收费',
                           '教程', '怎么用',
                           '用法', '用量', '用法用量', '剂量',
                           '吃什么', '吃多少', '吃几',
                           '多少钱', '价格',
                           '在哪', '哪里有',
                           '症状', '病因', '治疗',
                           '能治', '能不能',
                           '副作用', '不良反应',
                           '注意事项'}
        STRONG_SKILL_WORDS = {
            # Meta/系统综述
            '系统性综述', '系统综述', 'systematic review',
            'meta分析', 'meta-analysis', 'metaanalysis',
            # 文献计量
            '文献计量', 'bibliometric', 'bibliometrics',
            # 写作/论文
            '写', '撰写', '写作', '做一篇', '论文', 'manuscript',
            'draft', '草稿',
            # 文献
            '文献检索', '文献搜索', '文献', 'literature', 'pubmed',
            # 孟德尔/基因
            'prisma', '孟德尔', 'gwas', 'mendelian',
            'mendelian randomization', 'mr分析',
            # 图表
            '画图', '图表', '可视化', '森林图', 'visualization',
            'forest plot', 'figure',
            # 邮件
            '发邮件', 'send', 'email',
            # 系统操作
            '挂载', 'mount', '重启', 'restart',
            'debug', '审查代码', '重构', 'git',
            '列出来', '列出', '显示', '展示', 'list',
            # 任务管理
            '任务', 'task', '待办', 'todo',
            # 评审
            '评审', '审稿', 'review',
            # 引用
            '引用', '参考文献', 'citation', 'vancouver',
            # 指南
            '指南', 'guideline',
            # 分析
            '差异表达', '转录组', '基因',
        }

        # 垃圾输入检测
        if len(q_lower.strip()) < 2:
            print(f"[router_top5] 输入过短 → 拒选 (len={len(q_lower.strip())})", file=sys.stderr)
            return True
        # 纯符号/数字
        alpha_chars = sum(1 for c in q_lower if c.isalpha())
        if alpha_chars < 2:
            print(f"[router_top5] 无有效文字 → 拒选 (alpha={alpha_chars})", file=sys.stderr)
            return True

        has_info_seek = any(w in q_lower for w in INFO_SEEK_WORDS)
        has_strong_skill = any(w in q_lower for w in STRONG_SKILL_WORDS)

        # 信息查询意图拒选 + 强技能意图豁免
        # 规则: 有 info_seek 且没有 strong_skill → 拒选
        # 例外: query中包含明确的动作指令（列出来、帮我、发给、发邮件等）
        ACTION_PHRASES = {
            '列出来', '列一下', '列个', '列表', '显示', '展示', '列出', '给我',
            '发给', '发邮件', '发一封', '做一篇', '做一',
            '画一个', '画个', '写一篇', '写一',
            '查一下', '查一',
            '把',
        }
        has_action = any(ap in q_lower for ap in ACTION_PHRASES)

        if has_info_seek and not has_strong_skill and not has_action:
            # 移除年份数字（2024/2025/2026）导致的误拦
            year_pattern = re.search(r'20\d{2}', q_lower)
            # 如果仅因为年份触发 + 有任务词 → 放行
            if year_pattern:
                year_only_info_seek = all(w not in q_lower or w in ('最新','最近') or len(w) < 2 
                                          for w in INFO_SEEK_WORDS)
                if year_only_info_seek:
                    pass  # 放行
                else:
                    print(f"[router_top5] 信息查询意图 → 拒选L1，直接L0工具保底 "
                          f"(info_seek={has_info_seek} strong_skill={has_strong_skill})", file=sys.stderr)
                    return True
            else:
                print(f"[router_top5] 信息查询意图 → 拒选L1，直接L0工具保底 "
                      f"(info_seek={has_info_seek} strong_skill={has_strong_skill})", file=sys.stderr)
                return True

        # ================================================================
        # 防线A+B: 任务词重叠 + 分数阈值
        # ================================================================

        top1_score = rrf_scores[top_indices[0]]
        top_k = len(top_indices)
        top_k_score = rrf_scores[top_indices[-1]] if top_k >= 3 else 0
        gap_ratio = (top1_score - top_k_score) / (top1_score + 1e-10)

        # 翻译后 RRF 分数偏低（0.01-0.06），阈值相应调低
        MIN_ABSOLUTE_RRF = 0.008
        MIN_GAP_RATIO = 0.06

        # 防线B: 分数阈值
        score_noise = (top1_score < MIN_ABSOLUTE_RRF
                       or gap_ratio < MIN_GAP_RATIO)

        # 防线A: 任务关键词重叠（中英文双轨匹配）
        TASK_WORDS = {
            # Meta/系统综述
            '系统性综述', '系统综述', '系统评价', 'systematic review',
            'meta分析', 'meta-analysis', 'metaanalysis',
            # 文献计量
            '文献计量', 'bibliometric', 'bibliometrics',
            # 文献检索
            '文献检索', '文献搜索', '搜索文献', '检索文献',
            'literature search', 'literature review',
            'prisma', 'meta',
            # 指南
            '指南', 'guideline', 'guidelines',
            # 综述
            '综述', 'review',
            # 论文
            '论文', 'article', 'paper', 'manuscript',
            # 写作
            '写', '撰写', '写作', 'draft',
            # 分析（中英文双轨, 因为 capabilities 已全英文化）
            '分析', 'analysis', '统计', 'statistics',
            # 图表
            '图表', 'chart', '可视化', 'visualization',
            '画图', '森林图', 'forest plot',
            'figure', 'plot',
            # 邮件
            '邮件', 'email', 'mail',
            # 下载
            '下载', 'download', 'pdf',
            # 基因/孟德尔
            'gwas', '基因', 'gene', 'genetics',
            '孟德尔', 'mr分析', 'mendelian',
            # R/代码（'r '太短易漏过，已移除，用 r语言/rstats/rstudio 代替）
            'r语言', 'rstats', 'rstudio',
            # 引用
            '引用', '参考文献', 'reference', 'citation',
            'apa', 'vancouver', 'bibliography',
            # 评审
            '评审', '审稿', 'peer review',
            # 系统操作
            '挂载', 'mount', '重启', 'restart',
            # 调试
            'debug', '调试', 'bug',
            # 任务管理
            '任务', 'task', 'todo', '待办', '优先级', 'priority',
            '分配', 'schedule',
            # 列表
            '列出来', '列出', '列表', 'list',
            # Git
            'git', 'push', 'commit',
            # 审查
            '审查', 'review code',
        }
        # TASK_WORDS 已屏蔽过度宽泛的条目（如 'r ' 太短易漏过）
        _TASK_MIN_LEN = 2

        q_lower = query.lower()
        # 找出query中命中的任务词（过滤过短词，避免 'r ' 等短路词漏过）
        hit_tasks = {tw for tw in TASK_WORDS if tw in q_lower and len(tw) >= _TASK_MIN_LEN}

        # 防线A: query中有可识别的任务词 → 直接确认有任务意图
        # skill库已全英文化，中文任务词在英文capabilities中找不到子串匹配
        # 因此 hit_tasks 非空本身就说明用户有明确的任务意图
        # Phase 3 修正：hit_tasks非空→视为有任务重叠，不再检查skill_text
        has_overlap = True if hit_tasks else False

        # 必须双重触发：零任务重叠 AND 分数也低
        if (not has_overlap) and score_noise:
            print(f"[router_top5] 噪声拒绝: 任务词无重叠(hit={hit_tasks}) "
                  f"+ top1_rrf={top1_score:.4f} gap={gap_ratio:.3f}", file=sys.stderr)
            return True

        return False

    def route_formatted(self, query: str, top_k: int = 5) -> str:
        """返回格式化的文本块，直接给 LLM"""
        results = self.route(query, top_k=top_k)
        
        if not results:
            return "--- 无可选技能（噪声拒绝：所有候选与查询无关）---\n\n请回复 '拒选'。"
        
        from schemas import SkillProfile

        blocks = []
        for i, r in enumerate(results):
            profile = SkillProfile(
                name=r["name"],
                description=r.get("description", ""),
                dir=r.get("dir", ""),
                path=r.get("path", ""),
                capabilities=r.get("capabilities", []),
                brief_guide=r.get("brief_guide", ""),
                body_start=r.get("body_start", ""),
                relevance=r.get("relevance", 0.0),
            )
            blocks.append(profile.to_prompt_block(i + 1))

        header = f"--- 可选技能（语义匹配度 top{top_k}）---\n"
        footer = "\n--- 请从以上技能中选择最匹配的一个或多个。如果不匹配任何技能，请回复 '拒选'。---"
        return header + "\n\n".join(blocks) + footer


# ============================================================================
# 独立调用入口
# ============================================================================
if __name__ == "__main__":
    import sys

    router = Top5Router()
    if not router.embeddings:
        print("首次运行：构建 embeddings...")
        router.build_embeddings()

    if len(sys.argv) > 1:
        query = " ".join(sys.argv[1:])
    else:
        query = "做一篇脓毒症脑病的文献计量分析"

    print(f"Query: {query}\n")
    print(router.route_formatted(query, top_k=5))
