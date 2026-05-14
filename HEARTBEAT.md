# HEARTBEAT.md — Periodic Checks

## Memory Hygiene (memory-tdai)
**Trigger**: Weekly (every 168 heartbeat cycles ≈ weekly)
**Action**: Check memory-tdai data size, verify FTS5 indexes
**Note**: memory-tdai auto-cleans disabled (retentionDays not configured)

## Gateway Token Monitor (gateway_token_monitor)
**Trigger**: Every 3 heartbeat cycles (≈ every 90 min active usage)
**Action**:
1. Run `python3 lobster-evolver/scripts/gateway_token_monitor.py protect`
2. If YELLOW (>150k tokens) → auto-compress session (keep 200 msgs)
3. If RED (>200k tokens) → auto-compress + trigger `emergency_rescue.py --reason TOKEN_OVERFLOW_AUTO`
4. Alert if any YELLOW/RED found
**Purpose**: Proactive token management — don't wait for 200k overflow. Compress early.

## Skill Internalizer (skill-internalizer)
**Trigger**: Every 10 heartbeat cycles (≈ every 5 hours active usage)
**Action**:
1. Run `node skills/skill-internalizer/skill-internalizer.js analyze` — check for CANDIDATEs
2. If candidates exist, run `node skills/skill-internalizer/skill-internalizer.js notify` — push to evolution systems
3. capability-evolver reads `memory/skill-internalizer-candidates.md` on next run
4. self-improving reads `lessons.md` for error-driven evolution rules
**Purpose**: Track skill usage effectiveness and promote high-scoring skills from SKILL.md → MEMORY.md (muscle memory), inspired by Skill0's "Skills at training, zero at inference" paradigm.

## Tool Audit Anomaly Check (tool_audit)
**Trigger**: Every 6 heartbeat cycles (≈ every 3 hours active usage)
**Action**:
1. Run `python3 lobster-evolver/scripts/tool_audit.py summarize 24`
2. If anomalies detected → print alert to session
3. If anomaly count > 5 → write to `lobster-evolver/logs/anomaly_alert.md`
4. **Immediately follow with**: `python3 lobster-evolver/scripts/process_anomalies.py`
   - Reads `error_log.json` (new unprocessed errors)
   - Writes ATOMIC lessons → `lobster-evolver/lessons/ATOMIC/`
   - Infers and writes FUNCTIONAL lessons → `lobster-evolver/lessons/FUNCTIONAL/`
   - Marks error IDs as processed in `lobster-evolver/logs/.processed_anomalies`
5. **Bridge to self-improving**: `python3 lobster-evolver/scripts/bridge_to_self_improving.py`
   - Reads ATOMIC/FUNCTIONAL lessons → writes to `~/.self-improving/memory.md`
   - Reads error_log.json → writes to `~/.self-improving/corrections.md`
   - Two systems now share knowledge automatically
**Purpose**: Closes the loop: detect anomaly → extract rules → persist to skill library → bridge to self-improving-agent. Error-driven evolution now runs automatically, no manual "submit" needed.

## Skill Metrics Check (skill_metrics) ← NEW from OpenSpace GDPVal insight
**Trigger**: Every 6 heartbeat cycles (≈ every 3 hours active usage, same as tool_audit)
**Action**:
1. Run `python3 lobster-evolver/scripts/skill_metrics.py summarize 168` — weekly skill health
2. Run `python3 lobster-evolver/scripts/skill_metrics.py degradation` — check all skills for success rate regression
3. Run `python3 lobster-evolver/scripts/skill_metrics.py capture_check` — find CAPTURED lesson opportunities
4. If degradation detected → write to `lobster-evolver/logs/skill_degradation_alert.md`
5. If CAPTURED opportunities found → write to `lobster-evolver/lessons/CAPTURED/`
**Purpose**: OpenSpace GDPVal发现 165 skills 中，44个是File I/O fallback、29个是Execution Recovery——龙虾以前没有技能级质量监控，这补上了"技能退化检测"的最后一块板。比tool_audit高一个层级：tool_audit看工具调用，skill_metrics看技能效果。

## 🦞 Resource Leak Detection (leak_detector)
**Trigger**: Every 6 heartbeat cycles (≈ every 3 hours active usage, same as tool_audit)
**Action**:
1. Run `python3 lobster-evolver/scripts/leak_detector.py`
2. If CRITICAL → print alert to session
3. If WARNING → write to `lobster-evolver/logs/leak_alerts.jsonl`
4. **Bridge to self-evolution**: `python3 lobster-evolver/scripts/leak_bridge.py`
   - 读取最后2次 leak_alerts.jsonl
   - 连续 CRITICAL → 写入 ATOMIC lesson 到 evolution 管道
**Purpose**: 进程泄漏（如CodeBrain MCP泄漏模式）会导致skill调用链路全崩。集成到心跳能早发现早处理，不走skill路由层面。

## No other memory skills auto-run
- memory-strategy: manual-only (triggered by "记下来" etc.)
- memory-manager: manual-only (detect.sh, organize.sh)
- memory-setup: manual configuration helper
- memory-hygiene: manual LanceDB cleanup
- mempalace: disabled (cron removed, 重复 memory-tdai)

## Reverse Prompting (proactive)
**Trigger**: Every 6 heartbeat cycles (≈ every 3 hours active usage)
**Action**:
1. Read `USER.md` + recent `memory/YYYY-MM-DD.md`
2. Ask: "有什么我可以主动帮你的？基于我目前对你的了解"
3. 发现机会 → 主动提议（但不直接执行外部操作）
4. 将提议写入 `memory/proactive-proposals.md`

**触发词**: 当人类说"你有什么想法吗"或"你能帮我做什么"时，优先从 proactive-proposals.md 读取

## Working Buffer Check (danger zone)
**Trigger**: Every heartbeat
**Action**:
1. Run `python3 lobster-evolver/scripts/working_buffer.py --check`
2. If `in_danger: true` → 之后每条消息自动追加到 buffer:
   - Human: `python3 working_buffer.py --append --role human --message "..."`
   - Agent: `python3 working_buffer.py --append --role agent --summary "..."`
3. If `in_danger: false` → 正常模式，不追加

## Disk Alert Check (disk_alert)
**Trigger**: Every heartbeat (every ~30 min, real-time check)
**Action**:
1. Check `logs/.disk_alert_pending` file
2. If exists → read content and send QQ alert via `message` tool
3. After alert sent → delete `logs/.disk_alert_pending`
**Why**: cron runs in isolated session, cannot send QQ directly — heartbeat bridges the gap.
