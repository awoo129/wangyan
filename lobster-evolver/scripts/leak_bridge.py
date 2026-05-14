#!/usr/bin/env python3
"""
leak_detector → self-evolution 桥接
读取最后 2 次 leak_alerts.jsonl，连续 CRITICAL 则写入 ATOMIC lesson
"""

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
LOGS_DIR = SCRIPT_DIR.parent / "logs"
LEAK_LOG = LOGS_DIR / "leak_alerts.jsonl"


def read_last_n(n: int = 2) -> list[dict]:
    """读取最后 n 条 leak 日志"""
    if not LEAK_LOG.exists():
        return []
    lines = LEAK_LOG.read_text(encoding="utf-8").strip().split("\n")
    last_lines = [l for l in lines if l.strip()][-n:]
    records = []
    for line in last_lines:
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return records


def write_atomic_lesson(record: dict):
    """将 leak 告警写入 ATOMIC lesson"""
    from evolution_event import ATOMICLessonStore, EvolutionEvent

    timestamp = record.get("timestamp", "unknown")
    checks = record.get("checks", {})

    evt = EvolutionEvent(
        category="process_leak",
        summary=f"连续进程泄漏告警: {checks.get('gateway_processes', {}).get('status', '?')}/{checks.get('process_surge', {}).get('status', '?')}",
        description=json.dumps({
            "gateway": checks.get("gateway_processes", {}),
            "session": checks.get("session_leaks", {}),
            "orphan_python": checks.get("orphan_python", {}),
            "process_surge": checks.get("process_surge", {}),
        }, ensure_ascii=False),
        severity="critical",
        source="leak_detector",
    )

    store = ATOMICLessonStore()
    store.emit(evt)
    print(f"[leak_bridge] ATOMIC lesson written: process_leak@{timestamp}", file=sys.stderr)


def main():
    records = read_last_n(2)

    if len(records) < 2:
        print("[leak_bridge] 不足2条记录，跳过", file=sys.stderr)
        return

    # 检查最后2次是否都是 alert
    last_alerts = [r for r in records if r.get("alert")]
    criticals = [r for r in records if r.get("alert_level") == "critical"]

    if len(criticals) >= 2:
        print("[leak_bridge] 连续2次CRITICAL → 写入ATOMIC lesson", file=sys.stderr)
        write_atomic_lesson(criticals[-1])
    elif len(last_alerts) >= 2:
        print("[leak_bridge] 连续2次告警（含WARNING），写入ATOMIC lesson", file=sys.stderr)
        write_atomic_lesson(last_alerts[-1])
    else:
        print("[leak_bridge] 无连续告警，跳过", file=sys.stderr)


if __name__ == "__main__":
    main()
