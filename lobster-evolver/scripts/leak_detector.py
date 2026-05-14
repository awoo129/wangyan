#!/usr/bin/env python3
"""
资源泄漏检测器 — 参考 the-openclaw-optimizer 的 leak detection 思路
检测 OpenClaw 生态中的进程/会话泄漏

监控维度：
1. OpenClaw gateway 子进程数（MCP、CodeBrain 等）
2. session 文件残留（僵尸会话）
3. subprocess run 超时/僵尸
4. 进程数变化趋势

输出：leak 报告 → stdout + logs/leak_alerts.jsonl
"""

import json
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

LOGS_DIR = Path(__file__).resolve().parent.parent / "logs"


# =============================================================================
# 配置
# =============================================================================

THRESHOLDS = {
    "gateway_processes_max": 15,    # gateway 主进程 + 子进程上限
    "mcp_processes_max": 5,         # MCP 子进程上限
    "zombie_sessions_max": 3,       # 僵尸 session 文件上限
    "orphan_python_max": 10,        # 孤立 Python 进程上限
    "subprocess_timeout": 300,      # subprocess.run 超时秒数
    "process_surge_delta": 10,      # 进程数突增阈值（与上次对比）
}


# =============================================================================
# 检测函数
# =============================================================================

def check_gateway_processes() -> dict:
    """检测 gateway 进程树"""
    result = {
        "gateway_count": 0,
        "mcp_count": 0,
        "zombie_count": 0,
        "status": "ok",
    }

    try:
        # 查找 openclaw-gateway 主进程
        ps_out = subprocess.run(
            ["pgrep", "-a", "openclaw"],
            capture_output=True, text=True, timeout=10
        )
        lines = ps_out.stdout.strip().split("\n") if ps_out.stdout.strip() else []
        
        gateway_pids = []
        for line in lines:
            if "gateway" in line.lower():
                parts = line.split(None, 1)
                if parts:
                    gateway_pids.append(int(parts[0]))
        result["gateway_count"] = len(gateway_pids)
        
        # 检查 MCP 进程（按路径模式匹配，避免误匹配 scp、microphone 等）
        try:
            # 从 /proc 读取 cmdline 匹配 mcp 关键词，更精确
            ps_out = subprocess.run(
                ["ps", "-eo", "pid,args", "--no-headers"],
                capture_output=True, text=True, timeout=5
            )
            mcp_pids = []
            for line in ps_out.stdout.strip().split("\n"):
                if not line.strip():
                    continue
                parts = line.strip().split(None, 1)
                if len(parts) < 2:
                    continue
                pid_str, args = parts
                # 匹配 mcp 作为独立单词或路径段
                if re.search(r'(?:\s|/|\b)mcp(?:[-_]|\b)', args, re.IGNORECASE):
                    try:
                        mcp_pids.append(int(pid_str))
                    except ValueError:
                        pass
            result["mcp_count"] = len(mcp_pids)
        except (subprocess.TimeoutExpired, ValueError):
            pass
        
        # 检查僵尸进程（准确方法：看 STAT 列）
        try:
            zombie = subprocess.run(
                ["ps", "-eo", "stat,pid,comm", "--no-headers"],
                capture_output=True, text=True, timeout=5
            )
            if zombie.stdout.strip():
                zombie_lines = [l for l in zombie.stdout.strip().split("\n") 
                              if l.strip() and l.strip()[0] == 'Z']
                result["zombie_count"] = len(zombie_lines)
        except subprocess.TimeoutExpired:
            pass
        
        # 状态判定
        total = result["gateway_count"] + result["mcp_count"]
        if total > THRESHOLDS["gateway_processes_max"]:
            result["status"] = "critical"
        elif result["mcp_count"] > THRESHOLDS["mcp_processes_max"]:
            result["status"] = "warning"
        elif result["zombie_count"] > 0:
            result["status"] = "warning"
        else:
            result["status"] = "ok"
            
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
    
    return result


def check_session_leaks() -> dict:
    """检测 session 文件残留（僵尸会话）"""
    result = {
        "session_count": 0,
        "zombie_sessions": [],
        "status": "ok",
    }
    
    # 查找 OpenClaw session 目录
    session_dirs = [
        Path.home() / ".openclaw" / "session",
        Path.home() / ".local" / "share" / "openclaw" / "sessions",
        Path("/tmp") / "openclaw-sessions",
    ]
    
    for session_dir in session_dirs:
        if session_dir.exists():
            try:
                sessions = list(session_dir.iterdir())
                result["session_count"] += len(sessions)
                
                # 检测超过24小时无活动的 session
                # 对目录：取内部文件的最新 mtime（目录 mtime 不随内部文件更新）
                now = time.time()
                for s in sessions:
                    if s.is_file():
                        mtime = s.stat().st_mtime
                        stale_hours = (now - mtime) / 3600
                        if stale_hours > 24:
                            result["zombie_sessions"].append(f"{s.name}({stale_hours:.0f}h)")
                    elif s.is_dir():
                        # 取目录内最新文件的 mtime（最多遍历 5000 文件防 hang）
                        latest_mtime = s.stat().st_mtime
                        max_files = 5000
                        for i, f in enumerate(s.rglob("*")):
                            if i >= max_files:
                                break
                            try:
                                ft = f.stat().st_mtime
                                if ft > latest_mtime:
                                    latest_mtime = ft
                            except OSError:
                                pass
                        stale_hours = (now - latest_mtime) / 3600
                        if stale_hours > 24:
                            result["zombie_sessions"].append(f"{s.name}/({stale_hours:.0f}h)")
            except (PermissionError, OSError):
                pass
    
    if len(result["zombie_sessions"]) > THRESHOLDS["zombie_sessions_max"]:
        result["status"] = "warning"
    
    return result


def check_orphan_python() -> dict:
    """检测孤立 Python 进程（不是 gateway 子进程的 python）"""
    result = {
        "orphan_count": 0,
        "orphan_pids": [],
        "status": "ok",
    }
    
    # 获取自身 PID，排除自己
    import os
    self_pid = os.getpid()
    
    try:
        # 获取所有 python3 进程
        ps_out = subprocess.run(
            ["ps", "-eo", "pid,ppid,comm", "--no-headers"],
            capture_output=True, text=True, timeout=10
        )
        
        # 获取 gateway 进程树的所有 pid
        gateway_pid = None
        try:
            gw = subprocess.run(
                ["pgrep", "-f", "openclaw-gateway"],
                capture_output=True, text=True, timeout=5
            )
            if gw.stdout.strip():
                gateway_pid = int(gw.stdout.strip().split("\n")[0])
        except (ValueError, IndexError, subprocess.TimeoutExpired):
            pass
        
        gateway_tree = set()
        if gateway_pid:
            try:
                tree = subprocess.run(
                    ["pstree", "-p", str(gateway_pid)],
                    capture_output=True, text=True, timeout=5
                )
                gateway_tree = set(int(p) for p in re.findall(r'\((\d+)\)', tree.stdout))
            except (subprocess.TimeoutExpired, FileNotFoundError):
                pass
        
        # 遍历所有 python3 进程
        for line in ps_out.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                pid, ppid, comm = int(parts[0]), int(parts[1]), parts[2]
                if "python" in comm.lower() and pid not in gateway_tree and pid != self_pid:
                    result["orphan_pids"].append(pid)
            except (ValueError, IndexError):
                pass
        
        result["orphan_count"] = len(result["orphan_pids"])
        
        if result["orphan_count"] > THRESHOLDS["orphan_python_max"]:
            result["status"] = "warning"
            
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
    
    return result


def check_process_surge(state_file: Path) -> dict:
    """
    进程数突增检测：与上次记录对比
    用于检测短时间内进程数暴增（CodeBrain 泄漏模式）
    """
    result = {
        "current_count": 0,
        "last_count": 0,
        "delta": 0,
        "status": "ok",
    }
    
    try:
        # 获取当前总进程数
        ps_out = subprocess.run(
            ["ps", "-e", "--no-headers"],
            capture_output=True, text=True, timeout=5
        )
        current_count = len([l for l in ps_out.stdout.strip().split("\n") if l])
        result["current_count"] = current_count
        
        # 读取上次记录。state 文件不存在 = 首次运行，只设基线不告警
        is_first_run = not state_file.exists()
        last_count = 0
        if not is_first_run:
            try:
                state = json.loads(state_file.read_text())
                last_count = state.get("process_count", 0)
            except (json.JSONDecodeError, OSError):
                pass
        
        result["last_count"] = last_count
        result["delta"] = current_count - last_count if not is_first_run else 0
        
        # 写回当前值
        state_file.parent.mkdir(parents=True, exist_ok=True)
        state_file.write_text(json.dumps({
            "process_count": current_count,
            "timestamp": datetime.now().isoformat(),
        }, ensure_ascii=False))
        
        if not is_first_run and result["delta"] > THRESHOLDS["process_surge_delta"]:
            result["status"] = "critical"
            
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
    
    return result


# =============================================================================
# 统一检测入口
# =============================================================================

def detect_all(quiet: bool = False) -> dict:
    """执行全部泄漏检测，返回完整报告"""
    report = {
        "timestamp": datetime.now().isoformat(),
        "checks": {
            "gateway_processes": check_gateway_processes(),
            "session_leaks": check_session_leaks(),
            "orphan_python": check_orphan_python(),
            "process_surge": check_process_surge(
                LOGS_DIR / ".leak_detector_state.json"
            ),
        },
        "alert": False,
    }
    
    # 全局 alert 判定
    criticals = [c for c in report["checks"].values() if c.get("status") == "critical"]
    warnings = [c for c in report["checks"].values() if c.get("status") == "warning"]
    
    if criticals:
        report["alert"] = True
        report["alert_level"] = "critical"
    elif warnings:
        report["alert"] = True
        report["alert_level"] = "warning"
    
    # 记录日志
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / "leak_alerts.jsonl"
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(report, ensure_ascii=False) + "\n")
    
    if not quiet:
        _print_report(report)
    
    return report


def _print_report(report: dict):
    """打印可读报告"""
    checks = report["checks"]
    
    if report.get("alert_level") == "critical":
        print("🔴 泄漏检测: CRITICAL!", file=sys.stderr)
    elif report.get("alert_level") == "warning":
        print("🟡 泄漏检测: WARNING!", file=sys.stderr)
    else:
        print("🟢 泄漏检测: OK", file=sys.stderr)
    
    print(f"  gateway 进程: {checks['gateway_processes']['gateway_count']} "
          f"(MCP: {checks['gateway_processes']['mcp_count']}, "
          f"僵尸: {checks['gateway_processes']['zombie_count']})",
          file=sys.stderr)
    
    print(f"  session 文件: {checks['session_leaks']['session_count']} "
          f"(僵尸: {len(checks['session_leaks']['zombie_sessions'])})",
          file=sys.stderr)
    
    print(f"  孤立 Python: {checks['orphan_python']['orphan_count']}",
          file=sys.stderr)
    
    print(f"  进程突增: {checks['process_surge']['delta']:+d} "
          f"(上次={checks['process_surge']['last_count']}, "
          f"当前={checks['process_surge']['current_count']})",
          file=sys.stderr)
    
    if report["alert"]:
        print(f"  ⚠️ 告警级别: {report.get('alert_level', 'unknown')}", file=sys.stderr)


# =============================================================================
# CLI 入口
# =============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="资源泄漏检测器")
    parser.add_argument("--quiet", "-q", action="store_true", help="静默模式")
    parser.add_argument("--thresholds", "-t", help="阈值配置JSON")
    args = parser.parse_args()
    
    if args.thresholds:
        try:
            custom = json.loads(args.thresholds)
            # 校验：所有阈值必须为正整数
            known_keys = THRESHOLDS.keys()
            for k, v in custom.items():
                if k not in known_keys:
                    print(f"未知阈值: {k}", file=sys.stderr)
                    sys.exit(1)
                if not isinstance(v, (int, float)) or v <= 0:
                    print(f"阈值 {k} 必须是正数 (got {v})", file=sys.stderr)
                    sys.exit(1)
            THRESHOLDS.update(custom)
        except json.JSONDecodeError:
            print("阈值配置解析失败", file=sys.stderr)
            sys.exit(1)
    
    report = detect_all(quiet=args.quiet)
    
    if report["alert"]:
        sys.exit(1 if report.get("alert_level") == "critical" else 2)
