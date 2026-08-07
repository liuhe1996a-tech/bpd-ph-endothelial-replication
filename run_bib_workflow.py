"""Run the frozen base workflow and the virtual-cell benchmark extension."""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Execute; default prints both plans.")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--ai-config", default="AI_upgrade/workflow_config.example.json")
    parser.add_argument("--skip-base", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if not args.skip_base:
        cmd = [sys.executable, str(root / "run_all.py"), "--rscript", args.rscript]
        if args.execute:
            cmd.append("--execute")
        subprocess.run(cmd, cwd=root, check=True)
    cmd = [sys.executable, str(root / "AI_upgrade/scripts/00_run_ai_upgrade.py"), "--config", args.ai_config]
    if not args.execute:
        cmd.append("--dry-run")
    subprocess.run(cmd, cwd=root, check=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
