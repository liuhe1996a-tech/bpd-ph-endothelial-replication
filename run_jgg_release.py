"""Stable public entry point for the frozen JGG analysis release."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> int:
    command = [sys.executable, str(ROOT / "scripts" / "run_all_r10.py"), *sys.argv[1:]]
    return subprocess.call(command, cwd=ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
