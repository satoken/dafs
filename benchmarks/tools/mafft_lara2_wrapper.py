#!/usr/bin/env python3
"""Adapt MAFFT's LaRA helper invocation to the LaRA 2 command line."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    raw_bin_dir = os.environ.get("LARA2_BIN_DIR")
    if not raw_bin_dir:
        print("LARA2_BIN_DIR must point to the LaRA2 bin directory", file=sys.stderr)
        return 2
    binary = Path(raw_bin_dir) / "lara"
    if not binary.is_file() or not os.access(binary, os.X_OK):
        print(f"LaRA2 executable is not available: {binary}", file=sys.stderr)
        return 2

    command = [str(binary), "--version-check", "0", "--threads", "1", *sys.argv[1:]]
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
