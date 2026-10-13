#!/usr/bin/env python3
"""Build an MSA from LaRA2's T-Coffee library using the official workflow."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def executable_from_env(name: str, executable: str) -> Path:
    raw_directory = os.environ.get(name)
    if not raw_directory:
        raise ValueError(f"{name} must point to the {executable} bin directory")
    path = Path(raw_directory) / executable
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"{executable} is not executable: {path}")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    try:
        lara = executable_from_env("LARA2_BIN_DIR", "lara")
        t_coffee = executable_from_env("TCOFFEE_BIN_DIR", "t_coffee")
    except (ValueError, FileNotFoundError) as error:
        print(f"LaRA2/T-Coffee setup error: {error}", file=sys.stderr)
        return 2

    input_path = args.input.resolve()
    output_path = args.output.resolve()
    work_dir = output_path.parent
    work_dir.mkdir(parents=True, exist_ok=True)
    tcoffee_home = work_dir / "tcoffee-home"
    tcoffee_tmp = work_dir / "tcoffee-tmp"
    tcoffee_home.mkdir(exist_ok=True)
    tcoffee_tmp.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment["HOME_4_TCOFFEE"] = str(tcoffee_home)
    environment["TMP_4_TCOFFEE"] = str(tcoffee_tmp)

    library_path = work_dir / "lara2-library.tc_lib"
    lara_command = [
        str(lara), "--version-check", "0", "--threads", "1",
        "--infile", str(input_path), "--outformat", "lib",
        "--write", str(library_path),
    ]
    lara_result = subprocess.run(lara_command, cwd=work_dir,
                                 env=environment, check=False)
    if lara_result.returncode != 0:
        return lara_result.returncode
    if not library_path.is_file() or library_path.stat().st_size == 0:
        print("LaRA2 did not produce a T-Coffee library", file=sys.stderr)
        return 1

    tcoffee_command = [
        str(t_coffee), "-lib", str(library_path), "-seq", str(input_path),
        "-outfile", str(output_path), "-output", "fasta_aln",
        "-thread", "1", "-n_core", "1",
    ]
    result = subprocess.run(tcoffee_command, cwd=work_dir,
                            env=environment, check=False)
    if result.returncode != 0:
        return result.returncode
    if not output_path.is_file() or output_path.stat().st_size == 0:
        print("T-Coffee did not produce the requested MSA", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
