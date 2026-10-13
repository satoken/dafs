#!/usr/bin/env python3
"""Re-run the Murlet dataset1 condition grid with LinearAlifold enabled."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/configs/murlet-dataset1-condition-sweep-arm64-v272.json"
OUTPUT = ROOT / "benchmarks/configs/murlet-dataset1-linearalifold-iter100-beam100.json"
BEAMS = ("--align-beam", "--align-dd-beam", "--linfold-beam",
         "--fold-dd-beam", "--fold-final-beam")


def normalize_args(args: list[str], linear: bool) -> list[str]:
    result = [arg for arg in args if not arg.startswith("--max-iter=")]
    result.append("--max-iter=100")
    if linear:
        for option in BEAMS:
            result = [arg for arg in result if not arg.startswith(option + "=")]
            result.append(option + "=100")
    return result


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    conditions = []
    for old in source["conditions"]:
        condition = dict(old)
        linear = bool(condition["linear"])
        condition["args"] = normalize_args(list(condition["args"]), linear)
        conditions.append(condition)
        if linear:
            alifold = dict(condition)
            alifold["id"] = condition["id"] + "_alifold"
            alifold["label"] = condition["label"] + " (LinearAlifold)"
            alifold["args"] = ["--alifold" if arg == "--no-alifold" else arg
                               for arg in condition["args"]]
            conditions.append(alifold)
    source["conditions"] = conditions
    source["output_dir"] = "../results-murlet-dataset1-linearalifold-iter100-beam100"
    OUTPUT.write_text(json.dumps(source, indent=2) + "\n", encoding="utf-8")
    print(f"{len(conditions)} conditions x 85 datasets = {len(conditions) * 85} runs")
    print(OUTPUT)


if __name__ == "__main__":
    main()
