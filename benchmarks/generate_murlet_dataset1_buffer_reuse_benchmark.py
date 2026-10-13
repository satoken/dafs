#!/usr/bin/env python3
"""Generate the before/after benchmark for reusable linear workspaces."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/configs/murlet-dataset1-max-iter-sweep-arm64-v272.json"
OUTPUT = ROOT / "benchmarks/configs/murlet-dataset1-buffer-reuse-arm64-v272.json"

LINEAR_CONDITIONS = {
    "cl_linear_ath02",
    "lpc_contralign_ath02",
    "lpv_linear_fpct05",
    "lpv_probcons_fpct05",
}
MAX_ITERS = {100, 1000}


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    conditions = [
        copy.deepcopy(condition)
        for condition in source["conditions"]
        if condition.get("base_condition") in LINEAR_CONDITIONS
        and condition.get("max_iter") in MAX_ITERS
    ]
    expected = len(LINEAR_CONDITIONS) * len(MAX_ITERS)
    if len(conditions) != expected:
        raise SystemExit(f"expected {expected} conditions, got {len(conditions)}")

    output = copy.deepcopy(source)
    output["output_dir"] = "../results-murlet-dataset1-buffer-reuse-arm64-v272"
    output["conditions"] = conditions
    output["comparison"] = {
        "baseline_results": "../results-murlet-dataset1-max-iter-sweep-arm64-v272",
        "purpose": "verify reusable BeamAlign and LinearNussinov workspaces",
        "selected_conditions": sorted(LINEAR_CONDITIONS),
        "max_iter_values": sorted(MAX_ITERS),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"conditions={len(conditions)} datasets=85 runs={len(conditions) * 85}")
    print(OUTPUT)


if __name__ == "__main__":
    main()
