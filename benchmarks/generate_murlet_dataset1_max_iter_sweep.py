#!/usr/bin/env python3
"""Generate a max-iter sweep from the selected Murlet dataset1 conditions."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/configs/murlet-dataset1-condition-sweep-arm64-v272.json"
OUTPUT = ROOT / "benchmarks/configs/murlet-dataset1-max-iter-sweep-arm64-v272.json"

SELECTED = [
    "vp_alifold_w2",
    "cc_alifold_fth03",
    "cl_linear_ath02",
    "lpv_probcons_fpct05",
    "lpc_contralign_ath02",
    "lpv_linear_fpct05",
]
MAX_ITERS = range(100, 1001, 100)


def replace_max_iter(args: list[str], max_iter: int) -> list[str]:
    replaced = False
    result: list[str] = []
    for arg in args:
        if arg.startswith("--max-iter="):
            result.append(f"--max-iter={max_iter}")
            replaced = True
        else:
            result.append(arg)
    if not replaced:
        result.append(f"--max-iter={max_iter}")
    return result


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    by_id = {condition["id"]: condition for condition in source["conditions"]}
    missing = [condition_id for condition_id in SELECTED if condition_id not in by_id]
    if missing:
        raise SystemExit(f"selected conditions are missing: {missing}")

    conditions: list[dict[str, object]] = []
    for condition_id in SELECTED:
        base = by_id[condition_id]
        for max_iter in MAX_ITERS:
            condition = copy.deepcopy(base)
            condition["id"] = f"{condition_id}_iter{max_iter}"
            condition["phase"] = "max_iter_sweep"
            condition["base_condition"] = condition_id
            condition["max_iter"] = max_iter
            condition["args"] = replace_max_iter(condition["args"], max_iter)
            conditions.append(condition)

    output = copy.deepcopy(source)
    output["output_dir"] = "../results-murlet-dataset1-max-iter-sweep-arm64-v272"
    output["timeout_seconds"] = 3600
    output["conditions"] = conditions
    output["sweep"] = {
        "type": "max_iter",
        "selected_conditions": SELECTED,
        "values": list(MAX_ITERS),
        "selection_order": ["mcc", "sci", "sps"],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"conditions={len(conditions)} datasets=85 runs={len(conditions) * 85}")
    print(OUTPUT)


if __name__ == "__main__":
    main()
