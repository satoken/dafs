#!/usr/bin/env python3
"""Generate the max-iter sweep for the current RNAalifold profile module."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/configs/murlet-dataset1-profile-energy-comparison-arm64-v272.json"
OUTPUT = ROOT / "benchmarks/configs/murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272.json"
SELECTED = [
    "lpc_linear_rnaalifold",
    "lpc_contralign_rnaalifold",
    "lpv_linear_rnaalifold",
    "lpv_probcons_rnaalifold",
]
MAX_ITERS = range(100, 1001, 100)


def replace_max_iter(args: list[str], value: int) -> list[str]:
    result = []
    replaced = False
    for arg in args:
        if arg.startswith("--max-iter="):
            result.append(f"--max-iter={value}")
            replaced = True
        else:
            result.append(arg)
    if not replaced:
        result.append(f"--max-iter={value}")
    return result


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    by_id = {condition["id"]: condition for condition in source["conditions"]}
    missing = [condition_id for condition_id in SELECTED if condition_id not in by_id]
    if missing:
        raise SystemExit(f"selected conditions are missing: {missing}")

    conditions = []
    for condition_id in SELECTED:
        base = by_id[condition_id]
        base_name = condition_id.removesuffix("_rnaalifold")
        for max_iter in MAX_ITERS:
            condition = copy.deepcopy(base)
            condition["id"] = f"{base_name}_rnaalifold_iter{max_iter}"
            condition["label"] = (
                f"{base_name} (RNAalifold profile energy; max-iter={max_iter})"
            )
            condition["phase"] = "linearalifold_profile_max_iter_sweep"
            condition["base_condition"] = base_name
            condition["max_iter"] = max_iter
            condition["linear_profile_energy"] = "rnaalifold"
            condition["args"] = replace_max_iter(condition["args"], max_iter)
            conditions.append(condition)

    output = copy.deepcopy(source)
    output["output_dir"] = "../results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272"
    output["timeout_seconds"] = 3600
    output["conditions"] = conditions
    output["sweep"] = {
        "type": "max_iter",
        "module": "LinearAlifold profile energy (RNAalifold-compatible)",
        "selected_conditions": [condition.removesuffix("_rnaalifold") for condition in SELECTED],
        "values": list(MAX_ITERS),
        "fixed_beams": 100,
        "selection_order": ["mcc", "sci", "sps"],
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"conditions={len(conditions)} datasets=85 runs={len(conditions) * 85}")
    print(OUTPUT)


if __name__ == "__main__":
    main()
