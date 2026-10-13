#!/usr/bin/env python3
"""Prepare and summarize the progressive/final LinearAlifold ablation."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-linfold-optimization-iter100-beam100.json"

MODES = ("none", "progressive", "final", "both")
COMBINATIONS = {
    "lpc_linear": ["-a", "LinearAlign", "-s", "lpc",
                   "--ribosum-weight=0.5", "--final-ribosum-weight=0.01",
                   "--align-th=0.02"],
    "lpc_contralign": ["-a", "LinearAlign-CONTRAlign", "-s", "lpc",
                       "--ribosum-weight=0.3", "--final-ribosum-weight=0.005",
                       "--align-pct=0.5", "--align-th=0.02", "--weight=8"],
    "lpv_linear": ["-a", "LinearAlign", "-s", "lpv",
                   "--ribosum-weight=0.3", "--final-ribosum-weight=0",
                   "--align-pct=0.5", "--fold-pct=0.5"],
    "lpv_probcons": ["-a", "LinearAlign-ProbConsRNA", "-s", "lpv",
                     "--ribosum-weight=0.075", "--final-ribosum-weight=0",
                     "--weight=4"],
}
def prepare(args: argparse.Namespace) -> None:
    base = json.loads(BASE_CONFIG.read_text())
    modes = tuple(value for value in args.modes.split(",") if value)
    invalid_modes = sorted(set(modes) - set(MODES))
    if not modes or invalid_modes:
        raise ValueError(f"--modes must contain only {', '.join(MODES)}")
    if args.max_iters:
        iterations = tuple(int(value) for value in args.max_iters.split(",")
                           if value)
    else:
        iterations = (args.max_iter,)
    if not iterations or any(value <= 0 for value in iterations):
        raise ValueError("max-iter values must be positive")
    conditions = []
    for combination, prefix in COMBINATIONS.items():
        for max_iter in iterations:
            common = ["--dynamic-cbp", f"--max-iter={max_iter}",
                      "--align-beam=100", "--align-dd-beam=100",
                      "--linfold-beam=100", "--fold-dd-beam=100",
                      "--fold-final-beam=100"]
            for mode in modes:
                conditions.append({
                    "id": f"{combination}_{mode}_iter{max_iter}",
                    "label": f"{combination} ({mode}, iter={max_iter})",
                    "phase": "stage_ablation_max_iter_sweep",
                    "linear": True,
                    "alifold_stage": mode,
                    "max_iter": max_iter,
                    "args": prefix + [f"--alifold-stages={mode}"] + common,
                })
    config = {
        "binary": str(Path(args.binary).resolve()),
        "conditions": conditions,
        "dataset_filter": {"collections": ["dataset1"]},
        "dataset_manifests": base["dataset_manifests"],
        "environment": {"OMP_NUM_THREADS": "1"},
        "evaluation_scope": "alifold-stage-ablation",
        "evaluation_source_config": str(BASE_CONFIG.resolve()),
        "output_dir": str(Path(args.output_dir).resolve()),
        "random_seed": 42,
        "repetitions": 1,
        "sci": base["sci"],
        "timeout_seconds": 1800,
    }
    ids = [value for value in args.ids.split(",") if value]
    if ids:
        config["dataset_filter"] = {"collections": ["dataset1"], "ids": ids}
    output = Path(args.output)
    output.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n")
    print(output)


def summarize(args: argparse.Namespace) -> None:
    rows = list(csv.DictReader((Path(args.results) / "summary.csv").open()))
    groups = defaultdict(list)
    for row in rows:
        groups[row["condition"]].append(row)
    print("condition,n,mcc,sci,sps,wall_seconds,max_rss_kb,metrics_valid")
    for condition in sorted(groups):
        values = groups[condition]
        mean = lambda key: sum(float(row[key]) for row in values) / len(values)
        print(",".join([
            condition, str(len(values)), f"{mean('mcc'):.6f}",
            f"{mean('sci'):.6f}", f"{mean('sps'):.6f}",
            f"{mean('wall_seconds'):.6f}", f"{mean('max_rss_kb'):.1f}",
            str(all(row["metrics_valid"] == "True" for row in values)),
        ]))


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--binary", required=True)
    prepare_parser.add_argument("--output-dir", required=True)
    prepare_parser.add_argument("--output", required=True)
    prepare_parser.add_argument("--ids", default="")
    prepare_parser.add_argument("--max-iter", type=int, default=100)
    prepare_parser.add_argument("--max-iters", default="")
    prepare_parser.add_argument("--modes", default=','.join(MODES))
    prepare_parser.set_defaults(function=prepare)
    summarize_parser = subparsers.add_parser("summarize")
    summarize_parser.add_argument("--results", required=True)
    summarize_parser.set_defaults(function=summarize)
    arguments = parser.parse_args()
    arguments.function(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
