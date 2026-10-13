#!/usr/bin/env python3
"""Prepare and summarize the nonlinear Alifold max-iter sweep."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-linfold-optimization-iter100-beam100.json"
ITERATIONS = tuple(range(100, 1001, 100))
COMBINATIONS = {
    "contrafold_contralign": {
        "label": "CONTRAfold x CONTRAlign",
        "prefix": ["-a", "CONTRAlign", "-s", "CONTRAfold",
                   "--ribosum-weight=0.1", "--final-ribosum-weight=0",
                   "--fold-th=0.3", "--weight=2"],
    },
    "vienna_probcons": {
        "label": "ViennaRNA x ProbCons",
        "prefix": ["-a", "ProbCons", "-s", "Vienna",
                   "--ribosum-weight=0.02", "--final-ribosum-weight=0",
                   "--fold-th=0.3", "--weight=2"],
    },
}


def prepare(args: argparse.Namespace) -> None:
    base = json.loads(BASE_CONFIG.read_text())
    conditions = []
    for combination, definition in COMBINATIONS.items():
        for max_iter in ITERATIONS:
            for mode, alifold_flag in (("none", "--no-alifold"),
                                       ("both", "--alifold")):
                conditions.append({
                    "id": f"{combination}_{mode}_iter{max_iter}",
                    "label": f"{definition['label']} ({mode}, iter={max_iter})",
                    "phase": "nonlinear_alifold_max_iter_sweep",
                    "linear": False,
                    "alifold_stage": mode,
                    "max_iter": max_iter,
                    "args": definition["prefix"] + [alifold_flag,
                             "--dynamic-cbp", f"--max-iter={max_iter}"],
                })
    config = {
        "binary": str(Path(args.binary).resolve()),
        "conditions": conditions,
        "dataset_filter": {"collections": ["dataset1"]},
        "dataset_manifests": base["dataset_manifests"],
        "environment": {"OMP_NUM_THREADS": "1"},
        "evaluation_scope": "nonlinear-alifold-max-iter-sweep",
        "evaluation_source_config": str(BASE_CONFIG.resolve()),
        "output_dir": str(Path(args.output_dir).resolve()),
        "random_seed": 42,
        "repetitions": 1,
        "sci": base["sci"],
        "timeout_seconds": 3600,
    }
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
    prepare_parser.set_defaults(function=prepare)
    summarize_parser = subparsers.add_parser("summarize")
    summarize_parser.add_argument("--results", required=True)
    summarize_parser.set_defaults(function=summarize)
    arguments = parser.parse_args()
    arguments.function(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
