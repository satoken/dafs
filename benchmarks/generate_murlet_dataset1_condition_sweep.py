#!/usr/bin/env python3
"""Generate the full Murlet dataset1 DAFS condition sweep configuration."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / "benchmarks" / "configs" /
          "murlet-dataset1-condition-sweep-arm64-v272.json")


def set_option(args: list[str], name: str, value: object) -> list[str]:
    prefix = name + "="
    result = [arg for arg in args if not arg.startswith(prefix)]
    result.append(f"{name}={value:g}" if isinstance(value, float)
                  else f"{name}={value}")
    return result


def replace_ribosum(args: list[str], dd: float, final: float) -> list[str]:
    return set_option(set_option(args, "--ribosum-weight", dd),
                      "--final-ribosum-weight", final)


def value_id(value: float) -> str:
    return (f"{value:g}".replace(".", "p")
            .replace("-", "m"))


def main() -> None:
    beams = ["--align-beam=100", "--align-dd-beam=100",
             "--linfold-beam=100", "--fold-dd-beam=100",
             "--fold-final-beam=100"]

    # These are the operating points selected in the preceding 13-case study.
    # The two additional no-/alifold variants are retained because they were
    # explicit alternatives in the preceding condition sweeps.
    bases = [
        ("vp_alifold", "ViennaRNA x ProbCons (Alifold)", 0.02, 0.0,
         ["-a", "ProbCons", "-s", "Vienna", "--alifold", "--dynamic-cbp",
          "--ribosum-weight=0.02", "--final-ribosum-weight=0",
          "--fold-th=0.3", "--max-iter=80"]),
        ("vp_noalifold", "ViennaRNA x ProbCons (no Alifold)", 0.3, 0.0,
         ["-a", "ProbCons", "-s", "Vienna", "--no-alifold", "--dynamic-cbp",
          "--ribosum-weight=0.3", "--final-ribosum-weight=0",
          "--weight=2", "--max-iter=80"]),
        ("cc_alifold", "CONTRAfold x CONTRAlign (Alifold)", 0.1, 0.0,
         ["-a", "CONTRAlign", "-s", "CONTRAfold", "--alifold",
          "--dynamic-cbp", "--ribosum-weight=0.1",
          "--final-ribosum-weight=0", "--weight=2", "--max-iter=80"]),
        ("cc_noalifold", "CONTRAfold x CONTRAlign (no Alifold)", 0.1, 0.0,
         ["-a", "CONTRAlign", "-s", "CONTRAfold", "--no-alifold",
          "--dynamic-cbp", "--ribosum-weight=0.1",
          "--final-ribosum-weight=0", "--fold-pct=0", "--max-iter=80"]),
        ("cl_linear", "lpc x LinearAlign", 0.5, 0.01,
         ["-a", "LinearAlign", "-s", "lpc", "--no-alifold",
          "--dynamic-cbp", "--ribosum-weight=0.5",
          "--final-ribosum-weight=0.01", "--align-th=0.005",
          "--max-iter=80", *beams]),
        ("lpv_probcons", "lpv x LinearProbCons", 0.3, 0.0,
         ["-a", "LinearAlign-ProbConsRNA", "-s", "lpv", "--no-alifold",
          "--dynamic-cbp", "--ribosum-weight=0.3",
          "--final-ribosum-weight=0", "--weight=4", "--max-iter=80",
          *beams]),
        ("lpc_contralign", "lpc x LinearCONTRAlign", 0.3, 0.005,
         ["-a", "LinearAlign-CONTRAlign", "-s", "lpc", "--no-alifold",
          "--dynamic-cbp", "--ribosum-weight=0.3",
          "--final-ribosum-weight=0.005", "--align-pct=0.5",
          "--max-iter=80", *beams]),
        ("lpv_linear", "lpv x LinearAlign", 0.3, 0.0,
         ["-a", "LinearAlign", "-s", "lpv", "--no-alifold",
          "--dynamic-cbp", "--ribosum-weight=0.3",
          "--final-ribosum-weight=0", "--align-pct=0.5",
          "--max-iter=80", *beams]),
    ]

    conditions: list[dict[str, object]] = []
    seen_commands: set[tuple[str, ...]] = set()

    def add(condition_id: str, label: str, args: list[str], dd: float,
            final: float, phase: str) -> None:
        key = tuple(args)
        if key in seen_commands:
            return
        seen_commands.add(key)
        conditions.append({
            "id": condition_id,
            "label": label,
            "phase": phase,
            "linear": any("Linear" in arg for arg in args),
            "ribosum_weight": dd,
            "final_ribosum_weight": final,
            "args": args,
        })

    # The same one-factor grid used in the previous parameter studies.
    perturbations = [
        ("w2", "--weight", 2), ("w8", "--weight", 8),
        ("apct0", "--align-pct", 0), ("apct05", "--align-pct", 0.5),
        ("fpct0", "--fold-pct", 0), ("fpct05", "--fold-pct", 0.5),
        ("ath005", "--align-th", 0.005), ("ath02", "--align-th", 0.02),
        ("fth01", "--fold-th", 0.1), ("fth03", "--fold-th", 0.3),
    ]
    dd_values = [0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1,
                 0.15, 0.2, 0.3, 0.5, 0.75, 1, 2]
    final_values = dd_values[:]

    for prefix, label, dd, final, base_args in bases:
        add(f"{prefix}_base", label, base_args, dd, final, "base")
        for suffix, option, value in perturbations:
            args = set_option(base_args, option, value)
            add(f"{prefix}_{suffix}", label, args, dd, final, "parameter")

        for value in dd_values:
            args = replace_ribosum(base_args, value, 0)
            add(f"{prefix}_dd_{value_id(value)}", label, args, value, 0,
                "dd_ribosum")

        for value in final_values:
            args = replace_ribosum(base_args, dd, value)
            add(f"{prefix}_final_{value_id(value)}", label, args, dd, value,
                "final_ribosum")

    config = {
        "binary": "../../build-arm64/src/dafs",
        "output_dir": "../results-murlet-dataset1-condition-sweep-arm64-v272",
        "timeout_seconds": 1800,
        "repetitions": 1,
        "random_seed": 42,
        "environment": {"OMP_NUM_THREADS": "1"},
        "sci": {
            "bin_dir": "/home/sato-lab.org/satoken/.local/app/viennarna-2.7.2-arm64/bin",
            "version": "2.7.2",
        },
        "dataset_manifests": ["../data/murlet/manifest.json"],
        "dataset_filter": {"collections": ["dataset1"]},
        "conditions": conditions,
    }
    CONFIG.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"conditions={len(conditions)} datasets=85 runs={len(conditions) * 85}")
    print(CONFIG)


if __name__ == "__main__":
    main()
