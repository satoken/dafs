#!/usr/bin/env python3
"""Summarize the Murlet dataset1 condition sweep."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "benchmarks/results-murlet-dataset1-condition-sweep-arm64-v272"
CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-condition-sweep-arm64-v272.json"
OUTPUT = ROOT / "docs/murlet_dataset1_condition_sweep_2026-09-18_ja.md"
METRICS = ["sps", "alignment_ppv", "sensitivity", "ppv", "mcc",
           "cbp_f1", "sci", "wall_seconds", "max_rss_kb"]

GROUPS = {
    "ViennaRNA × ProbCons": [
        "ViennaRNA x ProbCons (Alifold)",
        "ViennaRNA x ProbCons (no Alifold)",
    ],
    "CONTRAfold × CONTRAlign": [
        "CONTRAfold x CONTRAlign (Alifold)",
        "CONTRAfold x CONTRAlign (no Alifold)",
    ],
    "lpc × LinearAlign": ["lpc x LinearAlign"],
    "lpv × LinearProbCons": ["lpv x LinearProbCons"],
    "lpc × LinearCONTRAlign": ["lpc x LinearCONTRAlign"],
    "lpv × LinearAlign": ["lpv x LinearAlign"],
}


def mean(rows: list[dict[str, str]], field: str) -> float:
    return sum(float(row[field]) for row in rows) / len(rows)


def rank(row: dict[str, object]) -> tuple[float, float, float]:
    return (-float(row["mcc"]), -float(row["sci"]), -float(row["sps"]))


def fmt(value: object) -> str:
    return f"{float(value):.4f}"


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    metadata = {item["id"]: item for item in config["conditions"]}
    labels = {item["label"]: group
              for group, names in GROUPS.items() for item in
              (metadata[cid] for cid in metadata)
              if item["label"] in names}

    with (SWEEP / "summary.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(config["conditions"]) * 85
    assert all(row["status"] == "success" for row in rows)
    assert all(row["metrics_valid"] == "True" for row in rows)
    assert all(row["quality_complete"] == "True" for row in rows)

    by_condition: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_condition[row["condition"]].append(row)

    condition_rows: list[dict[str, object]] = []
    for condition_id, condition_rows_raw in by_condition.items():
        item = metadata[condition_id]
        aggregate: dict[str, object] = {
            "condition": condition_id,
            "label": item["label"],
            "group": labels[item["label"]],
            "phase": item["phase"],
            "dd": item["ribosum_weight"],
            "final": item["final_ribosum_weight"],
            "args": item["args"],
        }
        aggregate.update({field: mean(condition_rows_raw, field)
                          for field in METRICS})
        condition_rows.append(aggregate)

    best_by_group = {}
    for group in GROUPS:
        candidates = [row for row in condition_rows
                      if row["group"] == group]
        best_by_group[group] = sorted(candidates, key=rank)[0]

    lines = [
        "# Murlet dataset1 条件スイープ結果",
        "",
        "Murlet `dataset1` 全85ケースについて、前回の条件探索範囲をarm64版DAFSで再評価した。8構成（VP/CONTRAfoldのAlifoldあり・なしを含む）、320条件、合計27,200 runである。各条件は各ケース1回、平均値は85ケースのmacro averageである。SCIはViennaRNA 2.7.2で計算した。",
        "",
        "条件の順位は、平均MCCを第一キー、平均SCIを第二キー、平均SPSを第三キーとする辞書式順序で決めた。MCCに許容差は設けていない。全27,200 runで `status=success`、`metrics_valid=True`、`quality_complete=True` だった。",
        "",
        "## 組み合わせ別の最良条件",
        "",
        "| 組み合わせ | 選択条件 | phase | DD RIBOSUM | 最終RIBOSUM | SPS | alignment PPV | 構造感度 | 構造PPV | MCC | CBP F1 | SCI |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for group in GROUPS:
        row = best_by_group[group]
        lines.append("| {} | `{}` | {} | {} | {} | {} | {} | {} | {} | **{}** | {} | {} |".format(
            group, row["condition"], row["phase"], row["dd"], row["final"],
            fmt(row["sps"]), fmt(row["alignment_ppv"]),
            fmt(row["sensitivity"]), fmt(row["ppv"]), fmt(row["mcc"]),
            fmt(row["cbp_f1"]), fmt(row["sci"])))

    lines += [
        "",
        "## phase別のMCC最良条件",
        "",
        "| 組み合わせ | base | parameter | DD RIBOSUM | 最終RIBOSUM |",
        "|---|---|---|---|---|",
    ]
    for group in GROUPS:
        cells = []
        for phase in ["base", "parameter", "dd_ribosum", "final_ribosum"]:
            candidates = [row for row in condition_rows
                          if row["group"] == group and row["phase"] == phase]
            if not candidates:
                cells.append("—")
                continue
            row = sorted(candidates, key=rank)[0]
            cells.append(f"`{row['condition']}` (MCC {fmt(row['mcc'])}, SCI {fmt(row['sci'])}, SPS {fmt(row['sps'])})")
        lines.append(f"| {group} | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} |")

    lines += [
        "",
        "## 選択条件の全引数",
        "",
    ]
    for group in GROUPS:
        row = best_by_group[group]
        lines += [f"### {group} — `{row['condition']}`", "", "```text",
                  " ".join(row["args"]), "```", ""]

    ltf_path = ROOT / "benchmarks/results-linearturbofold-murlet-dataset1-arm64-v272/summary.csv"
    if ltf_path.exists():
        with ltf_path.open(newline="", encoding="utf-8") as stream:
            ltf_rows = list(csv.DictReader(stream))
        lines += [
            "## LinearTurboFold baseline",
            "",
            "LinearTurboFoldは条件スイープ対象ではなく、前回と同じ固定条件（85ケース、hmm/CKY beam 100、3反復、threshold 0.3）をbaselineとして併記する。",
            "",
            "| 方法 | SPS | alignment PPV | 構造感度 | 構造PPV | MCC | CBP F1 | SCI |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        lines.append("| LinearTurboFold | {} | {} | {} | {} | {} | {} | {} |".format(
            *(fmt(mean(ltf_rows, field)) for field in
              ["sps", "alignment_ppv", "sensitivity", "ppv", "mcc", "cbp_f1", "sci"])))

    lines += [
        "",
        "## 実験ファイル",
        "",
        "- [スイープ設定](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/configs/murlet-dataset1-condition-sweep-arm64-v272.json)",
        "- [集計CSV](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-murlet-dataset1-condition-sweep-arm64-v272/summary.csv)",
        "- [LinearTurboFold結果](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-linearturbofold-murlet-dataset1-arm64-v272/summary.csv)",
        "",
    ]
    OUTPUT.write_text("\n".join(lines), encoding="utf-8")
    print(OUTPUT)
    for group, row in best_by_group.items():
        print(group, row["condition"], fmt(row["mcc"]), fmt(row["sci"]), fmt(row["sps"]))


if __name__ == "__main__":
    main()
