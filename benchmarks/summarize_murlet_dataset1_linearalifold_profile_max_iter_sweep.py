#!/usr/bin/env python3
"""Summarize the Murlet dataset1 LinearAlifold profile-energy max-iter sweep."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272.json"
RESULTS = ROOT / "benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272"
SUMMARY = RESULTS / "summary.csv"
OUTPUT_CSV = RESULTS / "max_iter_summary.csv"
OUTPUT_MD = ROOT / "docs/dafs_linearalifold_profile_max_iter_sweep_2026-09-21_ja.md"
BASE_LABELS = {
    "lpc_linear": "lpc × LinearAlign",
    "lpc_contralign": "lpc × LinearCONTRAlign",
    "lpv_linear": "lpv × LinearAlign",
    "lpv_probcons": "lpv × LinearProbConsRNA",
}


def mean(values: list[float]) -> float | None:
    return fmean(values) if values else None


def fmt(value: object, digits: int = 4) -> str:
    if value is None or value == "":
        return "—"
    return f"{float(value):.{digits}f}"


def load_trace(path: Path) -> dict[str, float | int | str | None]:
    events = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    summaries = [event for event in events if event.get("event") == "dd_summary"]
    iterations_by_merge: dict[int, list[dict[str, object]]] = defaultdict(list)
    for event in events:
        if event.get("event") == "dd_iteration":
            iterations_by_merge[int(event["merge_id"])].append(event)

    violations = [int(event.get("violated", 0)) for event in summaries]
    certified: list[int] = []
    for merge_id in sorted(iterations_by_merge):
        value = iterations_by_merge[merge_id][-1].get("certified_track_violated")
        if value is not None:
            certified.append(int(value))

    run_summary = next(
        (event for event in reversed(events) if event.get("event") == "run_summary"),
        {},
    )
    stop_reasons = [str(event.get("stop_reason")) for event in summaries]
    return {
        "merge_count": len(summaries),
        "violations_total": sum(violations),
        "violations_max_merge": max(violations, default=None),
        "certified_violations_total": sum(certified) if certified else None,
        "certified_violations_max_merge": max(certified, default=None),
        "run_dd_iterations": int(run_summary["dd_iterations"])
        if run_summary.get("dd_iterations") is not None else None,
        "stop_reason": Counter(stop_reasons).most_common(1)[0][0]
        if stop_reasons else None,
    }


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    metadata = {str(item["id"]): item for item in config["conditions"]}
    rows: list[dict[str, object]] = []
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            condition = metadata[row["condition"]]
            run_dir = RESULTS / "runs" / row["condition"] / row["dataset"] / "rep-01"
            trace = load_trace(run_dir / "metrics.jsonl")
            rows.append({
                "condition": row["condition"],
                "base_condition": condition["base_condition"],
                "label": condition["label"],
                "max_iter": int(condition["max_iter"]),
                "dataset": row["dataset"],
                "status": row["status"],
                "metrics_valid": row["metrics_valid"],
                "quality_complete": row["quality_complete"],
                "wall_seconds": float(row["wall_seconds"]),
                "internal_seconds": float(row["internal_seconds"]),
                "dd_iterations": int(row["dd_iterations"]),
                "sps": float(row["sps"]),
                "mcc": float(row["mcc"]),
                "sci": float(row["sci"]),
                **trace,
            })

    expected = 4 * 10 * 85
    if len(rows) != expected:
        raise SystemExit(f"expected {expected} rows, got {len(rows)}")

    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["base_condition"]), int(row["max_iter"]))].append(row)

    aggregate: list[dict[str, object]] = []
    for (base_condition, max_iter), items in sorted(grouped.items()):
        final_violations = [int(item["violations_total"]) for item in items]
        merge_counts = [int(item["merge_count"]) for item in items]
        certified = [item["certified_violations_total"] for item in items
                     if item["certified_violations_total"] is not None]
        certified_max = [item["certified_violations_max_merge"] for item in items
                         if item["certified_violations_max_merge"] is not None]
        aggregate.append({
            "base_condition": base_condition,
            "label": BASE_LABELS.get(base_condition, str(items[0]["label"])),
            "max_iter": max_iter,
            "n": len(items),
            "status_success": sum(item["status"] == "success" for item in items),
            "metrics_valid": sum(item["metrics_valid"] == "True" for item in items),
            "quality_complete": sum(item["quality_complete"] == "True" for item in items),
            "mcc": mean([float(item["mcc"]) for item in items]),
            "sci": mean([float(item["sci"]) for item in items]),
            "sps": mean([float(item["sps"]) for item in items]),
            "wall_seconds_mean": mean([float(item["wall_seconds"]) for item in items]),
            "wall_seconds_total": sum(float(item["wall_seconds"]) for item in items),
            "wall_seconds_max": max(float(item["wall_seconds"]) for item in items),
            "internal_seconds_mean": mean([float(item["internal_seconds"]) for item in items]),
            "dd_iterations_mean": mean([float(item["dd_iterations"]) for item in items]),
            "violations_total": sum(final_violations),
            "violations_mean_case": mean([float(value) for value in final_violations]),
            "violations_mean_merge": (
                sum(final_violations) / sum(merge_counts) if sum(merge_counts) else None
            ),
            "violations_max_merge": max(
                int(item["violations_max_merge"])
                for item in items if item["violations_max_merge"] is not None
            ),
            "certified_violations_mean_case": mean([float(value) for value in certified]),
            "certified_violations_mean_merge": (
                sum(float(item["certified_violations_total"]) for item in items
                    if item["certified_violations_total"] is not None)
                / sum(int(item["merge_count"]) for item in items
                      if item["certified_violations_total"] is not None)
                if certified else None
            ),
            "certified_violations_max_merge": max(certified_max, default=None),
        })

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0]))
        writer.writeheader()
        writer.writerows(aggregate)

    status = Counter(str(row["status"]) for row in rows)
    valid = Counter(str(row["metrics_valid"]) for row in rows)
    complete = Counter(str(row["quality_complete"]) for row in rows)
    lines = [
        "# Murlet dataset1 新LinearAlifold max-iter sweep",
        "",
        "新しいLinearAlifoldモジュールとRNAalifold互換プロファイル熱力学モードを固定し、`max-iter=100,200,...,1000`だけを変更した。",
        "Murlet dataset1の85例、4条件×10段階、計3,400 runを評価した。全beamは100固定。",
        "",
        "## 集計方法",
        "",
        "- MCC、SCI、SPSは85例のmacro average。計算時間はrun単位のwall timeとDD内部時間の平均。",
        "- actual DD iterはmetricsの`run_summary.dd_iterations`の平均。",
        "- 制約違反数は各progressive mergeの最終DD反復に記録された`dd_summary.violated`。",
        "  violation/caseは全mergeの違反数を1例単位で合計してから平均し、violation/mergeはmerge数で割った値。",
        "- certified violation/mergeはlinear版のcertified dual trackにおける各merge最終反復値。",
        "",
        f"全runの状態: status={dict(status)}, metrics_valid={dict(valid)}, quality_complete={dict(complete)}。",
        "",
        "## 結果",
        "",
    ]
    by_base: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in aggregate:
        by_base[str(item["base_condition"])].append(item)
    for base_condition, items in by_base.items():
        lines += [f"### {items[0]['label']} (`{base_condition}`)", ""]
        lines += [
            "| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for item in items:
            lines.append(
                f"| {item['max_iter']} | {fmt(item['mcc'])} | {fmt(item['sci'])} | {fmt(item['sps'])} | "
                f"{fmt(item['wall_seconds_mean'], 2)} | {fmt(item['internal_seconds_mean'], 2)} | "
                f"{fmt(item['dd_iterations_mean'], 1)} | {fmt(item['violations_mean_case'], 2)} | "
                f"{fmt(item['violations_mean_merge'], 2)} | {item['violations_max_merge']} | "
                f"{fmt(item['certified_violations_mean_merge'], 2)} |"
            )
        lines.append("")

    lines += ["## 固定した条件", ""]
    for base_condition in by_base:
        condition = next(item for item in metadata.values()
                         if item["base_condition"] == base_condition)
        args = " ".join(str(arg) for arg in condition["args"])
        lines.append(f"- `{base_condition}`: `{args}`（`--max-iter`のみ各行の値に置換）")
    lines += [
        "",
        "## ファイル",
        "",
        "- [再現用設定](../benchmarks/configs/murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272.json)",
        "- [全runのsummary.csv](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/summary.csv)",
        "- [集計CSV](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/max_iter_summary.csv)",
    ]
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUTPUT_CSV)
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
