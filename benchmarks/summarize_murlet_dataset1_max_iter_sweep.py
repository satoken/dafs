#!/usr/bin/env python3
"""Summarize the Murlet dataset1 max-iter sweep."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-max-iter-sweep-arm64-v272.json"
RESULTS = ROOT / "benchmarks/results-murlet-dataset1-max-iter-sweep-arm64-v272"
SUMMARY = RESULTS / "summary.csv"
OUTPUT_CSV = RESULTS / "max_iter_summary.csv"
OUTPUT_MD = ROOT / "docs/murlet_dataset1_max_iter_sweep_2026-09-18_ja.md"


def number(values: list[float]) -> float | None:
    return fmean(values) if values else None


def fmt(value: object, digits: int = 4) -> str:
    if value is None or value == "":
        return "—"
    return f"{float(value):.{digits}f}"


def load_trace(path: Path) -> dict[str, float | int | str | None]:
    events: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            events.append(json.loads(line))

    summaries = [event for event in events if event.get("event") == "dd_summary"]
    iterations_by_merge: dict[int, list[dict[str, object]]] = defaultdict(list)
    for event in events:
        if event.get("event") == "dd_iteration":
            iterations_by_merge[int(event["merge_id"])] .append(event)

    final_violations = [int(event.get("violated", 0)) for event in summaries]
    certified_final: list[int] = []
    for merge_id in sorted(iterations_by_merge):
        event = iterations_by_merge[merge_id][-1]
        value = event.get("certified_track_violated")
        if value is not None:
            certified_final.append(int(value))

    run_summary = next(
        (event for event in reversed(events) if event.get("event") == "run_summary"),
        {},
    )
    stop_reasons = [str(event.get("stop_reason")) for event in summaries]
    return {
        "merge_count": len(summaries),
        "violations_total": sum(final_violations),
        "violations_max_merge": max(final_violations, default=None),
        "certified_violations_total": sum(certified_final) if certified_final else None,
        "certified_violations_max_merge": max(certified_final, default=None),
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
            result_path = RESULTS / "runs" / row["condition"] / row["dataset"] / "rep-01" / "result.json"
            result = json.loads(result_path.read_text(encoding="utf-8"))
            trace = load_trace(result_path.parent / "metrics.jsonl")
            if not bool(condition.get("linear")):
                trace["certified_violations_total"] = None
                trace["certified_violations_max_merge"] = None
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

    if len(rows) != 5100:
        raise SystemExit(f"expected 5100 rows, got {len(rows)}")

    aggregate: list[dict[str, object]] = []
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["base_condition"]), int(row["max_iter"]))].append(row)
    for (base_condition, max_iter), items in sorted(grouped.items()):
        label = str(items[0]["label"])
        final_violation_values = [int(item["violations_total"]) for item in items]
        merge_counts = [int(item["merge_count"]) for item in items]
        certified_values = [item["certified_violations_total"] for item in items
                            if item["certified_violations_total"] is not None]
        certified_max_values = [item["certified_violations_max_merge"] for item in items
                                if item["certified_violations_max_merge"] is not None]
        aggregate.append({
            "base_condition": base_condition,
            "label": label,
            "max_iter": max_iter,
            "n": len(items),
            "status_success": sum(item["status"] == "success" for item in items),
            "metrics_valid": sum(item["metrics_valid"] == "True" for item in items),
            "quality_complete": sum(item["quality_complete"] == "True" for item in items),
            "mcc": number([float(item["mcc"]) for item in items]),
            "sci": number([float(item["sci"]) for item in items]),
            "sps": number([float(item["sps"]) for item in items]),
            "wall_seconds_mean": number([float(item["wall_seconds"]) for item in items]),
            "wall_seconds_total": sum(float(item["wall_seconds"]) for item in items),
            "wall_seconds_max": max(float(item["wall_seconds"]) for item in items),
            "internal_seconds_mean": number([float(item["internal_seconds"]) for item in items]),
            "dd_iterations_mean": number([float(item["dd_iterations"]) for item in items]),
            "violations_total": sum(final_violation_values),
            "violations_mean_case": number([float(value) for value in final_violation_values]),
            "violations_mean_merge": (
                sum(final_violation_values) / sum(merge_counts)
                if sum(merge_counts) else None
            ),
            "violations_max_merge": max(
                int(item["violations_max_merge"])
                for item in items if item["violations_max_merge"] is not None
            ),
            "certified_violations_mean_case": number(
                [float(value) for value in certified_values]
            ),
            "certified_violations_mean_merge": (
                sum(float(item["certified_violations_total"])
                    for item in items if item["certified_violations_total"] is not None)
                / sum(int(item["merge_count"]) for item in items
                      if item["certified_violations_total"] is not None)
                if certified_values else None
            ),
            "certified_violations_max_merge": max(certified_max_values, default=None),
        })

    fieldnames = list(aggregate[0])
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(aggregate)

    status = Counter(str(row["status"]) for row in rows)
    valid = Counter(str(row["metrics_valid"]) for row in rows)
    complete = Counter(str(row["quality_complete"]) for row in rows)
    lines = [
        "# Murlet dataset1 max-iter sweep",
        "",
        "前回のMCC→SCI→SPS選択条件を固定し、`max-iter=100,200,...,1000`だけを変更した。",
        "Murlet dataset1の85例、6条件×10段階、計5,100 runを評価した。",
        "",
        "## 集計方法",
        "",
        "- MCC、SCI、SPSは85例のmacro average。",
        "- 計算時間はGNU `time`のwall time（秒）。`internal DD`はmetrics JSONLのDD内部時間。",
        "- 制約違反数は、各progressive mergeの最終DD反復に記録された `dd_summary.violated`。",
        "  `violation mean/case`は全mergeの値を1例単位で合計してから平均し、",
        "  `violation mean/merge`は全mergeで割った値。",
        "- 線形版では、独立したcertified dual trackの最終反復値も併記する。",
        "",
        f"全runの状態: status={dict(status)}, metrics_valid={dict(valid)}, quality_complete={dict(complete)}。",
        "",
        "## 結果",
        "",
    ]
    grouped_aggregate: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in aggregate:
        grouped_aggregate[str(item["base_condition"])].append(item)
    for base_condition, items in grouped_aggregate.items():
        lines += [f"### {items[0]['label']} (`{base_condition}`)", ""]
        lines += [
            "| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for item in items:
            lines.append(
                f"| {item['max_iter']} | {fmt(item['mcc'])} | {fmt(item['sci'])} | "
                f"{fmt(item['sps'])} | {fmt(item['wall_seconds_mean'], 2)} | "
                f"{fmt(item['internal_seconds_mean'], 2)} | {fmt(item['dd_iterations_mean'], 1)} | "
                f"{fmt(item['violations_mean_case'], 2)} | {fmt(item['violations_mean_merge'], 2)} | "
                f"{item['violations_max_merge']} | {fmt(item['certified_violations_mean_merge'], 2)} |"
            )
        lines.append("")

    lines += [
        "## 固定したDAFSパラメータ",
        "",
    ]
    base_conditions = {str(item["base_condition"]): item for item in metadata.values()}
    for base_condition in grouped_aggregate:
        item = base_conditions[base_condition]
        args = " ".join(str(arg) for arg in item["args"])
        lines += [f"- `{base_condition}`: `{args}`（`--max-iter`のみ各行の値に置換）"]
    lines += [
        "",
        "## ファイル",
        "",
        "- [再現用設定](../benchmarks/configs/murlet-dataset1-max-iter-sweep-arm64-v272.json)",
        "- [全runのsummary.csv](../benchmarks/results-murlet-dataset1-max-iter-sweep-arm64-v272/summary.csv)",
        "- [集計CSV](../benchmarks/results-murlet-dataset1-max-iter-sweep-arm64-v272/max_iter_summary.csv)",
    ]
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUTPUT_CSV)
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
