#!/usr/bin/env python3
"""Compare the old and new LinearAlifold Murlet max-iter sweeps."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[1]
OLD_DIRS = {
    100: Path("/tmp/dafs-alifold-maxiter-100-500"),
    600: Path("/tmp/dafs-alifold-maxiter-600-1000"),
}
NEW_RESULTS = ROOT / "benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272"
NEW_SUMMARY = NEW_RESULTS / "max_iter_summary.csv"
OUTPUT_CSV = NEW_RESULTS / "old_new_comparison.csv"
OUTPUT_MD = ROOT / "docs/dafs_linearalifold_old_new_comparison_2026-09-21_ja.md"

LABELS = {
    "lpc_linear": "lpc × LinearAlign",
    "lpc_contralign": "lpc × LinearCONTRAlign",
    "lpv_linear": "lpv × LinearAlign",
    "lpv_probcons": "lpv × LinearProbConsRNA",
}


def avg(rows: list[dict[str, str]], field: str) -> float:
    return fmean(float(row[field]) for row in rows)


def trace(path: Path) -> tuple[int, int, int | None, int | None]:
    merge_count = 0
    violation_total = 0
    violation_max: int | None = None
    certified_by_merge: dict[int, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        event = json.loads(line)
        if event.get("event") == "dd_summary":
            value = int(event.get("violated", 0))
            merge_count += 1
            violation_total += value
            violation_max = value if violation_max is None else max(violation_max, value)
        elif event.get("event") == "dd_iteration":
            value = event.get("certified_track_violated")
            if value is not None:
                certified_by_merge[int(event["merge_id"])] = int(value)
    certified_total = sum(certified_by_merge.values()) if certified_by_merge else None
    certified_max = max(certified_by_merge.values()) if certified_by_merge else None
    return merge_count, violation_total, certified_total, certified_max


def old_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for first_iter, result_dir in OLD_DIRS.items():
        with (result_dir / "summary.csv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                condition = row["condition"]
                if "_both_iter" not in condition:
                    continue
                base, iter_text = condition.split("_both_iter")
                max_iter = int(iter_text)
                if (max_iter >= 600) != (first_iter >= 600):
                    continue
                run_dir = result_dir / "runs" / condition / row["dataset"] / "rep-01"
                merge_count, violations, certified, certified_max = trace(run_dir / "metrics.jsonl")
                rows.append({
                    "base_condition": base,
                    "max_iter": max_iter,
                    "dataset": row["dataset"],
                    "status": row["status"],
                    "mcc": float(row["mcc"]),
                    "sci": float(row["sci"]),
                    "sps": float(row["sps"]),
                    "wall_seconds": float(row["wall_seconds"]),
                    "internal_seconds": float(row["internal_seconds"]),
                    "dd_iterations": int(row["dd_iterations"]),
                    "merge_count": merge_count,
                    "violations_total": violations,
                    "certified_violations_total": certified,
                    "certified_violations_max_merge": certified_max,
                })
    if len(rows) != 4 * 10 * 85:
        raise SystemExit(f"expected 3400 old rows, got {len(rows)}")
    return rows


def aggregate_old(rows: list[dict[str, object]]) -> dict[tuple[str, int], dict[str, object]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["base_condition"]), int(row["max_iter"]))].append(row)
    result: dict[tuple[str, int], dict[str, object]] = {}
    for key, items in grouped.items():
        merge_count = sum(int(item["merge_count"]) for item in items)
        certified_items = [item for item in items if item["certified_violations_total"] is not None]
        certified_merges = sum(int(item["merge_count"]) for item in certified_items)
        result[key] = {
            "mcc": fmean(float(item["mcc"]) for item in items),
            "sci": fmean(float(item["sci"]) for item in items),
            "sps": fmean(float(item["sps"]) for item in items),
            "wall_seconds": fmean(float(item["wall_seconds"]) for item in items),
            "internal_seconds": fmean(float(item["internal_seconds"]) for item in items),
            "dd_iterations": fmean(float(item["dd_iterations"]) for item in items),
            "violations_mean_case": fmean(float(item["violations_total"]) for item in items),
            "violations_mean_merge": sum(int(item["violations_total"]) for item in items) / merge_count,
            "certified_violations_mean_merge": (
                sum(int(item["certified_violations_total"]) for item in certified_items) / certified_merges
                if certified_items else None
            ),
        }
    return result


def fmt(value: object, digits: int = 4) -> str:
    if value is None or value == "":
        return "—"
    return f"{float(value):.{digits}f}"


def main() -> None:
    old = aggregate_old(old_rows())
    with NEW_SUMMARY.open(newline="", encoding="utf-8") as stream:
        new = {(row["base_condition"], int(row["max_iter"])): row
               for row in csv.DictReader(stream)}

    rows: list[dict[str, object]] = []
    for base_condition in LABELS:
        for max_iter in range(100, 1001, 100):
            old_row = old[(base_condition, max_iter)]
            new_row = new[(base_condition, max_iter)]
            item: dict[str, object] = {
                "base_condition": base_condition,
                "max_iter": max_iter,
            }
            new_field_names = {
                "wall_seconds": "wall_seconds_mean",
                "internal_seconds": "internal_seconds_mean",
                "dd_iterations": "dd_iterations_mean",
            }
            for field in (
                "mcc", "sci", "sps", "wall_seconds", "internal_seconds",
                "dd_iterations", "violations_mean_case", "violations_mean_merge",
                "certified_violations_mean_merge",
            ):
                old_value = old_row[field]
                source_field = new_field_names.get(field, field)
                new_value = float(new_row[source_field]) if new_row[source_field] not in ("", "None") else None
                item[f"old_{field}"] = old_value
                item[f"new_{field}"] = new_value
                item[f"delta_{field}"] = new_value - old_value if new_value is not None else None
            rows.append(item)

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# LinearAlifold新旧比較（Murlet dataset1）",
        "",
        "同じMurlet dataset1の85例、4組合せ、beam=100、max-iter=100〜1000で、旧実装と新実装を比較した。",
        "",
        "- 旧: 2026-09-20実行、`--alifold-stages=both`、プロファイル熱力学は従来実装。",
        "- 新: 2026-09-21実行、`--alifold`、`--linear-profile-energy=rnaalifold`。",
        "- 旧バイナリSHA256: `da88c827aabf763aca8dc3c3576f9d2735eed35e85c467f64a4593e8abbebd04`。",
        "- 新バイナリSHA256: `ad60c39ef6def8bd51aa0fb7d92649137cac2aba067bd3d159269160190ccaf6`。",
        "- 新旧でバイナリは異なるため、差はプロファイル熱力学モード変更を主因と考える比較であり、他のバイナリ差分を完全に分離した因果比較ではない。",
        "",
        "## max-iterごとの比較",
        "",
    ]
    for base_condition in LABELS:
        lines += [f"### {LABELS[base_condition]} (`{base_condition}`)", ""]
        lines += [
            "| max-iter | 旧 MCC / SCI / SPS | 新 MCC / SCI / SPS | Δ MCC / Δ SCI / Δ SPS | 旧/新 wall秒 | 旧/新 violation/merge |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
        for row in rows:
            if row["base_condition"] != base_condition:
                continue
            lines.append(
                f"| {row['max_iter']} | {fmt(row['old_mcc'])} / {fmt(row['old_sci'])} / {fmt(row['old_sps'])} | "
                f"{fmt(row['new_mcc'])} / {fmt(row['new_sci'])} / {fmt(row['new_sps'])} | "
                f"{fmt(row['delta_mcc'], 4)} / {fmt(row['delta_sci'], 4)} / {fmt(row['delta_sps'], 4)} | "
                f"{fmt(row['old_wall_seconds'], 2)} / {fmt(row['new_wall_seconds'], 2)} | "
                f"{fmt(row['old_violations_mean_merge'], 2)} / {fmt(row['new_violations_mean_merge'], 2)} |"
            )
        lines.append("")

    lines += [
        "## max-iter=1000での要約",
        "",
        "| 組合せ | ΔMCC | ΔSCI | ΔSPS | wall増加 | violation/mergeの変化 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for base_condition in LABELS:
        row = next(item for item in rows if item["base_condition"] == base_condition and item["max_iter"] == 1000)
        wall_ratio = float(row["new_wall_seconds"]) / float(row["old_wall_seconds"])
        lines.append(
            f"| {LABELS[base_condition]} | {fmt(row['delta_mcc'], 4)} | {fmt(row['delta_sci'], 4)} | "
            f"{fmt(row['delta_sps'], 4)} | {wall_ratio:.2f}× | "
            f"{fmt(row['old_violations_mean_merge'], 2)} → {fmt(row['new_violations_mean_merge'], 2)} |"
        )
    lines += [
        "",
        "## ファイル",
        "",
        "- [新実装の再現用設定](../benchmarks/configs/murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272.json)",
        "- [新実装の全run summary](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/summary.csv)",
        "- [新旧比較CSV](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/old_new_comparison.csv)",
    ]
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUTPUT_CSV)
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
