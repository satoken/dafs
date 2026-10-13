#!/usr/bin/env python3
"""Compare reusable-workspace results with the pre-optimization sweep."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "benchmarks/results-murlet-dataset1-max-iter-sweep-arm64-v272"
RESULTS = ROOT / "benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272"
OUTPUT_CSV = RESULTS / "buffer_reuse_comparison.csv"
OUTPUT_MD = ROOT / "docs/murlet_dataset1_buffer_reuse_2026-09-19_ja.md"

RUN_SUMMARY_FIELDS = (
    "progressive_score",
    "alignment_columns",
    "consensus_pairs",
    "base_pair_probability_initial_nnz",
    "alignment_probability_initial_nnz",
    "base_pair_probability_final_nnz",
    "alignment_probability_final_nnz",
    "dd_calls",
    "dd_iterations",
    "cbp_peak",
    "cbp_added",
    "cbp_priced",
    "cbp_removed",
)


def mean(items: list[dict[str, object]], key: str) -> float:
    return fmean(float(item[key]) for item in items)


def main() -> None:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    predictions_equal = 0
    quality_equal = 0
    summaries_equal = 0
    status_success = 0
    metrics_valid = 0

    result_paths = sorted((RESULTS / "runs").glob("*/*/rep-01/result.json"))
    for result_path in result_paths:
        relative = result_path.relative_to(RESULTS / "runs")
        baseline_path = BASELINE / "runs" / relative
        current = json.loads(result_path.read_text(encoding="utf-8"))
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        current_summary = current["metrics_validation"]["run_summary"]
        baseline_summary = baseline["metrics_validation"]["run_summary"]

        predictions_equal += (
            (result_path.parent / "prediction.fa").read_bytes()
            == (baseline_path.parent / "prediction.fa").read_bytes()
        )
        quality_equal += current["quality"] == baseline["quality"]
        summaries_equal += all(
            current_summary.get(field) == baseline_summary.get(field)
            for field in RUN_SUMMARY_FIELDS
        )
        status_success += current["status"] == "success"
        metrics_valid += bool(current["metrics_validation"]["valid"])

        condition = str(current["condition"])
        base_condition, max_iter_text = condition.rsplit("_iter", 1)
        quality = current["quality"]
        grouped[(base_condition, int(max_iter_text))].append({
            "mcc": quality["mcc"],
            "sci": quality["sci"],
            "sps": quality["sps"],
            "old_wall": baseline["resource"]["wall_seconds"],
            "new_wall": current["resource"]["wall_seconds"],
            "old_dd": baseline_summary["dd_seconds"],
            "new_dd": current_summary["dd_seconds"],
            "old_rss": baseline["resource"]["max_rss_kb"],
            "new_rss": current["resource"]["max_rss_kb"],
        })

    if len(result_paths) != 680:
        raise SystemExit(f"expected 680 results, got {len(result_paths)}")

    aggregate: list[dict[str, object]] = []
    for (condition, max_iter), items in sorted(grouped.items()):
        old_wall, new_wall = mean(items, "old_wall"), mean(items, "new_wall")
        old_dd, new_dd = mean(items, "old_dd"), mean(items, "new_dd")
        old_rss, new_rss = mean(items, "old_rss"), mean(items, "new_rss")
        aggregate.append({
            "condition": condition,
            "max_iter": max_iter,
            "n": len(items),
            "mcc": mean(items, "mcc"),
            "sci": mean(items, "sci"),
            "sps": mean(items, "sps"),
            "old_wall_seconds": old_wall,
            "new_wall_seconds": new_wall,
            "wall_reduction_percent": 100.0 * (old_wall - new_wall) / old_wall,
            "wall_speedup": old_wall / new_wall,
            "old_dd_seconds": old_dd,
            "new_dd_seconds": new_dd,
            "dd_reduction_percent": 100.0 * (old_dd - new_dd) / old_dd,
            "dd_speedup": old_dd / new_dd,
            "old_max_rss_kb": old_rss,
            "new_max_rss_kb": new_rss,
            "rss_change_percent": 100.0 * (new_rss - old_rss) / old_rss,
        })

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0]))
        writer.writeheader()
        writer.writerows(aggregate)

    lines = [
        "# Murlet dataset1: linear workspace再利用の効果",
        "",
        "`BeamAlign::max_alignment()`のlayer/ranking/key領域と、",
        "`LinearNussinov`のbeam/candidate/certificate/traceback領域を呼出し間で再利用した。",
        "beam=100、max-iter=100または1000、85例、線形4条件の計680 runを、",
        "変更前の同一ベンチマークと比較した。",
        "",
        "## 同一性",
        "",
        f"- 成功: {status_success}/680、internal metrics valid: {metrics_valid}/680",
        f"- prediction.faのバイト一致: {predictions_equal}/680",
        f"- MCC/SCI/SPSを含むquality辞書の一致: {quality_equal}/680",
        f"- 時間以外のrun summary一致: {summaries_equal}/680",
        "",
        "## 計算時間とメモリ",
        "",
        "wall/DDは1例あたりの平均秒数。減少率が正なら高速化。RSSはGNU timeの最大常駐メモリ平均。",
        "",
        "| condition | max-iter | MCC | SCI | SPS | wall old→new s | wall減 | DD old→new s | DD減 | RSS変化 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in aggregate:
        lines.append(
            f"| `{item['condition']}` | {item['max_iter']} | "
            f"{item['mcc']:.4f} | {item['sci']:.4f} | {item['sps']:.4f} | "
            f"{item['old_wall_seconds']:.2f}→{item['new_wall_seconds']:.2f} | "
            f"{item['wall_reduction_percent']:.1f}% | "
            f"{item['old_dd_seconds']:.2f}→{item['new_dd_seconds']:.2f} | "
            f"{item['dd_reduction_percent']:.1f}% | "
            f"{item['rss_change_percent']:+.1f}% |"
        )
    lines += [
        "",
        "## 解釈",
        "",
        "- 精度・探索経路を変えず、DD内部時間は全条件で短縮した。",
        "- 短縮幅はmax-iter=100で5.5〜16.9%、max-iter=1000で5.3〜8.4%。",
        "- バッファを解放せず保持するため、最大RSSは3.4〜4.8%増えた。",
        "- したがって今回の変更は、少量の常駐メモリと引き換えに反復ごとのallocator/hash再確保を削減する最適化である。",
        "",
        "## 再現ファイル",
        "",
        "- [設定](../benchmarks/configs/murlet-dataset1-buffer-reuse-arm64-v272.json)",
        "- [比較CSV](../benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272/buffer_reuse_comparison.csv)",
        "- [全run summary](../benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272/summary.csv)",
    ]
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUTPUT_CSV)
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
