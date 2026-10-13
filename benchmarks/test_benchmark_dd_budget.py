#!/usr/bin/env python3
"""Focused regression tests for benchmark_dd_budget.py."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch


BENCHMARKS = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCHMARKS))
import benchmark_dd_budget as dd  # noqa: E402


def record(condition: str, dataset: str, mcc: float | None = 0.5) -> dict:
    row = {
        "condition": condition,
        "dataset": dataset,
        "repetition": 1,
        "status": "success",
        "metrics_valid": True,
        "quality_complete": True,
    }
    for field in dd.SUMMARY_FIELDS:
        row[field] = 1.0
    row["mcc"] = mcc
    row["sci"] = 0.5
    row["sps"] = 0.5
    return row


def synthetic_comparison(
    label: str,
    *,
    near_or_better: bool = True,
    per_condition_near_or_better: bool = True,
    sci_delta: float = 0.02,
) -> dict:
    side = {
        "n": 48,
        "metrics_valid": 48,
        "quality_complete": 48,
        "status_counts": {"success": 48},
    }
    return {
        "label": label,
        "paired_count": 48,
        "missing_candidate": 0,
        "missing_reference": 0,
        "candidate": side,
        "reference": side,
        "thresholds": {
            "loose": {"near_or_better": near_or_better},
            "strict": {"near_or_better": near_or_better},
        },
        "mean_deltas": {"mcc": 0.0, "sci": sci_delta, "sps": 0.0},
        "distributions": {field: {"n": 48} for field in dd.QUALITY_FIELDS},
        "per_condition": [
            {
                "condition": condition,
                "paired_count": 12,
                "thresholds": {
                    "loose": {"near_or_better": per_condition_near_or_better}
                },
            }
            for condition in dd.BASE_CONDITIONS
        ],
    }


class DDBudgetTests(unittest.TestCase):
    def test_boolean_flags_do_not_consume_following_argument(self) -> None:
        args = [
            "--dd-projected-norm",
            "--max-iter=100",
            "--seed",
            "7",
            "--refinement",
            "1",
        ]
        self.assertEqual(
            dd.remove_option(args, set(dd.FLAG_NAMES), {"--seed", "--refinement"}),
            ["--max-iter=100"],
        )

    def test_extrema_ignore_missing_deltas(self) -> None:
        candidate = [
            record("cl_linear_ath02", "missing", None),
            record("cl_linear_ath02", "present", 0.4),
        ]
        reference = [
            record("cl_linear_ath02", "missing", None),
            record("cl_linear_ath02", "present", 0.3),
        ]
        report = dd.paired_summary(candidate, reference, "synthetic")
        self.assertEqual(
            report["worst_cases"]["mcc"]["min_delta"]["dataset"], "present"
        )
        self.assertEqual(
            report["worst_cases"]["mcc"]["max_delta"]["dataset"], "present"
        )

    def test_condition_pairing_filters_340_reference_rows_to_85(self) -> None:
        datasets = [f"dataset-{index:03d}" for index in range(85)]
        candidate = [record("lpc_contralign_ath02", dataset) for dataset in datasets]
        reference = []
        for condition in dd.BASE_CONDITIONS:
            reference.extend(
                record(
                    condition,
                    dataset,
                    0.25 if condition == "lpc_contralign_ath02" else 0.75,
                )
                for dataset in datasets
            )
        report = dd.paired_summary(candidate, reference, "85-v-340")
        self.assertEqual(report["candidate_count"], 85)
        self.assertEqual(report["reference_input_count"], 340)
        self.assertEqual(report["reference_count"], 85)
        self.assertEqual(report["reference"]["n"], 85)
        self.assertAlmostEqual(report["reference"]["means"]["mcc"], 0.25)
        self.assertEqual(report["paired_count"], 85)
        self.assertEqual(report["missing_candidate"], 0)
        self.assertEqual(report["missing_reference"], 0)

    def test_selection_parses_four_part_label_and_allows_improvement(self) -> None:
        report = {
            "comparisons": [
                synthetic_comparison(
                    "pilot:projected-norm:iter100:vs-baseline500",
                    sci_delta=0.02,
                ),
                synthetic_comparison(
                    "pilot:column-bound:iter100:vs-baseline500",
                    near_or_better=False,
                ),
            ]
        }
        selection = dd.selection_from_report(report)
        self.assertEqual(selection["selected"]["mode"], "projected-norm")
        self.assertEqual(selection["selected"]["budget"], 100)

    def test_selection_requires_every_combo_to_pass(self) -> None:
        report = {
            "comparisons": [
                synthetic_comparison(
                    "pilot:projected-norm:iter100:vs-baseline500",
                    per_condition_near_or_better=False,
                )
            ]
        }
        self.assertIsNone(dd.selection_from_report(report)["selected"])

    def test_baseline_comparison_is_not_duplicated(self) -> None:
        rows = []
        for mode in dd.MODES:
            for budget in (100, 500):
                rows.append(record("cl_linear_ath02", f"{mode}-{budget}"))
                rows[-1].update({"scope": "pilot", "mode": mode, "budget": budget})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(dd, "load_records", return_value=rows):
                dd.summarize(Namespace(results_root=root))
            report = json.loads((root / "dd_budget_report.json").read_text())
        labels = [item["label"] for item in report["comparisons"]]
        self.assertEqual(labels.count("pilot:baseline:iter100:vs-baseline500"), 1)


if __name__ == "__main__":
    unittest.main()
