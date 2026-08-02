#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "benchmarks"))

from score import parse_alignment, score, structure_pairs  # noqa: E402


class BenchmarkScorerTest(unittest.TestCase):
    def test_identical_stockholm_and_dafs_output_score_one(self):
        reference = parse_alignment(REPOSITORY / "tests/data/tiny_reference.sto")
        prediction = REPOSITORY / "tests/data/tiny_prediction.fa"
        metrics = score(parse_alignment(prediction), reference)
        self.assertEqual(metrics["sps"], 1.0)
        self.assertEqual(metrics["alignment_ppv"], 1.0)
        self.assertEqual(metrics["sensitivity"], 1.0)
        self.assertEqual(metrics["ppv"], 1.0)
        self.assertEqual(metrics["mcc"], 1.0)
        self.assertEqual(metrics["cbp_f1"], 1.0)

    def test_missing_structure_is_reported_without_fabricated_scores(self):
        reference = parse_alignment(REPOSITORY / "tests/data/tiny_reference.sto")
        prediction = REPOSITORY / "tests/data/tiny_prediction_no_structure.fa"
        metrics = score(parse_alignment(prediction), reference)
        self.assertEqual(metrics["sps"], 1.0)
        self.assertIsNone(metrics["mcc"])
        self.assertIsNone(metrics["cbp_f1"])

    def test_wuss_pseudoknot_symbols(self):
        self.assertEqual(structure_pairs("<A..>a"), {(0, 4), (1, 5)})


if __name__ == "__main__":
    unittest.main()
