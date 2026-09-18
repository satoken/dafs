#!/usr/bin/env python3

import stat
import sys
import tempfile
from textwrap import dedent
import unittest
from collections import OrderedDict
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "benchmarks"))

from score import (SCIConfig, calculate_sci, parse_alignment,
                   parse_linearturbofold, match_prediction_names, score,
                   StructuralAlignment, structure_pairs)  # noqa: E402


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

    def test_match_prediction_names_restores_linear_turbofold_separator(self):
        prediction = StructuralAlignment(
            OrderedDict((("seq1_1", "GGGAAACCC"), ("seq2_1", "GGGAAUCCC"))),
            None,
        )
        reference = StructuralAlignment(
            OrderedDict((("seq1/1", "GGGAAACCC"), ("seq2/1", "GGGAAUCCC"))),
            None,
        )
        matched = match_prediction_names(prediction, reference)
        self.assertEqual(list(matched.sequences), ["seq1/1", "seq2/1"])
        self.assertEqual(score(matched, reference)["sps"], 1.0)

    def test_linearturbofold_alignment_and_individual_structures(self):
        reference = parse_alignment(REPOSITORY / "tests/data/tiny_reference.sto")
        with tempfile.TemporaryDirectory(prefix="ltf-score-test-") as directory:
            output_dir = Path(directory)
            (output_dir / "output.aln").write_text(dedent("""\
                >seq1
                GGGAAACCC
                >seq2
                GGGAAUCCC
            """), encoding="utf-8")
            for index, name, sequence in (
                    (1, "seq1", "GGGAAACCC"),
                    (2, "seq2", "GGGAAUCCC")):
                (output_dir / f"{index}_{name}.db").write_text(
                    f">{name}\n{sequence}\n(((...)))\n", encoding="utf-8")
            prediction = parse_linearturbofold(
                output_dir / "output.aln", output_dir)
            metrics = score(prediction, reference)

        self.assertEqual(metrics["sps"], 1.0)
        self.assertEqual(metrics["sensitivity"], 1.0)
        self.assertEqual(metrics["mcc"], 1.0)
        self.assertEqual(metrics["cbp_f1"], 1.0)

    def test_sci_uses_versioned_rnaalifold_and_rnafold(self):
        prediction = parse_alignment(REPOSITORY / "tests/data/tiny_prediction.fa")
        with tempfile.TemporaryDirectory(prefix="dafs-sci-test-") as directory:
            tools = Path(directory)
            rnaalifold = tools / "RNAalifold"
            rnafold = tools / "RNAfold"
            rnaalifold.write_text(dedent(f"""\
                #!{sys.executable}
                import sys
                if "--version" in sys.argv or "-V" in sys.argv:
                    print("RNAalifold 2.4.18")
                    raise SystemExit(0)
                print("2 sequences; length of alignment 9")
                print(">consensus")
                print("GGGAAACCC")
                print("(((...))) (-6.00 = -5.00 + -1.00)")
            """), encoding="utf-8")
            rnafold.write_text(dedent(f"""\
                #!{sys.executable}
                import sys
                if "--version" in sys.argv or "-V" in sys.argv:
                    print("RNAfold 2.4.18")
                    raise SystemExit(0)
                count = sum(line.startswith(">") for line in sys.stdin)
                for index in range(count):
                    print(f">seq{{index + 1}}")
                    print("GGGAAACCC")
                    print("(((...))) (-4.00)")
            """), encoding="utf-8")
            for tool in (rnaalifold, rnafold):
                tool.chmod(tool.stat().st_mode | stat.S_IXUSR)

            config = SCIConfig(rnaalifold, rnafold, "2.4.18", 10)
            metrics = calculate_sci(prediction, config)

        self.assertAlmostEqual(metrics["sci"], 1.5)
        self.assertEqual(metrics["sci_consensus_mfe"], -6.0)
        self.assertEqual(metrics["sci_mean_single_mfe"], -4.0)
        self.assertEqual(metrics["sci_single_mfes"], [-4.0, -4.0])
        self.assertEqual(metrics["sci_vienna_version"], {
            "RNAalifold": "RNAalifold 2.4.18",
            "RNAfold": "RNAfold 2.4.18",
        })

    def test_sci_rejects_unexpected_vienna_version(self):
        prediction = parse_alignment(REPOSITORY / "tests/data/tiny_prediction.fa")
        with tempfile.TemporaryDirectory(prefix="dafs-sci-version-test-") as directory:
            tools = Path(directory)
            for name in ("RNAalifold", "RNAfold"):
                tool = tools / name
                tool.write_text(dedent(f"""\
                    #!{sys.executable}
                    print("{{name}} 2.4.18")
                """), encoding="utf-8")
                tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
            config = SCIConfig(tools / "RNAalifold", tools / "RNAfold", "2.5.0", 10)
            with self.assertRaisesRegex(RuntimeError, "version"):
                calculate_sci(prediction, config)


if __name__ == "__main__":
    unittest.main()
