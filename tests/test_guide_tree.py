"""Test guide-tree selection and preservation of the original sequence input."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = str(Path(sys.argv.pop(1)).resolve())
FASTA = ">first name,:\nGGGAAACCC\n>second\nGGGAAUCCC\n>third\nGGGACACCC\n"


class GuideTreeTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory(prefix="dafs tree test ")
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.input = self.root / "input with spaces.fa"
        self.input.write_text(FASTA)

    def run_dafs(self, *arguments, data=None):
        if data is not None:
            self.input.write_text(data)
        # Prove that every tree works without any executable on PATH.
        return subprocess.run(
            [BINARY, "-a", "CONTRAlign", "-s", "CONTRAfold", "--max-iter=2",
             *arguments, str(self.input)], env=dict(os.environ, PATH=""),
            text=True, capture_output=True, timeout=10)

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(">SS_cons\n", result.stdout)

    def test_default_is_mafft_linkage_with_dafs_distances(self):
        metrics = self.root / "metrics.jsonl"
        implicit = self.run_dafs("--metrics-jsonl=" + str(metrics))
        explicit = self.run_dafs("--guide-tree=mafft-dafs")
        self.assert_success(implicit)
        self.assert_success(explicit)
        self.assertEqual(implicit.stdout, explicit.stdout)
        records = [json.loads(line) for line in metrics.read_text().splitlines()]
        self.assertEqual(records[0]["guide_tree_method"], "mafft-dafs")
        stage = next(item for item in records if item.get("name") == "guide_tree")
        self.assertEqual(stage["method"], "mafft-dafs")
        self.assertGreater(stage["seconds"], 0)

    def test_all_trees_on_linear_and_nonlinear_paths(self):
        for method in ["mafft-dafs", "mafft", "mafft-linear", "dafs", "raf"]:
            for model in [(), ("--alifold",),
                          ("-a", "LinearAlign", "-s", "lpc", "--dynamic-cbp")]:
                with self.subTest(method=method, model=model):
                    result = self.run_dafs("--guide-tree=" + method, *model)
                    self.assert_success(result)
                    records = {}
                    for line in result.stdout.splitlines()[1:]:
                        if line.startswith(">"):
                            name = line[1:].strip()
                            records[name] = ""
                        else:
                            records[name] += line
                    self.assertEqual(records["first name,:"].replace("-", ""), "GGGAAACCC")
                    self.assertEqual(records["second"].replace("-", ""), "GGGAAUCCC")
                    self.assertEqual(records["third"].replace("-", ""), "GGGACACCC")
                    self.assertEqual(list(records), ["SS_cons", "first name,:", "second", "third"])
                    self.assertEqual(self.input.read_text(), FASTA)

    def test_dafs_distances_unchanged_by_linkage_and_isolated_option(self):
        for model in [(), ("-a", "LinearAlign-CONTRAlign", "-s", "lpc", "--dynamic-cbp")]:
            matrices = []
            for method in ["dafs", "mafft-dafs"]:
                for isolated in [(), ("--disallow-isolated",)]:
                    result = self.run_dafs("--guide-tree=" + method,
                                           "--guide-tree-diagnostics", *model, *isolated)
                    self.assert_success(result)
                    pairs = [line for line in result.stderr.splitlines()
                             if line.startswith("[guide-tree-pair] ")]
                    self.assertEqual(len(pairs), 3)
                    matrices.append(pairs)
            self.assertTrue(all(matrix == matrices[0] for matrix in matrices))

    def test_linear_contralign_keeps_five_probability_states(self):
        metrics = self.root / "five-state.jsonl"
        result = self.run_dafs("-a", "LinearAlign-CONTRAlign", "-s", "lpc",
                               "--raw-align-cutoff=0.01", "--align-th=0.001",
                               "--metrics-jsonl=" + str(metrics))
        self.assert_success(result)
        start = json.loads(metrics.read_text().splitlines()[0])
        self.assertEqual(start["alignment_probability_states"], 5)
        self.assertAlmostEqual(start["raw_alignment_cutoff"], 0.01)
        self.assertAlmostEqual(start["alignment_threshold"], 0.001)

    def test_raw_cutoff_default_and_explicit_value_match(self):
        implicit = self.run_dafs()
        explicit = self.run_dafs("--raw-align-cutoff=0.01")
        self.assert_success(implicit)
        self.assert_success(explicit)
        self.assertEqual(implicit.stdout, explicit.stdout)
        for cutoff in ["-2", "nan", "inf", "1.1"]:
            with self.subTest(cutoff=cutoff):
                self.assertNotEqual(self.run_dafs("--raw-align-cutoff=" + cutoff).returncode, 0)

    def test_raf_tree_uses_raw_posteriors_before_pct(self):
        default = self.run_dafs("--guide-tree=raf")
        changed = self.run_dafs("--guide-tree=raf", "--align-pct=0", "--fold-pct=0",
                                "--fourway-pct=0.1")
        self.assert_success(default)
        self.assert_success(changed)
        self.assertEqual(default.stdout.splitlines()[0], changed.stdout.splitlines()[0])

    def test_wide_beam_matches_exact_and_records_width(self):
        metrics = self.root / "beam.jsonl"
        exact = self.run_dafs("--guide-tree=mafft")
        beam = self.run_dafs("--guide-tree=mafft-linear", "--guide-tree-beam=20",
                             "--metrics-jsonl=" + str(metrics))
        self.assert_success(beam)
        self.assertEqual(exact.stdout, beam.stdout)
        records = [json.loads(line) for line in metrics.read_text().splitlines()]
        self.assertEqual(records[0]["guide_tree_method"], "mafft-linear")
        self.assertEqual(records[0]["guide_tree_beam"], 20)
        stage = next(item for item in records if item.get("name") == "guide_tree")
        self.assertEqual(stage["beam"], 20)

    def test_linear_default_beam_is_200(self):
        metrics = self.root / "linear-default.jsonl"
        result = self.run_dafs("--guide-tree=mafft-linear", "--metrics-jsonl=" + str(metrics))
        self.assert_success(result)
        records = [json.loads(line) for line in metrics.read_text().splitlines()]
        self.assertEqual(records[0]["guide_tree_beam"], 200)

    def test_duplicate_names_preserved(self):
        result = self.run_dafs(data=FASTA.replace("second", "third"))
        self.assert_success(result)
        self.assertEqual(result.stdout.count("> third\n"), 2)

    def test_single_sequence_on_all_trees(self):
        for method in ["mafft-dafs", "mafft", "mafft-linear", "dafs", "raf"]:
            with self.subTest(method=method):
                self.assert_success(self.run_dafs("--guide-tree=" + method,
                                                 data=">single\nGGGAAACCC\n"))

    def test_empty_input_and_invalid_method_fail(self):
        invalid = self.run_dafs("--guide-tree=typo")
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("--guide-tree", invalid.stderr)
        for width in ["0", "-1", "2147483648", "abc"]:
            with self.subTest(width=width):
                self.assertNotEqual(self.run_dafs("--guide-tree=mafft-linear",
                                                 "--guide-tree-beam=" + width).returncode, 0)
        for data in ["", ">empty\n", ">empty\n>other\nACGU\n"]:
            with self.subTest(data=data):
                self.assertNotEqual(self.run_dafs(data=data).returncode, 0)


if __name__ == "__main__":
    unittest.main()
