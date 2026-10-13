#!/usr/bin/env python3

import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "benchmarks"))
from run_competitor import (adapt_prediction, resolve_environment,
                            resolve_executable, unsupported_input_symbols)
from score import SCIConfig, StructuralAlignment


class CompetitorRunnerTest(unittest.TestCase):
    def test_unsupported_symbols_are_reported_without_changing_input(self):
        with tempfile.TemporaryDirectory(prefix="competitor-alphabet-test-") as directory:
            fasta = Path(directory) / "input.fa"
            fasta.write_text(">seq1\nacguNRY-.\n", encoding="utf-8")

            self.assertEqual(
                unsupported_input_symbols(fasta, "ACGTU"), ["-", ".", "N", "R", "Y"])
            self.assertEqual(unsupported_input_symbols(fasta, None), [])

    def test_incompatible_input_is_recorded_as_nonfatal_skip(self):
        with tempfile.TemporaryDirectory(prefix="competitor-skip-test-") as directory:
            root = Path(directory)
            fasta = root / "ambiguous.fa"
            fasta.write_text(">seq1\nACGUN\n>seq2\nACGU\n", encoding="utf-8")
            executable = root / "unused-aligner"
            executable.write_text("#!/bin/sh\nexit 99\n", encoding="utf-8")
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            config_path = root / "config.json"
            config_path.write_text(json.dumps({
                "output_dir": str(root / "results"),
                "datasets": [{"id": "ambiguous", "input": str(fasta)}],
                "methods": [{
                    "id": "guarded-aligner",
                    "binary": str(executable),
                    "argv": ["{binary}", "{input}"],
                    "input_alphabet": "ACGTU",
                }],
            }), encoding="utf-8")

            completed = subprocess.run(
                [sys.executable, str(REPOSITORY / "benchmarks/run_competitor.py"),
                 "--config", str(config_path)],
                cwd=REPOSITORY, capture_output=True, text=True, check=False)

            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            result_path = (root / "results/runs/guarded-aligner/ambiguous/rep-01"
                           / "result.json")
            result = json.loads(result_path.read_text(encoding="utf-8"))
            self.assertEqual(result["status"], "skipped_incompatible_input")
            self.assertEqual(result["metrics_validation"]["unsupported_symbols"], ["N"])
            self.assertFalse((result_path.parent / "prediction.fa").exists())

    def test_path_environment_is_relative_to_config_file(self):
        with tempfile.TemporaryDirectory(prefix="competitor-env-test-") as directory:
            root = Path(directory)
            config_dir = root / "configs"
            config_dir.mkdir()
            expected = root / "tools" / "mafft"
            expected.mkdir(parents=True)

            environment = resolve_environment(
                {"environment": {"MAFFT_BINARIES": "../tools/mafft"}},
                config_dir,
            )

            self.assertEqual(environment["MAFFT_BINARIES"], str(expected))

    def test_sci_vienna_bin_is_added_to_child_path(self):
        with tempfile.TemporaryDirectory(prefix="competitor-vienna-path-") as directory:
            vienna_bin = Path(directory) / "bin"
            vienna_bin.mkdir()
            sci = SCIConfig(vienna_bin / "RNAalifold", vienna_bin / "RNAfold")

            environment = resolve_environment({}, Path(directory), sci)

            self.assertEqual(environment["PATH"].split(os.pathsep, 1)[0],
                             str(vienna_bin))

    def test_output_adapter_drops_tool_consensus_record(self):
        prediction = StructuralAlignment(
            OrderedDict([
                ("seq1", "AC-G"),
                ("consensus", "(..)"),
                ("seq2", "ACUG"),
            ]),
            None,
        )

        adapted = adapt_prediction(
            {"drop_sequence_names": ["consensus"]}, prediction)

        self.assertEqual(list(adapted.sequences), ["seq1", "seq2"])
        self.assertEqual(adapted.sequences["seq1"], "AC-G")

    def test_output_adapter_requires_expected_record(self):
        prediction = StructuralAlignment(OrderedDict([("seq1", "ACG")]), None)

        with self.assertRaisesRegex(ValueError, "not present"):
            adapt_prediction({"drop_sequence_names": ["consensus"]}, prediction)

    def test_resolve_executable_preserves_method_symlink(self):
        with tempfile.TemporaryDirectory(prefix="competitor-symlink-test-") as directory:
            root = Path(directory)
            target = root / "mafft"
            target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            target.chmod(target.stat().st_mode | stat.S_IXUSR)
            alias = root / "mafft-qinsi"
            alias.symlink_to(target.name)

            self.assertEqual(resolve_executable(str(alias), root), alias)

    def test_smoke_harness_normalizes_external_outputs(self):
        with tempfile.TemporaryDirectory(prefix="competitor-runner-test-") as directory:
            root = Path(directory)
            tools = root / "tools"
            tools.mkdir()

            aligner = tools / "fake-aligner"
            aligner.write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path

records = []
name = None
for raw in Path(sys.argv[-1]).read_text().splitlines():
    if raw.startswith('>'):
        name = raw[1:].strip()
        records.append([name, ''])
    elif name is not None:
        records[-1][1] += raw.strip()
width = max(len(sequence) for _, sequence in records)
for name, sequence in records:
    print(f'>{name}')
    print(sequence.ljust(width, '-'))
""", encoding="utf-8")
            aligner.chmod(aligner.stat().st_mode | stat.S_IXUSR)

            locarna = tools / "fake-mlocarna"
            locarna.write_text("""#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

target = Path(sys.argv[sys.argv.index('--tgtdir') + 1])
target.joinpath('results').mkdir(parents=True, exist_ok=True)
source = Path(sys.argv[-1])
output = subprocess.check_output([sys.executable, str(Path(__file__).with_name('fake-aligner')), str(source)], text=True)
records = []
name = None
for raw in output.splitlines():
    if raw.startswith('>'):
        name = raw[1:]
        records.append([name, ''])
    elif name is not None:
        records[-1][1] += raw.strip()
with target.joinpath('results', 'result.aln').open('w') as stream:
    stream.write('CLUSTAL W (1.82) multiple sequence alignment\\n\\n')
    for name, sequence in records:
        stream.write(f'{name:<12} {sequence}\\n')
""", encoding="utf-8")
            locarna.chmod(locarna.stat().st_mode | stat.S_IXUSR)

            rnaalifold = tools / "RNAalifold"
            rnaalifold.write_text("""#!/usr/bin/env python3
import sys

if '--version' in sys.argv or '-V' in sys.argv:
    print('RNAalifold 2.7.2')
    raise SystemExit(0)
text = sys.stdin.read().splitlines()
sequences = [line.strip() for line in text if line.strip() and not line.startswith('>')]
width = len(sequences[0])
structure = '(' + '.' * (width - 2) + ')'
if '--sci' in sys.argv:
    print(f'{structure} (-6.00 = -5.00 + -1.00)')
else:
    print(f'{structure} (-6.00)')
""", encoding="utf-8")
            rnaalifold.chmod(rnaalifold.stat().st_mode | stat.S_IXUSR)

            rnafold = tools / "RNAfold"
            rnafold.write_text("""#!/usr/bin/env python3
import sys

if '--version' in sys.argv or '-V' in sys.argv:
    print('RNAfold 2.7.2')
    raise SystemExit(0)
records = sum(line.startswith('>') for line in sys.stdin)
for index in range(records):
    print(f'>seq{index + 1}')
    print('GGGAAACCC')
    print('(((...))) (-4.00)')
""", encoding="utf-8")
            rnafold.chmod(rnafold.stat().st_mode | stat.S_IXUSR)

            config = {
                "method_group": "test",
                "output_dir": str(root / "results"),
                "timeout_seconds": 60,
                "repetitions": 1,
                "dataset_manifests": [],
                "datasets": [{
                    "id": "tiny",
                    "input": str(REPOSITORY / "tests/data/tiny.fa"),
                    "reference": str(REPOSITORY / "tests/data/tiny_reference.sto"),
                }],
                "sci": {"bin_dir": str(tools), "version": "2.7.2"},
                "methods": [
                    {
                        "id": "locarna",
                        "binary": str(locarna),
                        "argv": ["{binary}", "--quiet", "--keep-sequence-order",
                                 "--tgtdir", "{run_dir}/tool", "{input}"],
                        "output_mode": "path",
                        "output_path": "{run_dir}/tool/results/result.aln",
                        "consensus_structure": "rnaalifold",
                    },
                    *({
                        "id": name,
                        "binary": str(aligner),
                        "argv": ["{binary}", "{input}"],
                        "output_mode": "stdout",
                        "consensus_structure": "rnaalifold",
                    } for name in ("raf", "centroidalign", "mafft-qinsi", "ginsi")),
                ],
            }
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(REPOSITORY / "benchmarks/run_competitor.py"),
                 "--config", str(config_path)],
                cwd=REPOSITORY, capture_output=True, text=True, check=False)

            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            summary = (root / "results" / "summary.jsonl").read_text().splitlines()
            self.assertEqual(len(summary), 5)
            self.assertTrue(all(json.loads(row)["status"] == "success" for row in summary))


if __name__ == "__main__":
    unittest.main()
