#!/usr/bin/env python3
"""Regression tests for resumable benchmark provenance."""

from __future__ import annotations

import json
import csv
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
RUNNER = REPOSITORY / "benchmarks" / "run.py"
SUMMARIZER = REPOSITORY / "benchmarks" / "summarize_linearalifold_sweep.py"
sys.path.insert(0, str(REPOSITORY / "benchmarks"))
from run import command_digest, flatten_summary, sha256_file


def mock_binary_source(marker: str) -> str:
    """Return a tiny executable that emits valid DAFS runner artifacts."""
    return f"""#!{sys.executable}
import json
import pathlib
import sys

marker = {marker!r}
metrics = pathlib.Path(sys.argv[sys.argv.index('--metrics-jsonl') + 1])
metrics.write_text(json.dumps({{'event': 'run_summary', 'seconds': 0.0, 'marker': marker}}) + '\\n')
input_path = pathlib.Path(sys.argv[-1])
sys.stdout.write(input_path.read_text())
"""


def make_summarizer_fixture(root: Path) -> tuple[Path, Path, Path, Path]:
    """Build the smallest complete pilot accepted by the sweep summarizer."""
    binary = Path(sys.executable).resolve()
    input_path = root / "input.fa"
    reference_path = root / "reference.sto"
    input_path.write_text("input\n", encoding="utf-8")
    reference_path.write_text("reference\n", encoding="utf-8")
    config_path = root / "config.json"
    output_dir = root / "results"

    conditions = []
    phases = ("base", "parameter", "dd_ribosum", "final_ribosum")
    common_args = [
        "--max-iter=100", "--align-beam=100", "--align-dd-beam=100",
        "--linfold-beam=100", "--fold-dd-beam=100",
        "--fold-final-beam=100",
    ]
    for index, phase in enumerate(phases):
        base_id = f"arch{index}"
        base_args = ["--no-alifold", *common_args]
        conditions.append({
            "id": base_id, "label": base_id, "linear": True,
            "phase": phase, "args": base_args,
        })
        conditions.append({
            "id": f"{base_id}_alifold",
            "label": f"{base_id} (LinearAlifold)", "linear": True,
            "phase": phase,
            "args": ["--alifold", *common_args],
        })
    config = {
        "binary": str(binary), "output_dir": str(output_dir),
        "conditions": conditions,
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    binary_sha256 = sha256_file(binary)
    config_sha256 = sha256_file(config_path)
    time_binary = str(binary)
    time_binary_sha256 = binary_sha256
    environment = {"OMP_NUM_THREADS": "1"}
    context_base = {
        "schema_version": 2, "binary": str(binary),
        "binary_sha256": binary_sha256,
        "config_path": str(config_path.resolve()),
        "config_sha256": config_sha256,
        "time_binary": time_binary,
        "time_binary_sha256": time_binary_sha256,
        "environment": environment, "sci": None,
    }
    manifest = {
        "schema_version": 2, "digest_schema_version": 2,
        "created_utc": "2020-01-01T00:00:00+00:00",
        "config": config, "config_path": str(config_path.resolve()),
        "config_sha256": config_sha256, "binary": str(binary),
        "binary_sha256": binary_sha256, "time_binary": time_binary,
        "time_binary_sha256": time_binary_sha256,
        "environment": environment, "sci": None,
        "digest_context": context_base,
    }
    output_dir.mkdir()
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    input_sha256 = sha256_file(input_path)
    reference_sha256 = sha256_file(reference_path)
    rows = []
    for condition_index, condition in enumerate(conditions):
        for dataset_index in range(3):
            dataset = f"dataset{dataset_index}"
            run_dir = output_dir / "runs" / condition["id"] / dataset / "rep-01"
            run_dir.mkdir(parents=True)
            metrics_path = run_dir / "metrics.jsonl"
            metrics_path.write_text("", encoding="utf-8")
            command = [str(binary), *condition["args"],
                       "--metrics-jsonl", str(metrics_path), str(input_path)]
            context = {
                **context_base, "condition": condition["id"],
                "dataset": dataset, "repetition": 1,
                "reference_sha256": reference_sha256,
            }
            result = {
                "schema_version": 2, "method": "DAFS",
                "condition": condition["id"], "ribosum": 0.0,
                "dataset": dataset, "repetition": 1, "status": "success",
                "return_code": 0, "timed_out": False,
                "started_utc": "2020-01-01T00:00:01+00:00",
                "harness_wall_seconds": 1.0, "command": command,
                "command_sha256": command_digest(command, input_sha256, context),
                "digest_context": context, "binary": str(binary),
                "binary_sha256": binary_sha256,
                "config_path": str(config_path.resolve()),
                "config_sha256": config_sha256, "input": str(input_path),
                "input_sha256": input_sha256, "reference": str(reference_path),
                "reference_sha256": reference_sha256,
                "resource": {"wall_seconds": 1.0, "user_seconds": 0.1,
                             "system_seconds": 0.01, "max_rss_kb": 1000},
                "metrics_validation": {"valid": True, "run_summary": {}},
                "quality": {
                    "sps": 0.5, "sensitivity": 0.5,
                    "alignment_ppv": 0.5, "ppv": 0.5,
                    "mcc": 0.4 + 0.01 * condition_index,
                    "cbp_f1": 0.5, "sci": 0.7,
                },
                "quality_complete": True, "score_error": None,
                "artifacts": {
                    "prediction": str(run_dir / "prediction.fa"),
                    "stderr": str(run_dir / "stderr.log"),
                    "metrics": str(metrics_path),
                    "resource": str(run_dir / "resource.json"),
                },
            }
            (run_dir / "result.json").write_text(
                json.dumps(result, indent=2) + "\n", encoding="utf-8")
            rows.append(flatten_summary(result))
    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return config_path, summary_path, input_path, reference_path


class BenchmarkRunnerProvenanceTest(unittest.TestCase):
    def run_benchmark(self, config: Path, *, check: bool = True) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(RUNNER), "--config", str(config)],
            cwd=REPOSITORY,
            check=check,
            text=True,
            capture_output=True,
        )

    def run_summarizer(self, config: Path, summary: Path,
                       output: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SUMMARIZER), "--pilot",
             "--config", str(config), "--summary", str(summary),
             "--output", str(output)],
            cwd=REPOSITORY, text=True, capture_output=True,
        )

    def test_summarizer_rejects_tampered_csv_and_accepts_unchanged_pilot(self) -> None:
        with tempfile.TemporaryDirectory(prefix="dafs-summarizer-csv-") as directory:
            root = Path(directory)
            config, summary, _, _ = make_summarizer_fixture(root)
            valid = self.run_summarizer(config, summary, root / "valid.md")
            self.assertEqual(valid.returncode, 0, valid.stderr)

            with summary.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            fieldnames = list(rows[0])
            rows[0]["mcc"] = "9.0"
            with summary.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)
            tampered = self.run_summarizer(config, summary, root / "tampered.md")
            self.assertNotEqual(tampered.returncode, 0)
            self.assertIn("summary row mismatch", tampered.stderr)

    def test_summarizer_rehashes_input_and_reference_files(self) -> None:
        with tempfile.TemporaryDirectory(prefix="dafs-summarizer-hash-") as directory:
            root = Path(directory)
            config, summary, input_path, reference_path = make_summarizer_fixture(root)
            valid = self.run_summarizer(config, summary, root / "valid.md")
            self.assertEqual(valid.returncode, 0, valid.stderr)

            input_path.write_text("changed input\n", encoding="utf-8")
            changed_input = self.run_summarizer(config, summary, root / "input.md")
            self.assertNotEqual(changed_input.returncode, 0)
            self.assertIn("input checksum mismatch", changed_input.stderr)

            # Restore the input so the reference check reaches the changed file.
            input_path.write_text("input\n", encoding="utf-8")
            reference_path.write_text("changed reference\n", encoding="utf-8")
            changed_reference = self.run_summarizer(
                config, summary, root / "reference.md")
            self.assertNotEqual(changed_reference.returncode, 0)
            self.assertIn("reference checksum mismatch", changed_reference.stderr)

    def test_resume_skips_same_binary_and_reruns_changed_binary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="dafs-runner-provenance-") as directory:
            root = Path(directory)
            binary = root / "mock-dafs"
            binary.write_text(mock_binary_source("first"), encoding="utf-8")
            binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
            input_path = root / "input.fa"
            input_path.write_text(">seq1\nAAAA\n", encoding="utf-8")
            output_dir = root / "results"
            config = root / "config.json"
            config.write_text(json.dumps({
                "binary": str(binary),
                "output_dir": str(output_dir),
                "timeout_seconds": 20,
                "repetitions": 1,
                "random_seed": 42,
                "environment": {"OMP_NUM_THREADS": "1"},
                "datasets": [{"id": "tiny", "input": str(input_path)}],
                "conditions": [{"id": "mock", "args": []}],
            }), encoding="utf-8")

            first = self.run_benchmark(config)
            self.assertIn("RUN", first.stdout)
            result_path = output_dir / "runs" / "mock" / "tiny" / "rep-01" / "result.json"
            manifest_path = output_dir / "manifest.json"
            first_result = json.loads(result_path.read_text(encoding="utf-8"))
            first_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(first_result["status"], "success")
            self.assertEqual(first_result["schema_version"], 2)
            self.assertEqual(first_result["binary_sha256"], first_manifest["binary_sha256"])
            with (output_dir / "summary.csv").open(newline="", encoding="utf-8") as stream:
                summary_rows = list(csv.DictReader(stream))
            self.assertEqual(len(summary_rows), 1)
            self.assertEqual(
                (summary_rows[0]["condition"], summary_rows[0]["dataset"]),
                ("mock", "tiny"),
            )
            first_started = first_result["started_utc"]
            first_manifest_created = first_manifest["created_utc"]

            second = self.run_benchmark(config)
            self.assertIn("SKIP mock/tiny/rep-01", second.stdout)
            second_result = json.loads(result_path.read_text(encoding="utf-8"))
            second_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(second_result["started_utc"], first_started)
            self.assertEqual(second_manifest["created_utc"], first_manifest_created)
            self.assertEqual(second_manifest["binary_sha256"], first_manifest["binary_sha256"])

            # Keep the path and config unchanged while changing executable bytes.
            binary.write_text(mock_binary_source("second"), encoding="utf-8")
            os.chmod(binary, binary.stat().st_mode | stat.S_IXUSR)
            third = self.run_benchmark(config)
            self.assertIn("RUN", third.stdout)
            self.assertNotIn("SKIP mock/tiny/rep-01", third.stdout)
            third_result = json.loads(result_path.read_text(encoding="utf-8"))
            third_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertNotEqual(third_result["binary_sha256"], first_result["binary_sha256"])
            self.assertEqual(third_result["binary_sha256"], third_manifest["binary_sha256"])
            self.assertNotEqual(third_result["started_utc"], first_started)
            self.assertNotEqual(third_manifest["created_utc"], first_manifest_created)

            # A failed expected job must remain visible and make the runner fail.
            binary.write_text(
                f"#!{sys.executable}\nraise SystemExit(7)\n",
                encoding="utf-8",
            )
            os.chmod(binary, binary.stat().st_mode | stat.S_IXUSR)
            failed = self.run_benchmark(config, check=False)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("status=failed", failed.stderr)
            failed_result = json.loads(result_path.read_text(encoding="utf-8"))
            self.assertEqual(failed_result["status"], "failed")


if __name__ == "__main__":
    unittest.main()
