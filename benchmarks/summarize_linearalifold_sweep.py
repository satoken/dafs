#!/usr/bin/env python3
"""Summarize the fixed-budget LinearAlifold condition sweep."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

try:
    # Running this file directly puts benchmarks/ on sys.path.
    from run import command_digest, flatten_summary, sha256_file
except ImportError:  # pragma: no cover - import path used by package callers
    from benchmarks.run import command_digest, flatten_summary, sha256_file


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json"
SUMMARY = ROOT / "benchmarks/results-murlet-dataset1-linearalifold-repair-v2-iter100-beam100/summary.csv"
OUTPUT = ROOT / "docs/murlet_dataset1_linearalifold_repair_v2_iter100_beam100_ja.md"
METRICS = ("mcc", "sps", "alignment_ppv", "sensitivity", "ppv", "cbp_f1",
           "sci", "wall_seconds", "max_rss_kb")
BEAMS = ("--align-beam=100", "--align-dd-beam=100",
         "--linfold-beam=100", "--fold-dd-beam=100",
         "--fold-final-beam=100")
DIGEST_SCHEMA_VERSION = 2
EXPECTED_CONDITIONS = 480
EXPECTED_LINEAR_CONDITIONS = 320
EXPECTED_LINEAR_ALIFOLD_CONDITIONS = 160
PILOT_EXPECTED_CONDITIONS = 8
PILOT_EXPECTED_LINEAR_CONDITIONS = 8
PILOT_EXPECTED_LINEAR_ALIFOLD_CONDITIONS = 4
PILOT_EXPECTED_DATASETS = 3


def aggregate(rows: list[dict[str, str]]) -> dict[str, float]:
    return {field: mean(float(row[field]) for row in rows) for field in METRICS}


def rank(row: dict[str, float]) -> tuple[float, float, float]:
    return row["mcc"], row["sci"], row["sps"]


def csv_value_matches(actual: str, expected: object) -> bool:
    """Compare one CSV cell with the value emitted by ``flatten_summary``."""
    if expected is None:
        return actual == ""
    if isinstance(expected, bool):
        return actual == ("True" if expected else "False")
    if isinstance(expected, (int, float)):
        try:
            value = float(actual)
        except (TypeError, ValueError):
            return False
        return math.isfinite(value) and value == expected
    return actual == str(expected)


def validate_summary_row(row: dict[str, str], result: dict[str, object],
                         result_path: Path) -> None:
    """Ensure aggregation consumes exactly the values stored in result.json."""
    expected = flatten_summary(result)
    if set(row) != set(expected):
        missing = sorted(set(expected) - set(row), key=str)
        extra = sorted(set(row) - set(expected), key=str)
        raise ValueError(
            f"summary columns differ from result schema at {result_path}: "
            f"missing={missing}, extra={extra}")
    for field, expected_value in expected.items():
        if not csv_value_matches(row[field], expected_value):
            raise ValueError(
                f"summary row mismatch at {result_path}, field {field}: "
                f"CSV={row[field]!r}, result={expected_value!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--datasets", type=int,
                        help="dataset count (85 for full, 3 for --pilot)")
    parser.add_argument(
        "--pilot", action="store_true",
        help="validate the fixed 8-condition, 3-dataset pilot")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    metadata = {condition["id"]: condition for condition in config["conditions"]}
    if len(metadata) != len(config["conditions"]):
        raise ValueError("duplicate condition ids")
    expected_conditions = (PILOT_EXPECTED_CONDITIONS if args.pilot
                           else EXPECTED_CONDITIONS)
    expected_linear_conditions = (PILOT_EXPECTED_LINEAR_CONDITIONS
                                  if args.pilot else EXPECTED_LINEAR_CONDITIONS)
    expected_linear_alifold_conditions = (
        PILOT_EXPECTED_LINEAR_ALIFOLD_CONDITIONS if args.pilot
        else EXPECTED_LINEAR_ALIFOLD_CONDITIONS)
    expected_datasets = (PILOT_EXPECTED_DATASETS if args.pilot
                         else (args.datasets if args.datasets is not None else 85))
    if args.pilot and args.datasets not in (None, PILOT_EXPECTED_DATASETS):
        raise ValueError("--pilot requires exactly 3 datasets")
    if expected_datasets <= 0:
        raise ValueError("dataset count must be positive")
    if len(metadata) != expected_conditions:
        raise ValueError(
            f"expected {expected_conditions} conditions, found {len(metadata)}")
    if sum(bool(condition["linear"]) for condition in metadata.values()) != expected_linear_conditions:
        raise ValueError(
            f"expected {expected_linear_conditions} linear conditions")
    if sum(bool(condition["linear"]) and condition["id"].endswith("_alifold")
           for condition in metadata.values()) != expected_linear_alifold_conditions:
        raise ValueError(
            f"expected {expected_linear_alifold_conditions} LinearAlifold conditions")
    for cid, condition in metadata.items():
        options = condition["args"]
        if options.count("--max-iter=100") != 1:
            raise ValueError(f"max-iter is not fixed at 100: {cid}")
        if sum(options.count(flag) for flag in ("--alifold", "--no-alifold")) != 1:
            raise ValueError(f"alifold mode is not specified exactly once: {cid}")
        if condition["linear"]:
            if any(options.count(beam) != 1 for beam in BEAMS):
                raise ValueError(f"linear beams are not fixed at 100: {cid}")
        elif any(beam in options for beam in BEAMS):
            raise ValueError(f"nonlinear condition has linear beam flags: {cid}")
        if cid.endswith("_alifold"):
            base_id = cid.removesuffix("_alifold")
            if base_id not in metadata:
                raise ValueError(f"missing control pair for {cid}")
            base = metadata[base_id]
            if not condition["linear"] or not base["linear"] or \
                    condition["args"] != ["--alifold" if arg == "--no-alifold" else arg
                                          for arg in base["args"]]:
                raise ValueError(f"invalid alifold control pair: {cid}")
    with args.summary.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    expected = len(metadata) * expected_datasets
    if len(rows) != expected:
        raise ValueError(f"expected {expected} runs, found {len(rows)}")
    if {row["condition"] for row in rows} != set(metadata):
        raise ValueError("summary conditions do not match the configuration")
    if {row["repetition"] for row in rows} != {"1"}:
        raise ValueError("the fixed sweep must contain exactly repetition 1")
    if len({(row["condition"], row["dataset"], row["repetition"])
            for row in rows}) != expected:
        raise ValueError("duplicate condition/dataset/repetition rows")
    bad = [row for row in rows if row["status"] != "success" or
           row["metrics_valid"] != "True" or row["quality_complete"] != "True"]
    if bad:
        raise ValueError(f"{len(bad)} runs failed integrity checks")

    manifest_path = args.summary.parent / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("schema_version") != 2 or
            manifest.get("digest_schema_version") != DIGEST_SCHEMA_VERSION):
        raise ValueError("benchmark manifest uses an unsupported digest schema")
    file_hash_cache: dict[Path, str] = {}

    def current_sha256(value: object, label: str, result_path: Path) -> str:
        """Hash a provenance file once, keyed by its canonical path."""
        if not isinstance(value, str) or not value:
            raise ValueError(f"missing {label} path: {result_path}")
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"{label} file is missing at {result_path}: {path}")
        if path not in file_hash_cache:
            file_hash_cache[path] = sha256_file(path)
        return file_hash_cache[path]

    binary = Path(manifest["binary"])
    if (not binary.is_file() or
            current_sha256(str(binary), "benchmark binary", manifest_path) !=
            manifest["binary_sha256"]):
        raise ValueError("benchmark binary differs from the manifest")
    time_binary = Path(manifest["time_binary"])
    if (not time_binary.is_file() or
            current_sha256(str(time_binary), "GNU time binary", manifest_path) !=
            manifest["time_binary_sha256"]):
        raise ValueError("GNU time binary differs from the manifest")
    if manifest.get("sci"):
        for key, label in (("rnaalifold", "RNAalifold"),
                           ("rnafold", "RNAfold")):
            sci_binary = Path(manifest["sci"][key])
            sha_key = f"{key}_sha256"
            if (not sci_binary.is_file() or
                    current_sha256(str(sci_binary), f"SCI {label} binary",
                                   manifest_path) !=
                    manifest["sci"][sha_key]):
                raise ValueError(f"SCI {label} binary differs from the manifest")
    if (not args.config.is_file() or
            current_sha256(str(args.config), "benchmark config", manifest_path) !=
            manifest["config_sha256"]):
        raise ValueError("benchmark config differs from the manifest")
    for row in rows:
        result_path = (args.summary.parent / "runs" / row["condition"] /
                       row["dataset"] / f"rep-{int(row['repetition']):02d}" /
                       "result.json")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("schema_version") != 2:
            raise ValueError(f"unsupported result schema: {result_path}")
        if (result.get("condition"), result.get("dataset"),
                str(result.get("repetition"))) != (
                    row["condition"], row["dataset"], row["repetition"]):
            raise ValueError(f"result identity mismatch: {result_path}")
        validate_summary_row(row, result, result_path)
        if (result.get("started_utc", "") < manifest["created_utc"] or
                result.get("command", [None])[0] != manifest["binary"]):
            raise ValueError(f"stale or mismatched run: {result_path}")
        if (result.get("binary_sha256") != manifest["binary_sha256"] or
                result.get("config_sha256") != manifest["config_sha256"]):
            raise ValueError(f"result provenance mismatch: {result_path}")
        context = result.get("digest_context")
        if not isinstance(context, dict) or context.get("schema_version") != DIGEST_SCHEMA_VERSION:
            raise ValueError(f"missing digest context: {result_path}")
        manifest_context = manifest.get("digest_context", {})
        for key in ("binary", "binary_sha256", "config_path", "config_sha256",
                    "time_binary", "time_binary_sha256", "environment", "sci"):
            if context.get(key) != manifest_context.get(key):
                raise ValueError(f"result run context mismatch: {result_path}")
        if (context.get("binary_sha256") != manifest["binary_sha256"] or
                context.get("config_sha256") != manifest["config_sha256"] or
                context.get("condition") != row["condition"] or
                context.get("dataset") != row["dataset"] or
                context.get("repetition") != int(row["repetition"])):
            raise ValueError(f"result digest context mismatch: {result_path}")
        input_sha256 = current_sha256(result.get("input"), "input", result_path)
        if input_sha256 != result.get("input_sha256"):
            raise ValueError(f"input checksum mismatch: {result_path}")
        reference = result.get("reference")
        if reference is None:
            reference_sha256 = None
        else:
            reference_sha256 = current_sha256(reference, "reference", result_path)
        if (reference_sha256 != result.get("reference_sha256") or
                context.get("reference_sha256") != reference_sha256):
            raise ValueError(f"reference checksum mismatch: {result_path}")
        command = result.get("command")
        artifacts = result.get("artifacts", {})
        expected_command = [
            manifest["binary"], *map(str, metadata[row["condition"]]["args"]),
            "--metrics-jsonl", str(artifacts.get("metrics")),
            str(result.get("input")),
        ]
        if command != expected_command:
            raise ValueError(f"result command mismatch: {result_path}")
        if result.get("command_sha256") != command_digest(
                command, input_sha256, context):
            raise ValueError(f"result digest mismatch: {result_path}")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_case: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        grouped[row["condition"]].append(row)
        by_case[row["condition"], row["dataset"]] = row
    dataset_ids = {row["dataset"] for row in rows}
    if len(dataset_ids) != expected_datasets or any(
            {row["dataset"] for row in grouped[cid]} != dataset_ids
            for cid in metadata):
        raise ValueError("one or more conditions have incomplete dataset coverage")
    summaries = {cid: aggregate(grouped[cid]) for cid in metadata}

    lines = [
        "# Murlet dataset1 LinearAlifold 条件スイープ",
        "",
        f"{expected_datasets}ケース、{len(metadata)}条件、{len(rows)} run。",
        "全条件で `--max-iter=100`、linear系の5種類のbeamを100に固定した。",
        f"順位は{expected_datasets}ケースのmacro平均MCC、SCI、SPSの辞書式順序で決めた。",
        "全runで `status=success`、`metrics_valid=True`、`quality_complete=True` を確認した。",
        "全runの開始時刻と実行ファイルのパスをmanifestと照合し、バイナリと設定のSHA256も確認した。",
        f"評価バイナリ SHA256: `{manifest['binary_sha256']}`。",
        f"設定 SHA256: `{manifest['config_sha256']}`。実行開始: `{manifest['created_utc']}`。",
        "同じケースを条件選択と結果報告に使う探索結果であり、独立したテスト集合の精度ではない。",
        "",
        "## 構成別の最良条件",
        "",
        "| 構成 | profile | 条件 | MCC | SPS | CBP F1 | SCI | wall秒 | RSS MiB |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    families: dict[tuple[str, str], list[str]] = defaultdict(list)
    for cid, item in metadata.items():
        label = str(item["label"]).removesuffix(" (LinearAlifold)")
        profile = "LinearAlifold" if cid.endswith("_alifold") else (
            "RNAalifold" if not item["linear"] and "(Alifold)" in label else
            "off" if "--no-alifold" in item["args"] else "on")
        families[label, profile].append(cid)
    winners = {}
    for key, ids in sorted(families.items()):
        cid = max(ids, key=lambda name: rank(summaries[name]))
        winners[key] = cid
        stats = summaries[cid]
        lines.append("| {} | {} | `{}` | {:.4f} | {:.4f} | {:.4f} | {:.4f} | {:.2f} | {:.1f} |".format(
            key[0], key[1], cid, stats["mcc"], stats["sps"],
            stats["cbp_f1"], stats["sci"], stats["wall_seconds"],
            stats["max_rss_kb"] / 1024))

    phases = ("base", "parameter", "dd_ribosum", "final_ribosum")
    lines += [
        "",
        "## phase 別の MCC 最良条件",
        "",
        f"各セルは条件 ID と{expected_datasets}ケース平均MCC。",
        "",
        "| 構成 | profile | base | parameter | DD RIBOSUM | 最終 RIBOSUM |",
        "|---|---|---|---|---|---|",
    ]
    for key, ids in sorted(families.items()):
        cells = []
        for phase in phases:
            candidates = [cid for cid in ids if metadata[cid]["phase"] == phase]
            if not candidates:
                cells.append("—")
                continue
            cid = max(candidates, key=lambda name: rank(summaries[name]))
            cells.append(f"`{cid}` ({summaries[cid]['mcc']:.4f})")
        lines.append(f"| {key[0]} | {key[1]} | " + " | ".join(cells) + " |")

    paired = []
    paired_by_architecture: dict[str, list[dict[str, float]]] = defaultdict(list)
    for cid in metadata:
        if not cid.endswith("_alifold"):
            continue
        base = cid.removesuffix("_alifold")
        for dataset in {row["dataset"] for row in grouped[cid]}:
            on = by_case[cid, dataset]
            off = by_case[base, dataset]
            delta = {field: float(on[field]) - float(off[field])
                     for field in ("mcc", "sps", "sci", "wall_seconds", "max_rss_kb")}
            paired.append(delta)
            paired_by_architecture[metadata[base]["label"]].append(delta)
    lines += [
        "",
        "## LinearAlifold の対応条件差",
        "",
        f"linear系の同一条件・同一ケースを対応させた{len(paired)}組の記述統計。",
        "条件間には共通ケースがあるため、独立標本としての有意性検定は行わない。",
        "非linear系のAlifold有効・無効群は他の引数も異なるため、Alifold単独の効果として扱わない。",
        "",
        "| 指標 | 有効−無効の平均差 | 有効が高い組数 |",
        "|---|---:|---:|",
    ]
    for field in ("mcc", "sps", "sci", "wall_seconds", "max_rss_kb"):
        delta = [row[field] for row in paired]
        scale = 1024 if field == "max_rss_kb" else 1
        label = "RSS MiB" if field == "max_rss_kb" else field
        lines.append(f"| {label} | {mean(delta) / scale:+.4f} | {sum(value > 0 for value in delta)}/{len(delta)} |")
    lines += [
        "",
        "### linear 構成ごとの対応条件差",
        "",
        f"各構成について条件数×{expected_datasets}ケースの対応差を集計した。条件とケースの反復は相関している。",
        "",
        "| 構成 | 組数 | ΔMCC | MCC改善組 | ΔSPS | ΔSCI | Δwall秒 | ΔRSS MiB |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for architecture, deltas in sorted(paired_by_architecture.items()):
        lines.append("| {} | {} | {:+.4f} | {}/{} | {:+.4f} | {:+.4f} | {:+.2f} | {:+.1f} |".format(
            architecture, len(deltas), mean(x["mcc"] for x in deltas),
            sum(x["mcc"] > 0 for x in deltas), len(deltas),
            mean(x["sps"] for x in deltas), mean(x["sci"] for x in deltas),
            mean(x["wall_seconds"] for x in deltas),
            mean(x["max_rss_kb"] for x in deltas) / 1024))

    lines += [
        "",
        "## 全条件の上位10件",
        "",
        "| 順位 | 条件 | MCC | SPS | SCI | wall秒 |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for index, cid in enumerate(sorted(metadata, key=lambda name: rank(summaries[name]),
                                       reverse=True)[:10], 1):
        stats = summaries[cid]
        lines.append("| {} | `{}` | {:.4f} | {:.4f} | {:.4f} | {:.2f} |".format(
            index, cid, stats["mcc"], stats["sps"], stats["sci"],
            stats["wall_seconds"]))

    lines += ["", "## 最良条件の全引数", ""]
    for key, cid in sorted(winners.items()):
        lines += [f"### {key[0]} — {key[1]} — `{cid}`", "", "```text",
                  " ".join(metadata[cid]["args"]), "```", ""]
    lines += ["## 成果物", "", f"- 設定: `{args.config}`",
              f"- 実測CSV: `{args.summary}`", ""]
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
