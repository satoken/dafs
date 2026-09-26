#!/usr/bin/env python3
"""Score a DAFS structural alignment against FASTA or Stockholm reference.

This dependency-free scorer reports SPS, secondary-structure sensitivity,
PPV and MCC, plus pair-pair (CBP) F1 for evaluating the RIBOSUM objective.
SCI can optionally be calculated with an explicitly selected ViennaRNA build.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


GAPS = frozenset("-.~_")
SCI_FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
MFE_IN_PARENS = re.compile(
    rf"\(\s*({SCI_FLOAT})\s*(?:=|\))")
MFE_LABEL = re.compile(
    rf"(?:minimum\s+free\s+energy|mfe)\s*=\s*({SCI_FLOAT})",
    re.IGNORECASE)


@dataclass(frozen=True)
class SCIConfig:
    """External ViennaRNA programs and reproducibility settings for SCI."""

    rnaalifold: str | Path
    rnafold: str | Path
    expected_version: str | None = None
    timeout_seconds: float = 3600.0

    def __post_init__(self) -> None:
        if self.timeout_seconds <= 0:
            raise ValueError("SCI timeout_seconds must be positive")


@dataclass(frozen=True)
class StructuralAlignment:
    sequences: OrderedDict[str, str]
    structure: str | None
    structures: OrderedDict[str, str] | None = None

    @property
    def columns(self) -> int:
        return len(next(iter(self.sequences.values())))


def _validate(sequences: OrderedDict[str, str], structure: str | None,
              source: Path,
              structures: OrderedDict[str, str] | None = None) -> StructuralAlignment:
    if not sequences:
        raise ValueError(f"no sequences found in {source}")
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"inconsistent alignment lengths in {source}: {sorted(lengths)}")
    columns = next(iter(lengths))
    if structure is not None and len(structure) != columns:
        raise ValueError(
            f"SS_cons length {len(structure)} != alignment length {columns} in {source}")
    if structures is not None:
        unknown = sorted(set(structures) - set(sequences))
        if unknown:
            raise ValueError(f"structures contain unknown sequences in {source}: {unknown}")
        for name, value in structures.items():
            if len(value) != columns:
                raise ValueError(
                    f"structure length for {name!r} is {len(value)}, "
                    f"not alignment length {columns} in {source}")
    return StructuralAlignment(sequences, structure, structures)


def parse_stockholm(path: Path) -> StructuralAlignment:
    fragments: OrderedDict[str, list[str]] = OrderedDict()
    structure: list[str] = []
    structure_fragments: OrderedDict[str, list[str]] = OrderedDict()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line == "//":
            continue
        if line.startswith("#=GC SS_cons"):
            parts = line.split(maxsplit=2)
            if len(parts) != 3:
                raise ValueError(f"malformed SS_cons line in {path}: {raw}")
            structure.append(parts[2])
        elif line.startswith("#=GR"):
            parts = line.split(maxsplit=3)
            if len(parts) == 4 and parts[2] == "SS":
                structure_fragments.setdefault(parts[1], []).append(parts[3])
        elif line.startswith("#"):
            continue
        else:
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"malformed Stockholm sequence line in {path}: {raw}")
            fragments.setdefault(parts[0], []).append(parts[1])
    sequences = OrderedDict((name, "".join(parts))
                            for name, parts in fragments.items())
    structures = OrderedDict(
        (name, "".join(parts)) for name, parts in structure_fragments.items())
    return _validate(sequences, "".join(structure) or None, path,
                     structures or None)


def parse_fasta_alignment(path: Path) -> StructuralAlignment:
    records: OrderedDict[str, list[str]] = OrderedDict()
    name: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            if not name:
                raise ValueError(f"empty FASTA header in {path}")
            if name in records:
                raise ValueError(f"duplicate FASTA record {name!r} in {path}")
            records[name] = []
        elif name is not None:
            records[name].append("".join(line.split()))
        else:
            # DAFS writes the guide tree before its FASTA-like prediction.
            continue
    structure_parts = records.pop("SS_cons", None)
    structure = "".join(structure_parts) if structure_parts is not None else None
    sequences = OrderedDict((key, "".join(parts)) for key, parts in records.items())
    return _validate(sequences, structure, path)


def parse_clustal_alignment(path: Path) -> StructuralAlignment:
    """Read the CLUSTAL-like alignment written by LocARNA and similar tools."""
    fragments: OrderedDict[str, list[str]] = OrderedDict()
    sequence_token = re.compile(r"[A-Za-z0-9_.:/|+\-~]+$")
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()
        if not line or line.upper().startswith(("CLUSTAL", "MUSCLE", "PROBCONS")):
            continue
        parts = line.split()
        if len(parts) < 2:
            # Consensus rows contain only symbols such as '*' or ':'.
            continue
        name, fragment = parts[0], parts[1]
        if not sequence_token.fullmatch(fragment):
            continue
        fragments.setdefault(name, []).append(fragment)
    sequences = OrderedDict((name, "".join(parts))
                            for name, parts in fragments.items())
    return _validate(sequences, None, path)


def parse_alignment(path: str | Path) -> StructuralAlignment:
    source = Path(path)
    first = next((line.strip() for line in source.read_text(encoding="utf-8").splitlines()
                  if line.strip()), "")
    if first.startswith("# STOCKHOLM"):
        return parse_stockholm(source)
    if first.upper().startswith(("CLUSTAL", "MUSCLE", "PROBCONS")):
        return parse_clustal_alignment(source)
    return parse_fasta_alignment(source)


def _linearturbofold_db_sort_key(path: Path) -> tuple[int, str]:
    prefix = path.name.split("_", 1)[0]
    try:
        return int(prefix), path.name
    except ValueError:
        return 10**9, path.name


def parse_linearturbofold(alignment_path: str | Path,
                          output_dir: str | Path) -> StructuralAlignment:
    """Read LinearTurboFold's output.aln and numbered per-sequence .db files."""
    source = Path(alignment_path)
    alignment = parse_fasta_alignment(source)
    db_paths = sorted(Path(output_dir).glob("*.db"),
                      key=_linearturbofold_db_sort_key)
    if len(db_paths) != len(alignment.sequences):
        raise ValueError(
            f"LinearTurboFold produced {len(db_paths)} .db files for "
            f"{len(alignment.sequences)} sequences in {source.parent}")

    aligned_structures: OrderedDict[str, str] = OrderedDict()
    for (name, aligned_sequence), db_path in zip(alignment.sequences.items(), db_paths):
        lines = [line.strip() for line in db_path.read_text(encoding="utf-8").splitlines()
                 if line.strip()]
        if len(lines) < 3:
            raise ValueError(f"malformed LinearTurboFold structure file: {db_path}")
        structure = lines[-1]
        residues = sum(symbol not in GAPS for symbol in aligned_sequence)
        if len(structure) != residues:
            raise ValueError(
                f"structure length {len(structure)} != {residues} residues "
                f"for {name!r} in {db_path}")
        index = 0
        projected: list[str] = []
        for symbol in aligned_sequence:
            if symbol in GAPS:
                projected.append("-")
            else:
                projected.append(structure[index])
                index += 1
        aligned_structures[name] = "".join(projected)
    return StructuralAlignment(alignment.sequences, None, aligned_structures)


def match_prediction_names(prediction: StructuralAlignment,
                           reference: StructuralAlignment) -> StructuralAlignment:
    """Match method output names to reference names when separators changed.

    LinearTurboFold replaces ``/`` in FASTA identifiers with ``_`` when it
    writes ``output.aln``.  Stockholm references retain the original names,
    so restore that one-to-one mapping before reference-based scoring.
    Exact matches take precedence; normalized matches must be unique.
    """
    reference_names = set(reference.sequences)
    normalized: dict[str, list[str]] = {}
    for name in reference_names:
        normalized.setdefault(name.replace("/", "_"), []).append(name)

    mapping: dict[str, str] = {}
    used: set[str] = set()
    for name in prediction.sequences:
        if name in reference_names:
            target = name
        else:
            candidates = normalized.get(name, [])
            if len(candidates) != 1:
                continue
            target = candidates[0]
        if target in used:
            raise ValueError(
                f"prediction sequence names map to duplicate reference {target!r}")
        mapping[name] = target
        used.add(target)

    if not mapping:
        return prediction
    sequences = OrderedDict(
        (mapping.get(name, name), sequence)
        for name, sequence in prediction.sequences.items())
    structures = None
    if prediction.structures is not None:
        structures = OrderedDict(
            (mapping.get(name, name), structure)
            for name, structure in prediction.structures.items())
    return StructuralAlignment(sequences, prediction.structure, structures)


def _resolve_executable(value: str | Path) -> str:
    executable = shutil.which(str(value))
    if executable is None:
        raise FileNotFoundError(
            f"ViennaRNA executable not found or not executable: {value}")
    return executable


def _run_sci_command(executable: str, arguments: list[str], input_text: str,
                     cwd: Path, timeout_seconds: float) -> tuple[str, str]:
    command = [executable, *arguments]
    try:
        completed = subprocess.run(
            command, input=input_text, text=True, capture_output=True,
            cwd=cwd, timeout=timeout_seconds, check=False)
    except subprocess.TimeoutExpired as error:
        raise TimeoutError(
            f"SCI command timed out after {timeout_seconds:g}s: "
            f"{' '.join(command)}") from error
    except OSError as error:
        raise RuntimeError(
            f"could not execute SCI command {' '.join(command)}: {error}") from error
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        if len(detail) > 1000:
            detail = detail[-1000:]
        raise RuntimeError(
            f"SCI command failed with exit code {completed.returncode}: "
            f"{' '.join(command)}\n{detail}")
    return completed.stdout, completed.stderr


def _tool_version(executable: str, timeout_seconds: float) -> str:
    """Return a concise version string from a ViennaRNA executable."""
    for flag in ("--version", "-V"):
        try:
            completed = subprocess.run(
                [executable, flag], input="", text=True, capture_output=True,
                timeout=timeout_seconds, check=False)
        except (OSError, subprocess.TimeoutExpired):
            continue
        output = "\n".join(
            part.strip() for part in (completed.stdout, completed.stderr)
            if part.strip())
        if completed.returncode == 0 and output:
            return output
    return "unknown"


def _mfe_from_line(line: str) -> float | None:
    labeled = MFE_LABEL.search(line)
    if labeled:
        return float(labeled.group(1))
    in_parentheses = MFE_IN_PARENS.search(line)
    return float(in_parentheses.group(1)) if in_parentheses else None


def _mfe_values(output: str, expected: int | None = None) -> list[float]:
    values = [value for line in output.splitlines()
              if (value := _mfe_from_line(line)) is not None]
    if expected is not None and len(values) != expected:
        tail = "\n".join(output.splitlines()[-12:])
        raise ValueError(
            f"expected {expected} MFE values from RNAfold, found {len(values)}\n"
            f"{tail}")
    return values


def _sci_fasta(alignment: StructuralAlignment, aligned: bool) -> str:
    records: list[str] = []
    for index, (name, sequence) in enumerate(alignment.sequences.items(), 1):
        if aligned:
            value = "".join("-" if symbol in GAPS else symbol
                            for symbol in sequence)
        else:
            value = "".join(symbol for symbol in sequence if symbol not in GAPS)
        if not value:
            raise ValueError(f"sequence {name!r} is empty after removing gaps")
        # Individual records use simple IDs; alignment IDs are retained because
        # RNAalifold may include them in diagnostics and output annotations.
        header = name.replace("\n", " ").strip() if aligned else f"seq{index}"
        if not header:
            header = f"seq{index}"
        records.extend((f">{header}", value))
    return "\n".join(records) + "\n"


def calculate_sci(alignment: StructuralAlignment,
                  config: SCIConfig) -> dict[str, object]:
    """Calculate the standard RNAalifold/RNAfold structure conservation index.

    SCI is the MFE of the aligned consensus divided by the mean of the
    individual-sequence MFEs.  ViennaRNA is deliberately invoked here rather
    than imported as a Python module so the exact external executables and
    their version can be recorded in benchmark results.
    """
    if not alignment.sequences:
        raise ValueError("cannot calculate SCI for an empty alignment")

    rnaalifold = _resolve_executable(config.rnaalifold)
    rnafold = _resolve_executable(config.rnafold)
    versions = {
        "RNAalifold": _tool_version(rnaalifold, config.timeout_seconds),
        "RNAfold": _tool_version(rnafold, config.timeout_seconds),
    }
    if config.expected_version:
        mismatches = {
            name: version for name, version in versions.items()
            if config.expected_version not in version
        }
        if mismatches:
            details = "; ".join(
                f"{name}={version!r}" for name, version in mismatches.items())
            raise RuntimeError(
                f"ViennaRNA version {config.expected_version!r} was requested, "
                f"but the executable reported {details}")

    with tempfile.TemporaryDirectory(prefix="dafs-sci-") as directory:
        workdir = Path(directory)
        consensus_output, _ = _run_sci_command(
            rnaalifold, ["--sci", "--noPS", "--input-format=F"],
            _sci_fasta(alignment, aligned=True), workdir,
            config.timeout_seconds)
        consensus_values = _mfe_values(consensus_output)
        if not consensus_values:
            tail = "\n".join(consensus_output.splitlines()[-12:])
            raise ValueError(f"could not parse RNAalifold consensus MFE\n{tail}")
        consensus_mfe = consensus_values[0]

        individual_output, _ = _run_sci_command(
            rnafold, ["--noPS"],
            _sci_fasta(alignment, aligned=False), workdir,
            config.timeout_seconds)
        individual_mfes = _mfe_values(
            individual_output, expected=len(alignment.sequences))

    mean_single_mfe = sum(individual_mfes) / len(individual_mfes)
    # This agrees with ViennaRNA's explicit zero convention for a zero
    # denominator and avoids manufacturing an infinity/NaN in JSON output.
    sci = (0.0 if math.isclose(mean_single_mfe, 0.0, abs_tol=1e-12)
           else consensus_mfe / mean_single_mfe)
    return {
        "sci": sci,
        "sci_consensus_mfe": consensus_mfe,
        "sci_mean_single_mfe": mean_single_mfe,
        "sci_single_mfes": individual_mfes,
        "sci_vienna_version": versions,
        "sci_vienna_expected_version": config.expected_version,
        "sci_rnaalifold": rnaalifold,
        "sci_rnafold": rnafold,
        "sci_reason": None,
    }


def derive_consensus_structure(alignment: StructuralAlignment,
                               config: SCIConfig) -> str:
    """Derive a common dot-bracket structure with the configured RNAalifold.

    External aligners such as MAFFT intentionally write only an alignment.  The
    original DAFS paper evaluated those alignments after a common structure was
    predicted, so the comparison runner uses this helper to make that step
    explicit and version-pinned.
    """
    rnaalifold = _resolve_executable(config.rnaalifold)
    output, _ = _run_sci_command(
        rnaalifold, ["--noPS", "--input-format=F"],
        _sci_fasta(alignment, aligned=True), Path.cwd(), config.timeout_seconds)
    structure_symbols = set(".()[]{}<>")
    for line in output.splitlines():
        for token in line.strip().split():
            if (len(token) == alignment.columns and
                    any(symbol in structure_symbols for symbol in token) and
                    set(token) <= structure_symbols | set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")):
                # Validate WUSS/pseudoknot brackets before returning the token.
                structure_pairs(token)
                return token
    tail = "\n".join(output.splitlines()[-12:])
    raise ValueError(f"could not parse RNAalifold consensus structure\n{tail}")


def with_consensus_structure(alignment: StructuralAlignment,
                             config: SCIConfig) -> StructuralAlignment:
    """Return an alignment with an RNAalifold consensus when none is present."""
    if alignment.structure is not None:
        return alignment
    return StructuralAlignment(
        alignment.sequences,
        derive_consensus_structure(alignment, config),
        alignment.structures,
    )


def structure_pairs(structure: str) -> set[tuple[int, int]]:
    open_to_close = {"(": ")", "[": "]", "{": "}", "<": ">"}
    open_to_close.update({chr(code): chr(code + 32)
                          for code in range(ord("A"), ord("Z") + 1)})
    close_to_open = {close: opening for opening, close in open_to_close.items()}
    stacks: dict[str, list[int]] = {opening: [] for opening in open_to_close}
    pairs: set[tuple[int, int]] = set()
    for column, symbol in enumerate(structure):
        if symbol in open_to_close:
            stacks[symbol].append(column)
        elif symbol in close_to_open:
            opening = close_to_open[symbol]
            if not stacks[opening]:
                raise ValueError(f"unmatched structure symbol {symbol!r} at column {column}")
            pairs.add((stacks[opening].pop(), column))
        elif symbol not in ".,:_-~":
            raise ValueError(f"unsupported structure symbol {symbol!r} at column {column}")
    unmatched = [(opening, positions) for opening, positions in stacks.items() if positions]
    if unmatched:
        raise ValueError(f"unmatched structure openings: {unmatched[:3]}")
    return pairs


def residue_map(aligned_sequence: str) -> list[int | None]:
    mapping: list[int | None] = []
    residue = 0
    for symbol in aligned_sequence:
        if symbol in GAPS:
            mapping.append(None)
        else:
            mapping.append(residue)
            residue += 1
    return mapping


def aligned_residue_pairs(alignment: StructuralAlignment, first: str,
                          second: str) -> set[tuple[int, int]]:
    first_map = residue_map(alignment.sequences[first])
    second_map = residue_map(alignment.sequences[second])
    return {(i, j) for i, j in zip(first_map, second_map)
            if i is not None and j is not None}


def mapped_pairs(alignment: StructuralAlignment, name: str) -> set[tuple[int, int]]:
    structure = _structure_string(alignment, name)
    if structure is None:
        return set()
    mapping = residue_map(alignment.sequences[name])
    result = set()
    for left, right in structure_pairs(structure):
        i, j = mapping[left], mapping[right]
        if i is not None and j is not None:
            result.add((i, j))
    return result


def _structure_string(alignment: StructuralAlignment,
                      name: str) -> str | None:
    if alignment.structures is not None and name in alignment.structures:
        return alignment.structures[name]
    return alignment.structure


def cbp_matches(alignment: StructuralAlignment,
                names: Iterable[str]) -> set[tuple[str, int, int, str, int, int]]:
    ordered_names = sorted(names)
    mappings = {name: residue_map(alignment.sequences[name])
                for name in ordered_names}
    pairs_by_name = {
        name: (set() if (structure := _structure_string(alignment, name)) is None
               else structure_pairs(structure))
        for name in ordered_names
    }
    matches = set()
    for index, first in enumerate(ordered_names):
        for second in ordered_names[index + 1:]:
            first_map, second_map = mappings[first], mappings[second]
            for left, right in pairs_by_name[first] & pairs_by_name[second]:
                values = (first_map[left], first_map[right],
                          second_map[left], second_map[right])
                if all(value is not None for value in values):
                    i, j, k, l = values
                    matches.add((first, i, j, second, k, l))
    return matches


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _f1(tp: int, fp: int, fn: int) -> float | None:
    denominator = 2 * tp + fp + fn
    return 2 * tp / denominator if denominator else None


def score(predicted: StructuralAlignment,
          reference: StructuralAlignment,
          sci_config: SCIConfig | None = None) -> dict[str, object]:
    names = sorted(set(predicted.sequences) & set(reference.sequences))
    missing_predicted = sorted(set(reference.sequences) - set(predicted.sequences))
    extra_predicted = sorted(set(predicted.sequences) - set(reference.sequences))
    if len(names) < 2:
        raise ValueError("prediction and reference must share at least two sequences")

    alignment_tp = alignment_predicted = alignment_reference = 0
    for index, first in enumerate(names):
        for second in names[index + 1:]:
            predicted_pairs = aligned_residue_pairs(predicted, first, second)
            reference_pairs = aligned_residue_pairs(reference, first, second)
            alignment_tp += len(predicted_pairs & reference_pairs)
            alignment_predicted += len(predicted_pairs)
            alignment_reference += len(reference_pairs)

    result: dict[str, object] = {
        "schema_version": 1,
        "shared_sequences": len(names),
        "missing_predicted_sequences": missing_predicted,
        "extra_predicted_sequences": extra_predicted,
        "alignment_tp": alignment_tp,
        "alignment_predicted": alignment_predicted,
        "alignment_reference": alignment_reference,
        "sps": _ratio(alignment_tp, alignment_reference),
        "alignment_ppv": _ratio(alignment_tp, alignment_predicted),
        "sci": None,
        "sci_reason": ("not requested; provide a version-pinned ViennaRNA "
                       "SCI configuration"),
    }

    if sci_config is not None:
        result.update(calculate_sci(predicted, sci_config))

    if any(_structure_string(predicted, name) is None or
           _structure_string(reference, name) is None for name in names):
        result.update({
            "structure_tp": None, "structure_fp": None, "structure_fn": None,
            "structure_tn": None, "sensitivity": None, "ppv": None,
            "mcc": None, "cbp_tp": None, "cbp_fp": None, "cbp_fn": None,
            "cbp_precision": None, "cbp_recall": None, "cbp_f1": None,
        })
        return result

    structure_tp = structure_fp = structure_fn = structure_tn = 0
    for name in names:
        predicted_pairs = mapped_pairs(predicted, name)
        reference_pairs = mapped_pairs(reference, name)
        tp = len(predicted_pairs & reference_pairs)
        fp = len(predicted_pairs - reference_pairs)
        fn = len(reference_pairs - predicted_pairs)
        residues = sum(symbol not in GAPS for symbol in reference.sequences[name])
        possible = residues * (residues - 1) // 2
        structure_tp += tp
        structure_fp += fp
        structure_fn += fn
        structure_tn += possible - tp - fp - fn

    denominator = math.sqrt(
        (structure_tp + structure_fp) * (structure_tp + structure_fn) *
        (structure_tn + structure_fp) * (structure_tn + structure_fn))
    predicted_cbp = cbp_matches(predicted, names)
    reference_cbp = cbp_matches(reference, names)
    cbp_tp = len(predicted_cbp & reference_cbp)
    cbp_fp = len(predicted_cbp - reference_cbp)
    cbp_fn = len(reference_cbp - predicted_cbp)
    result.update({
        "structure_tp": structure_tp,
        "structure_fp": structure_fp,
        "structure_fn": structure_fn,
        "structure_tn": structure_tn,
        "sensitivity": _ratio(structure_tp, structure_tp + structure_fn),
        "ppv": _ratio(structure_tp, structure_tp + structure_fp),
        "mcc": ((structure_tp * structure_tn - structure_fp * structure_fn) /
                denominator) if denominator else None,
        "cbp_tp": cbp_tp,
        "cbp_fp": cbp_fp,
        "cbp_fn": cbp_fn,
        "cbp_precision": _ratio(cbp_tp, cbp_tp + cbp_fp),
        "cbp_recall": _ratio(cbp_tp, cbp_tp + cbp_fn),
        "cbp_f1": _f1(cbp_tp, cbp_fp, cbp_fn),
    })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--reference", type=Path,
                        help="reference alignment; optional for SCI-only scoring")
    parser.add_argument("--vienna-bin-dir", type=Path,
                        help="directory containing RNAalifold and RNAfold")
    parser.add_argument("--rnaalifold", type=Path,
                        help="RNAalifold executable for SCI")
    parser.add_argument("--rnafold", type=Path,
                        help="RNAfold executable for SCI")
    parser.add_argument("--vienna-version",
                        help="expected version substring for both SCI tools")
    parser.add_argument("--sci-timeout", type=float, default=3600.0,
                        help="per-tool SCI timeout in seconds (default: 3600)")
    parser.add_argument("--output", type=Path,
                        help="write JSON here instead of stdout")
    args = parser.parse_args()
    if args.vienna_bin_dir and (args.rnaalifold or args.rnafold):
        parser.error("--vienna-bin-dir cannot be combined with --rnaalifold/--rnafold")
    if bool(args.rnaalifold) != bool(args.rnafold):
        parser.error("--rnaalifold and --rnafold must be specified together")
    sci_config = None
    if args.vienna_bin_dir:
        sci_config = SCIConfig(
            args.vienna_bin_dir / "RNAalifold",
            args.vienna_bin_dir / "RNAfold",
            args.vienna_version,
            args.sci_timeout)
    elif args.rnaalifold and args.rnafold:
        sci_config = SCIConfig(
            args.rnaalifold, args.rnafold, args.vienna_version,
            args.sci_timeout)
    if args.reference is None and sci_config is None:
        parser.error("--reference is required unless SCI tools are configured")
    predicted = parse_alignment(args.prediction)
    reference = parse_alignment(args.reference) if args.reference else None
    result = (score(predicted, reference, sci_config)
              if reference is not None else calculate_sci(predicted, sci_config))
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
