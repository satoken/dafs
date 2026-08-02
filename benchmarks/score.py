#!/usr/bin/env python3
"""Score a DAFS structural alignment against FASTA or Stockholm reference.

This dependency-free scorer reports SPS, secondary-structure sensitivity,
PPV and MCC, plus pair-pair (CBP) F1 for evaluating the RIBOSUM objective.
SCI is left to a separately version-pinned ViennaRNA evaluation.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


GAPS = frozenset("-.~_")


@dataclass(frozen=True)
class StructuralAlignment:
    sequences: OrderedDict[str, str]
    structure: str | None

    @property
    def columns(self) -> int:
        return len(next(iter(self.sequences.values())))


def _validate(sequences: OrderedDict[str, str], structure: str | None,
              source: Path) -> StructuralAlignment:
    if not sequences:
        raise ValueError(f"no sequences found in {source}")
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"inconsistent alignment lengths in {source}: {sorted(lengths)}")
    columns = next(iter(lengths))
    if structure is not None and len(structure) != columns:
        raise ValueError(
            f"SS_cons length {len(structure)} != alignment length {columns} in {source}")
    return StructuralAlignment(sequences, structure)


def parse_stockholm(path: Path) -> StructuralAlignment:
    fragments: OrderedDict[str, list[str]] = OrderedDict()
    structure: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line == "//":
            continue
        if line.startswith("#=GC SS_cons"):
            parts = line.split(maxsplit=2)
            if len(parts) != 3:
                raise ValueError(f"malformed SS_cons line in {path}: {raw}")
            structure.append(parts[2])
        elif line.startswith("#"):
            continue
        else:
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"malformed Stockholm sequence line in {path}: {raw}")
            fragments.setdefault(parts[0], []).append(parts[1])
    sequences = OrderedDict((name, "".join(parts))
                            for name, parts in fragments.items())
    return _validate(sequences, "".join(structure) or None, path)


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


def parse_alignment(path: str | Path) -> StructuralAlignment:
    source = Path(path)
    first = next((line.strip() for line in source.read_text(encoding="utf-8").splitlines()
                  if line.strip()), "")
    return (parse_stockholm(source) if first.startswith("# STOCKHOLM")
            else parse_fasta_alignment(source))


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
    if alignment.structure is None:
        return set()
    mapping = residue_map(alignment.sequences[name])
    result = set()
    for left, right in structure_pairs(alignment.structure):
        i, j = mapping[left], mapping[right]
        if i is not None and j is not None:
            result.add((i, j))
    return result


def cbp_matches(alignment: StructuralAlignment,
                names: Iterable[str]) -> set[tuple[str, int, int, str, int, int]]:
    if alignment.structure is None:
        return set()
    pairs = structure_pairs(alignment.structure)
    ordered_names = sorted(names)
    mappings = {name: residue_map(alignment.sequences[name])
                for name in ordered_names}
    matches = set()
    for index, first in enumerate(ordered_names):
        for second in ordered_names[index + 1:]:
            first_map, second_map = mappings[first], mappings[second]
            for left, right in pairs:
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
          reference: StructuralAlignment) -> dict[str, object]:
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
        "sci_reason": "requires separately version-pinned RNAalifold/RNAfold",
    }

    if predicted.structure is None or reference.structure is None:
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
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", type=Path,
                        help="write JSON here instead of stdout")
    args = parser.parse_args()
    result = score(parse_alignment(args.prediction), parse_alignment(args.reference))
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
