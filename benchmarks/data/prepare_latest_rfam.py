#!/usr/bin/env python3
"""Build reference-backed benchmarks from the current official Rfam seed."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
from collections import Counter, OrderedDict
from pathlib import Path

from prepare import (is_pseudoknotted_wuss, parse_rfam_entry, safe_id,
                     split_stockholm, _write_rfam_case)


RELEASE = "15.1"
URL = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{RELEASE}/Rfam.seed.gz"
SAMPLE_SIZE = 10
MIN_LONG_NT = 500
DATA_ROOT = Path(__file__).resolve().parent


def sample_names(accession: str, names: list[str]) -> list[str]:
    seed = int(hashlib.sha256(f"Rfam-{RELEASE}:{accession}".encode()).hexdigest()[:16], 16)
    return sorted(random.Random(seed).sample(names, SAMPLE_SIZE), key=names.index)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path,
                        default=DATA_ROOT / "raw" / f"Rfam.seed.{RELEASE}.gz")
    args = parser.parse_args()
    seed_path = args.seed.resolve()
    if not seed_path.is_file():
        raise FileNotFoundError(seed_path)
    source_sha = hashlib.sha256(seed_path.read_bytes()).hexdigest()
    entries = {}
    excluded = []
    import gzip
    with gzip.open(seed_path, "rt", encoding="latin-1") as stream:
        for block in split_stockholm(stream.read()):
            structure_track_count = sum(
                line.startswith("#=GC SS_cons") for line in block.splitlines())
            if structure_track_count != 1:
                accession = next((line.split(maxsplit=2)[2]
                                  for line in block.splitlines()
                                  if line.startswith("#=GF AC")), "unknown")
                excluded.append({
                    "accession": accession,
                    "reason": "expected exactly one #=GC SS_cons track",
                    "track_count": structure_track_count,
                })
                continue
            entry = parse_rfam_entry(block)
            entries[str(entry["accession"])] = entry

    current_root = DATA_ROOT / f"rfam{RELEASE}-current"
    long_root = DATA_ROOT / f"rfam{RELEASE}-long"
    cases = []
    long_cases = []
    eligible_counts = Counter()
    length_bins = Counter()
    for accession in sorted(entries):
        entry = entries[accession]
        sequences = entry["sequences"]
        assert isinstance(sequences, OrderedDict)
        if len(sequences) < SAMPLE_SIZE:
            continue
        if "Published" not in str(entry["structure_source"]):
            continue
        kind = "pk" if is_pseudoknotted_wuss(str(entry["structure"])) else "pkfree"
        eligible_counts[kind] += 1
        names = list(sequences)
        name_lengths = {
            name: len(re.sub(r"[.\-~_]", "", sequences[name]))
            for name in names
        }
        selected = sample_names(accession, names)
        case_id = safe_id(f"rfam{RELEASE}-{kind}-{accession}-sample10")
        case = _write_rfam_case(current_root, entry, selected, case_id, kind, 0)
        case["collection"] = kind
        cases.append(case)
        observed_max = max(name_lengths[name] for name in selected)
        length_bins[
            "<500" if observed_max < 500 else
            "500-999" if observed_max < 1000 else
            "1000-1999" if observed_max < 2000 else ">=2000"
        ] += 1

        longest_name = max(names, key=lambda name: (name_lengths[name], name))
        if name_lengths[longest_name] >= MIN_LONG_NT:
            selected_long = [longest_name]
            remainder = [name for name in names if name != longest_name]
            seed_value = int(hashlib.sha256(
                f"Rfam-{RELEASE}:{accession}:long".encode()).hexdigest()[:16], 16)
            selected_long.extend(random.Random(seed_value).sample(
                remainder, SAMPLE_SIZE - 1))
            selected_long = [name for name in names if name in set(selected_long)]
            long_id = safe_id(f"rfam{RELEASE}-long-{kind}-{accession}-n{ name_lengths[longest_name] }")
            long_case = _write_rfam_case(
                long_root, entry, selected_long, long_id, kind, 0)
            long_case["collection"] = kind
            long_case["longest_source_sequence_nt"] = name_lengths[longest_name]
            long_cases.append(long_case)

    if not cases:
        raise ValueError("Rfam seed produced no published-structure families")
    common = {
        "schema_version": 1,
        "rfam_release": RELEASE,
        "source": URL,
        "source_sha256": source_sha,
        "selection_rule": (
            f"Published SS_cons, at least {SAMPLE_SIZE} seed sequences; "
            "one deterministic sample of ten sequences per family"),
        "reference_based": True,
    }
    current_manifest = {
        **common,
        "name": f"Rfam {RELEASE} current published-structure benchmark",
        "eligible_family_counts": dict(sorted(eligible_counts.items())),
        "sampled_max_length_bins": dict(sorted(length_bins.items())),
        "excluded_malformed_structure_tracks": excluded,
        "datasets": cases,
    }
    long_manifest = {
        **common,
        "name": f"Rfam {RELEASE} long-sequence benchmark",
        "selection_rule": (
            f"Published SS_cons, at least {SAMPLE_SIZE} seed sequences, "
            f"longest ungapped sequence >= {MIN_LONG_NT} nt; each ten-sequence "
            "sample includes that longest sequence"),
        "minimum_longest_sequence_nt": MIN_LONG_NT,
        "datasets": long_cases,
    }
    for directory, manifest in ((current_root, current_manifest),
                                 (long_root, long_manifest)):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
    print(f"Rfam {RELEASE} seed SHA256: {source_sha}")
    print(f"Published-structure families: {dict(sorted(eligible_counts.items()))}")
    print(f"Current reference cases: {len(cases)}")
    print(f"Long reference cases (>= {MIN_LONG_NT} nt): {len(long_cases)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
