#!/usr/bin/env python3
"""Import the original DAFS benchmark archives supplied by the user.

The importer reads ZIP and tar archives without extracting archive paths. It
writes normalized benchmark inputs/references below ``original-paper/`` and
records source archive checksums in a manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import tarfile
import zipfile
from collections import Counter, OrderedDict
from pathlib import Path, PurePosixPath
from typing import BinaryIO


ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = Path.home() / "dafs-data"
OUTPUT = ROOT / "original-paper"
GAPS = ".-~_"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_id(value: str) -> str:
    result = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-.")
    if not result:
        raise ValueError(f"empty dataset id from {value!r}")
    return result


def parse_fasta(data: bytes, source: str) -> OrderedDict[str, str]:
    records: OrderedDict[str, list[str]] = OrderedDict()
    name: str | None = None
    for line in data.decode("utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            if not name or name in records:
                raise ValueError(f"empty or duplicate FASTA name in {source}: {name!r}")
            records[name] = []
        elif name is None:
            raise ValueError(f"sequence before FASTA header in {source}")
        else:
            records[name].append("".join(line.split()))
    sequences = OrderedDict((name, "".join(parts)) for name, parts in records.items())
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError(f"no sequences or empty sequence in {source}")
    return sequences


def parse_mxscarna(data: bytes, source: str) -> tuple[OrderedDict[str, str], str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    structure = None
    for line in data.decode("utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        if fields[0] == "ss_cons":
            structure = fields[1]
        else:
            if fields[0] in sequences:
                raise ValueError(f"duplicate MXSCARNA sequence in {source}: {fields[0]}")
            sequences[fields[0]] = fields[1]
    if not sequences or structure is None:
        raise ValueError(f"MXSCARNA reference lacks sequences or ss_cons: {source}")
    lengths = {len(value) for value in sequences.values()} | {len(structure)}
    if len(lengths) != 1:
        raise ValueError(f"inconsistent MXSCARNA alignment columns in {source}: {lengths}")
    return sequences, structure


def write_fasta(path: Path, records: OrderedDict[str, str],
                structure: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        if structure is not None:
            stream.write(f">SS_cons\n{structure}\n")
        for name, sequence in records.items():
            stream.write(f">{name}\n{sequence}\n")


def add_dataset(datasets: list[dict[str, object]], category: str, source_id: str,
                input_records: OrderedDict[str, str], reference_records: OrderedDict[str, str] | None,
                reference_structure: str | None, source_archive: str,
                source_member: str, output_prefix: str = "original-paper") -> None:
    dataset_id = safe_id(f"{category}-{source_id}")
    input_path = OUTPUT / "inputs" / category / f"{dataset_id}.fa"
    reference_path = OUTPUT / "references" / category / f"{dataset_id}.fa"
    if reference_records is not None:
        if list(input_records) != list(reference_records):
            if set(input_records) != set(reference_records):
                raise ValueError(f"input/reference sequence names differ for {dataset_id}")
            reference_records = OrderedDict(
                (name, reference_records[name]) for name in input_records)
        if len({len(value) for value in reference_records.values()}) != 1:
            raise ValueError(f"reference alignment has unequal columns: {dataset_id}")
        if reference_structure is not None and len(reference_structure) != len(
                next(iter(reference_records.values()))):
            raise ValueError(f"SS_cons length differs from alignment: {dataset_id}")
        for name, raw in input_records.items():
            aligned = reference_records[name]
            ungapped = "".join(ch for ch in aligned if ch not in GAPS).upper()
            if ungapped != "".join(ch for ch in raw if ch not in GAPS).upper():
                raise ValueError(f"input and reference residues differ for {dataset_id}/{name}")
        write_fasta(reference_path, reference_records, reference_structure)
        reference_rel = str(reference_path.relative_to(OUTPUT))
    else:
        reference_rel = None
    write_fasta(input_path, input_records)
    lengths = [len("".join(ch for ch in value if ch not in GAPS))
               for value in input_records.values()]
    datasets.append({
        "id": dataset_id,
        "collection": category,
        "input": str(input_path.relative_to(OUTPUT)),
        "reference": reference_rel,
        "sequence_count": len(input_records),
        "min_sequence_length": min(lengths),
        "max_sequence_length": max(lengths),
        "source_archive": source_archive,
        "source_member": source_member,
    })


def import_rfam_archive(source: Path, archive_name: str, category: str,
                        datasets: list[dict[str, object]]) -> None:
    with zipfile.ZipFile(source / archive_name) as archive:
        names = set(archive.namelist())
        raw_names = sorted(name for name in names
                           if name.startswith(f"{category}/raw/") and name.endswith(".fa"))
        reference_ids = {
            name.removeprefix(f"{category}/ref/").removesuffix(".fa")
            for name in names if name.startswith(f"{category}/ref/") and name.endswith(".fa")
        }
        raw_ids = {
            name.removeprefix(f"{category}/raw/").removesuffix(".fa")
            for name in raw_names
        }
        if raw_ids != reference_ids:
            raise ValueError(f"raw/reference case mismatch in {archive_name}")
        expected = 691 if category == "pkfree" else 82
        if len(raw_names) != expected:
            raise ValueError(f"expected {expected} {category} cases, found {len(raw_names)}")
        for raw_name in raw_names:
            source_id = raw_name.removeprefix(f"{category}/raw/").removesuffix(".fa")
            reference_name = f"{category}/ref/{source_id}.fa"
            raw = parse_fasta(archive.read(raw_name), raw_name)
            reference = parse_fasta(archive.read(reference_name), reference_name)
            ss_cons = reference.pop("SS_cons", None)
            if ss_cons is None:
                raise ValueError(f"reference lacks SS_cons: {reference_name}")
            add_dataset(datasets, category, source_id, raw, reference, ss_cons,
                        archive_name, raw_name)


def import_mastr(source: Path, datasets: list[dict[str, object]]) -> None:
    archive_name = "datasets_mastr.tar.gz"
    with tarfile.open(source / archive_name, "r:gz") as archive:
        members = sorted((member for member in archive.getmembers()
                          if member.isfile() and member.name.endswith(".fa")),
                         key=lambda member: member.name)
        if len(members) != 52:
            raise ValueError(f"expected 52 MASTR data files, found {len(members)}")
        for member in members:
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"could not read {member.name}")
            records = parse_fasta(stream.read(), member.name)
            add_dataset(datasets, "mastr", PurePosixPath(member.name).stem,
                        records, None, None, archive_name, member.name)


def import_mxscarna(source: Path, datasets: list[dict[str, object]]) -> None:
    archive_name = "dataset_mxscarna.tar.gz"
    with tarfile.open(source / archive_name, "r:gz") as archive:
        members = sorted((member for member in archive.getmembers()
                          if member.isfile() and member.name.endswith(".fa")),
                         key=lambda member: member.name)
        if len(members) != 1693:
            raise ValueError(f"expected 1693 MXSCARNA alignments, found {len(members)}")
        for member in members:
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"could not read {member.name}")
            aligned, structure = parse_mxscarna(stream.read(), member.name)
            raw = OrderedDict((name, "".join(ch for ch in sequence if ch not in GAPS))
                              for name, sequence in aligned.items())
            add_dataset(datasets, "mxscarna", PurePosixPath(member.name).stem,
                        raw, aligned, structure, archive_name, member.name)


def import_bralibase21(source: Path, datasets: list[dict[str, object]]) -> dict[str, object]:
    archive_name = "bralibase2.1-selected.zip"
    archive_path = source / archive_name
    file_total = 0
    with zipfile.ZipFile(archive_path) as outer:
        for nested_name in sorted(outer.namelist()):
            if not nested_name.endswith(".tgz"):
                continue
            collection = PurePosixPath(nested_name).stem
            with tarfile.open(fileobj=io.BytesIO(outer.read(nested_name)), mode="r:gz") as nested:
                files = {member.name: member for member in nested.getmembers()
                         if member.isfile()}
                raw_members = sorted(name for name in files if name.endswith(".raw.fa"))
                reference_names = {name.removesuffix(".raw.fa") + ".ref.fa"
                                   for name in raw_members}
                present_references = {name for name in reference_names if name in files}
                if present_references != reference_names:
                    raise ValueError(f"BRAliBase raw/reference mismatch in {nested_name}")
                for member_name in raw_members:
                    member = files[member_name]
                    reference_member = files[member_name.removesuffix(".raw.fa") + ".ref.fa"]
                    raw_stream = nested.extractfile(member)
                    ref_stream = nested.extractfile(reference_member)
                    if raw_stream is None or ref_stream is None:
                        raise ValueError(f"could not read BRAliBase case {member_name}")
                    raw = parse_fasta(raw_stream.read(), member_name)
                    aligned = parse_fasta(ref_stream.read(), reference_member.name)
                    source_id = safe_id(f"{collection}-{PurePosixPath(member_name).stem}")
                    add_dataset(datasets, "bralibase2.1", source_id, raw,
                                aligned, None, archive_name,
                                f"{nested_name}:{member_name}")
                    file_total += 1
    return {
        "dataset_count": file_total,
        "source_archive": archive_name,
        "sha256": sha256(archive_path.read_bytes()),
        "description": "BRAliBase 2.1 selected cases with raw/ref FASTA pairs",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--skip-bralibase21", action="store_true",
                        help="skip importing the large 6-archive BRAliBase 2.1 collection")
    args = parser.parse_args()
    source = args.source_dir.expanduser().resolve()
    required = ["pkfree.zip", "pk.zip", "datasets_mastr.tar.gz",
                "dataset_mxscarna.tar.gz"]
    if not args.skip_bralibase21:
        required.extend(["bralibase2.1-selected.zip", "bralibase2.zip"])
    missing = [name for name in required if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing source archives in {source}: {missing}")

    datasets: list[dict[str, object]] = []
    import_rfam_archive(source, "pkfree.zip", "pkfree", datasets)
    import_rfam_archive(source, "pk.zip", "pk", datasets)
    import_mastr(source, datasets)
    import_mxscarna(source, datasets)
    bralibase = None
    if not args.skip_bralibase21:
        bralibase = import_bralibase21(source, datasets)

    counts = Counter(str(record["collection"]) for record in datasets)
    sources = {}
    for archive_name in required:
        path = source / archive_name
        sources[archive_name] = {"sha256": sha256(path.read_bytes()), "bytes": path.stat().st_size}
    manifest = {
        "schema_version": 1,
        "name": "Original DAFS paper and supplied structural RNA datasets",
        "source_directory": str(source),
        "source_archives": sources,
        "dataset_counts": dict(sorted(counts.items())),
        "bralibase2.1_summary": bralibase,
        "bralibase2_zip_note": (
            "The supplied BRAliBase 2.0 archive contains .ct/.struct annotation files, "
            "not standalone sequence/alignment FASTA cases; it is checksum-recorded "
            "but not expanded into runnable alignment datasets." if not args.skip_bralibase21 else None),
        "datasets": datasets,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Imported dataset cases: {len(datasets)}")
    print(f"Collection counts: {dict(sorted(counts.items()))}")
    print(f"Manifest: {OUTPUT / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
