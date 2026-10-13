#!/usr/bin/env python3
"""Download and prepare the public datasets used by the DAFS paper."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import random
import shutil
import tarfile
import time
import urllib.request
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"

SOURCES = {
    "murlet": {
        "url": ("https://web.archive.org/web/20071017090934id_/"
                "http://www.ncrna.org/software/murlet/murlet_dataset_tgz"),
        "path": RAW / "murlet_dataset.tgz",
        "sha256": "87fff699745f49c25af0876a65d07d1b0395eebb36235066a3869ee5697635b0",
        "paper": "10.1093/bioinformatics/btm146",
    },
    "rfam_11_seed": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/Rfam/11.0/Rfam.seed.gz",
        "path": RAW / "Rfam.seed.11.0.gz",
        "sha256": "900a2ec8dd9dcf7ceb49067ffb5f41e320baf43e2445a7056ecd79a673e92737",
        "paper": "10.1093/bioinformatics/bts612",
    },
}

LINEARTURBOFOLD = {
    "repository": "https://github.com/LinearFold/LinearTurboFold",
    "commit": "445808976857df01c57e5d9302a0bbee10515bb5",
    "data_tree": "e6dd90779d2273240bfa621d7ea3739dbeac92cc",
    "paper": "10.1073/pnas.2116269118",
    "collections": ("16S-25", "23S-5", "RNaseP-20", "SRP-10", "telomerase-5"),
}

DAFS_PAPER = {
    "doi": "10.1093/bioinformatics/bts612",
    "article": "https://academic.oup.com/bioinformatics/article/28/24/3218/246194",
    "availability_url": "http://www.ncrna.org/software/dafs/",
    "rfam_release": "11.0",
    "original_pkfree_count": 691,
    "original_pk_count": 82,
    "sequences_per_alignment": 10,
    "selection_rule": (
        "manually curated Rfam seed alignments with consensus structures "
        "published in the literature"
    ),
    "established_datasets": [
        {
            "id": "bralibase2",
            "paper": "Gardner et al. 2005",
            "expected_alignment_count": 599,
            "status": "prepared_in_murlet_archive",
            "manifest_collection_prefix": "BRAliBaseII-",
        },
        {
            "id": "murlet",
            "paper": "Kiryu et al. 2007",
            "expected_alignment_count": 273,
            "status": "prepared_in_murlet_archive",
            "manifest_collections": ["dataset1", "dataset2"],
        },
        {
            "id": "mastr",
            "paper": "Lindgreen et al. 2007",
            "expected_alignment_count": 52,
            "status": "historical_archive_not_recovered",
            "source_url": "http://mastr.binf.ku.dk/",
            "archive_url": (
                "https://web.archive.org/web/20151107002200/"
                "http://servers.binf.ku.dk/mastr/"),
        },
        {
            "id": "mxscarna",
            "paper": "Tabei et al. 2008",
            "expected_alignment_count": 1693,
            "status": "historical_archive_not_recovered",
            "source_url": "http://mxscarna.ncrna.org/",
            "archive_url": (
                "https://web.archive.org/web/20150526045159/"
                "http://mxscarna.ncrna.org/"),
        },
    ],
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def obtain(source: dict[str, object], supplied: Path | None = None) -> Path:
    destination = Path(source["path"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    if supplied is not None:
        shutil.copyfile(supplied, destination)
    if not destination.is_file():
        temporary = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(str(source["url"]), headers={
            "User-Agent": "DAFS benchmark dataset preparation/1",
        })
        with urllib.request.urlopen(request) as response, temporary.open("wb") as stream:
            shutil.copyfileobj(response, stream)
        temporary.replace(destination)
    actual = digest(destination)
    if actual != source["sha256"]:
        raise ValueError(f"checksum mismatch for {destination}: {actual}")
    return destination


def git_blob_digest(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def download_bytes(url: str, attempts: int = 5) -> bytes:
    request = urllib.request.Request(url, headers={
        "User-Agent": "DAFS benchmark dataset preparation/1",
    })
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except OSError:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def split_stockholm(text: str) -> list[str]:
    blocks = []
    for fragment in re.split(r"\n//\s*\n", text.replace("\r\n", "\n")):
        fragment = fragment.strip()
        if fragment:
            if not fragment.startswith("# STOCKHOLM 1.0"):
                raise ValueError("unexpected content outside a Stockholm alignment")
            blocks.append(fragment + "\n//\n")
    return blocks


def parse_stockholm(block: str) -> tuple[OrderedDict[str, str], str, str | None]:
    sequences: OrderedDict[str, str] = OrderedDict()
    structure = []
    accession = None
    for line in block.splitlines():
        if line.startswith("#=GF AC"):
            accession = line.split(maxsplit=2)[2]
        elif line.startswith("#=GC SS_cons"):
            structure.append(line.split(maxsplit=2)[2])
        elif line and not line.startswith("#") and line != "//":
            fields = line.split()
            if len(fields) >= 2:
                sequences[fields[0]] = sequences.get(fields[0], "") + fields[1]
    if not sequences or not structure:
        raise ValueError("Stockholm alignment lacks sequences or SS_cons")
    lengths = {len(sequence) for sequence in sequences.values()}
    lengths.add(len("".join(structure)))
    if len(lengths) != 1:
        raise ValueError(f"inconsistent Stockholm alignment lengths: {sorted(lengths)}")
    return sequences, "".join(structure), accession


def parse_rfam_entry(block: str) -> dict[str, object]:
    """Parse one Rfam 11.0 seed block, retaining its curation metadata."""
    sequences: OrderedDict[str, str] = OrderedDict()
    structure: list[str] = []
    accession = None
    family_id = None
    structure_source: list[str] = []
    for line in block.splitlines():
        if line.startswith("#=GF AC"):
            accession = line.split(maxsplit=2)[2]
        elif line.startswith("#=GF ID"):
            family_id = line.split(maxsplit=2)[2]
        elif line.startswith("#=GF SS"):
            structure_source.append(line.split(maxsplit=2)[2])
        elif line.startswith("#=GC SS_cons"):
            structure.append(line.split(maxsplit=2)[2])
        elif line and not line.startswith("#") and line != "//":
            fields = line.split()
            if len(fields) >= 2:
                sequences[fields[0]] = sequences.get(fields[0], "") + fields[1]
    if not accession or not sequences or not structure:
        raise ValueError("Rfam seed block lacks accession, sequences, or SS_cons")
    consensus = "".join(structure)
    lengths = {len(sequence) for sequence in sequences.values()}
    lengths.add(len(consensus))
    if len(lengths) != 1:
        raise ValueError(
            f"inconsistent Rfam seed alignment lengths for {accession}: {sorted(lengths)}")
    return {
        "accession": accession,
        "family_id": family_id,
        "structure_source": "; ".join(structure_source),
        "sequences": sequences,
        "structure": consensus,
    }


def is_pseudoknotted_wuss(structure: str) -> bool:
    """Recognize the non-crossing WUSS pairs used by Rfam 11.0.

    In the Rfam seed files ``<>`` are ordinary paired symbols.  Alphabetic
    WUSS pair labels (and the additional bracket classes) denote crossings.
    """
    return bool(re.search(r"[A-Za-z\[\]{}]", structure))


def _rfam_entries(seed_path: Path) -> list[dict[str, object]]:
    with gzip.open(seed_path, "rt", encoding="latin-1") as stream:
        text = stream.read()
    return [parse_rfam_entry(block) for block in split_stockholm(text)]


def _selected_stockholm(entry: dict[str, object], names: list[str]) -> str:
    sequences = entry["sequences"]
    assert isinstance(sequences, OrderedDict)
    lines = ["# STOCKHOLM 1.0", f"#=GF AC   {entry['accession']}"]
    if entry.get("family_id"):
        lines.append(f"#=GF ID   {entry['family_id']}")
    if entry.get("structure_source"):
        lines.append(f"#=GF SS   {entry['structure_source']}")
    lines.append("")
    for name in names:
        lines.append(f"{name:<36} {sequences[name]}")
    lines.extend([
        f"#=GC SS_cons{'':>25} {entry['structure']}",
        "//",
        "",
    ])
    return "\n".join(lines)


def _write_rfam_case(target: Path, entry: dict[str, object], names: list[str],
                     dataset_id: str, kind: str, ordinal: int) -> dict[str, object]:
    inputs = target / "inputs"
    references = target / "references"
    inputs.mkdir(parents=True, exist_ok=True)
    references.mkdir(parents=True, exist_ok=True)
    input_path = inputs / f"{dataset_id}.fa"
    reference_path = references / f"{dataset_id}.sto"
    sequences = entry["sequences"]
    assert isinstance(sequences, OrderedDict)
    with input_path.open("w", encoding="utf-8") as stream:
        for name in names:
            unaligned = re.sub(r"[.\-~_]", "", sequences[name]).upper()
            stream.write(f">{name}\n{unaligned}\n")
    reference_path.write_text(
        _selected_stockholm(entry, names), encoding="utf-8")
    lengths = [len(re.sub(r"[.\-~_]", "", sequences[name])) for name in names]
    return {
        "id": dataset_id,
        "input": str(input_path.relative_to(target)),
        "reference": str(reference_path.relative_to(target)),
        "family_accession": entry["accession"],
        "family_id": entry.get("family_id"),
        "kind": kind,
        "sample_ordinal": ordinal,
        "sequence_names": names,
        "sequence_count": len(names),
        "min_sequence_length": min(lengths),
        "max_sequence_length": max(lengths),
    }


def _sample_rfam_names(entry: dict[str, object], ordinal: int,
                       sample_size: int) -> list[str]:
    sequences = entry["sequences"]
    assert isinstance(sequences, OrderedDict)
    names = list(sequences)
    if len(names) < sample_size:
        raise ValueError(
            f"Rfam family {entry['accession']} has only {len(names)} sequences")
    # The original ten-sequence draw is no longer present in the public
    # distribution.  Keep the reconstruction deterministic and independent of
    # Python's process hash randomization.
    seed_bytes = f"{entry['accession']}:{ordinal}".encode("ascii")
    seed = int(hashlib.sha256(seed_bytes).hexdigest()[:16], 16)
    selected = sorted(random.Random(seed).sample(names, sample_size),
                      key=names.index)
    return selected


def prepare_rfam_paper(seed_path: Path,
                       selection_path: Path | None = None) -> dict[str, object]:
    """Prepare a transparent Rfam reconstruction and a reference-backed long set.

    A JSON selection file can replace the reconstruction when the historical
    DAFS ten-sequence list is recovered.  Its records are objects with
    ``id``, ``kind`` (``pkfree`` or ``pk``), ``accession`` and
    ``sequence_names`` fields.
    """
    entries = {str(entry["accession"]): entry for entry in _rfam_entries(seed_path)}
    selected_records = None
    if selection_path is not None:
        selected_records = json.loads(selection_path.read_text(encoding="utf-8"))
        if not isinstance(selected_records, list):
            raise ValueError("paper selection JSON must contain a list of records")
        target = ROOT / "rfam11-paper-exact"
        status = "exact_selection_file"
    else:
        target = ROOT / "rfam11-paper-reconstructed"
        status = "reconstructed_one_sample_per_eligible_family"

    eligible_entries: dict[str, str] = {}
    eligible = {"pkfree": 0, "pk": 0}
    for accession in sorted(entries):
        entry = entries[accession]
        sequences = entry["sequences"]
        assert isinstance(sequences, OrderedDict)
        source = str(entry["structure_source"])
        if "Published" not in source or len(sequences) < DAFS_PAPER["sequences_per_alignment"]:
            continue
        kind = "pk" if is_pseudoknotted_wuss(str(entry["structure"])) else "pkfree"
        eligible_entries[accession] = kind
        eligible[kind] += 1

    if selected_records is not None:
        expected_counts = {
            "pkfree": int(DAFS_PAPER["original_pkfree_count"]),
            "pk": int(DAFS_PAPER["original_pk_count"]),
        }
        if any(not isinstance(record, dict) for record in selected_records):
            raise ValueError("every paper selection record must be an object")
        seen_ids: set[str] = set()
        selected_counts = {"pkfree": 0, "pk": 0}
        for record in selected_records:
            missing = {"id", "kind", "accession", "sequence_names"} - set(record)
            if missing:
                raise ValueError(
                    f"paper selection record is missing fields: {sorted(missing)}")
            kind = str(record["kind"])
            if kind not in selected_counts:
                raise ValueError(f"unknown paper selection kind: {kind}")
            record_id = safe_id(str(record["id"]))
            if record_id in seen_ids:
                raise ValueError(f"duplicate paper selection id: {record_id}")
            seen_ids.add(record_id)
            accession = str(record["accession"])
            if accession not in eligible_entries:
                raise ValueError(
                    f"paper selection accession is not an eligible Rfam entry: {accession}")
            if eligible_entries[accession] != kind:
                raise ValueError(
                    f"selection kind mismatch for {accession}: "
                    f"{kind} != {eligible_entries[accession]}")
            selected_counts[kind] += 1
        if selected_counts != expected_counts:
            raise ValueError(
                "exact DAFS paper selection must contain "
                f"{expected_counts['pkfree']} pkfree and "
                f"{expected_counts['pk']} pk records; got {selected_counts}")

    datasets = []
    for accession in sorted(eligible_entries):
        entry = entries[accession]
        sequences = entry["sequences"]
        assert isinstance(sequences, OrderedDict)
        kind = eligible_entries[accession]
        if selected_records is not None:
            records = [record for record in selected_records
                       if str(record.get("accession")) == accession]
            for ordinal, record in enumerate(records):
                names = [str(name) for name in record.get("sequence_names", [])]
                if str(record.get("kind")) != kind:
                    raise ValueError(
                        f"selection kind mismatch for {accession}: {record.get('kind')} != {kind}")
                if len(names) != DAFS_PAPER["sequences_per_alignment"]:
                    raise ValueError(f"selection {record.get('id')} must contain ten sequences")
                if len(set(names)) != len(names) or any(name not in sequences for name in names):
                    raise ValueError(f"invalid sequence names in selection {record.get('id')}")
                dataset_id = safe_id(str(record["id"]))
                datasets.append(_write_rfam_case(
                    target, entry, names, dataset_id, kind, ordinal))
        else:
            names = _sample_rfam_names(
                entry, 0, DAFS_PAPER["sequences_per_alignment"])
            dataset_id = safe_id(f"rfam11-{kind}-{accession}-sample01")
            datasets.append(_write_rfam_case(target, entry, names, dataset_id, kind, 0))

    manifest = {
        "schema_version": 1,
        "name": "DAFS paper Rfam 11.0 benchmark",
        "source": SOURCES["rfam_11_seed"]["url"],
        "source_sha256": SOURCES["rfam_11_seed"]["sha256"],
        "paper": DAFS_PAPER,
        "selection_status": status,
        "exact_selection_available": selection_path is not None,
        "selection_file": (str(selection_path.resolve())
                            if selection_path is not None else None),
        "eligible_family_counts": eligible,
        "generated_alignment_counts": {
            "pkfree": sum(dataset["kind"] == "pkfree" for dataset in datasets),
            "pk": sum(dataset["kind"] == "pk" for dataset in datasets),
        },
        "datasets": datasets,
    }
    target.mkdir(parents=True, exist_ok=True)
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    long_target = ROOT / "rfam11-long"
    long_datasets = []
    for accession in sorted(entries):
        entry = entries[accession]
        sequences = entry["sequences"]
        assert isinstance(sequences, OrderedDict)
        if ("Published" not in str(entry["structure_source"]) or
                len(sequences) < DAFS_PAPER["sequences_per_alignment"]):
            continue
        names = list(sequences)
        lengths = [len(re.sub(r"[.\-~_]", "", value)) for value in sequences.values()]
        if max(lengths) < 800:
            continue
        kind = "pk" if is_pseudoknotted_wuss(str(entry["structure"])) else "pkfree"
        sample = _sample_rfam_names(
            entry, 0, DAFS_PAPER["sequences_per_alignment"])
        dataset_id = safe_id(f"rfam11-long-{kind}-{accession}-sample01")
        long_datasets.append(_write_rfam_case(
            long_target, entry, sample, dataset_id, kind, 0))
    long_manifest = {
        "schema_version": 1,
        "name": "Rfam 11.0 published-structure long-sequence set",
        "source": SOURCES["rfam_11_seed"]["url"],
        "source_sha256": SOURCES["rfam_11_seed"]["sha256"],
        "selection_rule": "published structures, at least ten seed sequences, max ungapped length >= 800",
        "reference_based": True,
        "datasets": long_datasets,
    }
    long_target.mkdir(parents=True, exist_ok=True)
    (long_target / "manifest.json").write_text(
        json.dumps(long_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"paper": manifest, "long": long_manifest}


def safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-.")
    if not value:
        raise ValueError("empty dataset identifier")
    return value


def prepare_murlet(archive_path: Path) -> dict[str, object]:
    target = ROOT / "murlet"
    inputs = target / "inputs"
    references = target / "references"
    inputs.mkdir(parents=True, exist_ok=True)
    references.mkdir(parents=True, exist_ok=True)
    datasets = []
    collection_counts: dict[str, int] = {}
    seen_ids: set[str] = set()

    with tarfile.open(archive_path, "r:gz") as archive:
        members = sorted((
            member for member in archive.getmembers()
            if member.isfile() and (
                member.name.endswith("dataset1_stockholm.txt")
                or member.name.endswith("dataset2_stockholm.txt")
                or "/BRAliBaseII/" in member.name and member.name.endswith("_stockholm.txt")
            )
        ), key=lambda member: member.name)
        for member in members:
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read {member.name}")
            text = extracted.read().decode("utf-8")
            if "/BRAliBaseII/" in member.name:
                collection = "BRAliBaseII-" + Path(member.name).stem.removesuffix("_stockholm")
            else:
                collection = Path(member.name).stem.removesuffix("_stockholm")
            blocks = split_stockholm(text)
            collection_counts[collection] = len(blocks)
            for index, block in enumerate(blocks, 1):
                sequences, _structure, accession = parse_stockholm(block)
                suffix = accession or f"{index:03d}"
                dataset_id = safe_id(f"murlet-{collection}-{suffix}")
                if dataset_id in seen_ids:
                    dataset_id = safe_id(f"{dataset_id}-{index:03d}")
                seen_ids.add(dataset_id)
                input_path = inputs / f"{dataset_id}.fa"
                reference_path = references / f"{dataset_id}.sto"
                with input_path.open("w", encoding="utf-8") as stream:
                    for name, aligned in sequences.items():
                        unaligned = re.sub(r"[.\-~_]", "", aligned)
                        stream.write(f">{name}\n{unaligned}\n")
                reference_path.write_text(block, encoding="utf-8")
                datasets.append({
                    "id": dataset_id,
                    "input": str(input_path.relative_to(target)),
                    "reference": str(reference_path.relative_to(target)),
                    "collection": collection,
                    "accession": accession,
                    "sequence_count": len(sequences),
                })

    manifest = {
        "schema_version": 1,
        "name": "Murlet paper benchmark archive",
        "source": SOURCES["murlet"]["url"],
        "source_sha256": SOURCES["murlet"]["sha256"],
        "paper": SOURCES["murlet"]["paper"],
        "collection_counts": collection_counts,
        "datasets": datasets,
    }
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def parse_linearturbofold_seq(data: bytes, source: str) -> str:
    chunks = [
        line.strip() for line in data.decode("utf-8").splitlines()
        if re.fullmatch(r"[ACGUTNacgutn]+1?", line.strip())
    ]
    if not chunks:
        raise ValueError(f"no RNA sequence in {source}")
    sequence = "".join(chunks).upper()
    if sequence.endswith("1"):
        sequence = sequence[:-1]
    sequence = sequence.replace("T", "U")
    if not sequence or re.search(r"[^ACGUN]", sequence):
        raise ValueError(f"invalid RNA sequence in {source}")
    return sequence


def prepare_linearturbofold_nonviral() -> dict[str, object]:
    """Fetch only official nonviral sequence groups and make DAFS FASTAs."""
    commit = LINEARTURBOFOLD["commit"]
    raw_root = RAW / "linearturbofold-nonviral"
    target = ROOT / "linearturbofold-nonviral"

    # Keep repeated preparation runs offline once the immutable GitHub tree
    # and its blobs have already been verified.  This matters on clusters
    # where compute nodes can read the checked-in/cache volume but cannot
    # reach api.github.com.
    cached_path = target / "manifest.json"
    if cached_path.is_file() and raw_root.is_dir():
        try:
            cached = json.loads(cached_path.read_text(encoding="utf-8"))
            if (cached.get("source_commit") == commit and
                    cached.get("source_data_tree") == LINEARTURBOFOLD["data_tree"]):
                valid = True
                for dataset in cached["datasets"]:
                    input_path = target / str(dataset["input"])
                    if not input_path.is_file():
                        valid = False
                        break
                    for source_path, expected_sha in zip(
                            dataset["source_files"], dataset["source_git_blob_sha1"]):
                        local_path = raw_root / Path(source_path).relative_to("data")
                        if (not local_path.is_file() or
                                git_blob_digest(local_path.read_bytes()) != expected_sha):
                            valid = False
                            break
                    if not valid:
                        break
                if valid:
                    return cached
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            # A stale or incomplete cache is repaired by the network path
            # below when network access is available.
            pass

    api = f"https://api.github.com/repos/LinearFold/LinearTurboFold/git/trees/{commit}?recursive=1"
    tree = json.loads(download_bytes(api))
    prefixes = tuple(f"data/{name}/" for name in LINEARTURBOFOLD["collections"])
    entries = sorted(
        (entry for entry in tree["tree"]
         if entry["type"] == "blob" and entry["path"].startswith(prefixes)
         and entry["path"].endswith(".seq")),
        key=lambda entry: entry["path"],
    )
    if not entries:
        raise ValueError("LinearTurboFold tree contains no selected .seq files")

    inputs = target / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for entry in entries:
        relative = Path(entry["path"]).relative_to("data")
        destination = raw_root / relative
        if destination.is_file():
            data = destination.read_bytes()
        else:
            url = f"https://raw.githubusercontent.com/LinearFold/LinearTurboFold/{commit}/{entry['path']}"
            data = download_bytes(url)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
        if git_blob_digest(data) != entry["sha"]:
            raise ValueError(f"Git blob checksum mismatch for {entry['path']}")
        collection, group = relative.parts[:2]
        grouped.setdefault((collection, group), []).append({
            "name": Path(relative.name).stem,
            "sequence": parse_linearturbofold_seq(data, entry["path"]),
            "source_path": entry["path"],
            "git_blob_sha1": entry["sha"],
        })

    datasets = []
    for (collection, group), records in sorted(grouped.items()):
        dataset_id = safe_id(f"linearturbofold-{collection}-{group}")
        input_path = inputs / f"{dataset_id}.fa"
        with input_path.open("w", encoding="utf-8") as stream:
            for record in records:
                stream.write(f">{record['name']}\n{record['sequence']}\n")
        lengths = [len(str(record["sequence"])) for record in records]
        datasets.append({
            "id": dataset_id,
            "input": str(input_path.relative_to(target)),
            "collection": collection,
            "group": group,
            "sequence_count": len(records),
            "min_sequence_length": min(lengths),
            "max_sequence_length": max(lengths),
            "total_nucleotides": sum(lengths),
            "source_files": [record["source_path"] for record in records],
            "source_git_blob_sha1": [record["git_blob_sha1"] for record in records],
        })
    manifest = {
        "schema_version": 1,
        "name": "LinearTurboFold official nonviral evaluation sequences",
        "repository": LINEARTURBOFOLD["repository"],
        "source_commit": commit,
        "source_data_tree": LINEARTURBOFOLD["data_tree"],
        "paper": LINEARTURBOFOLD["paper"],
        "excluded_collections": ["HIV-5", "sars-cov-2_data"],
        "reference_alignment_in_repository": False,
        "datasets": datasets,
    }
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--murlet-archive", type=Path,
                        help="Import an already downloaded Murlet archive")
    parser.add_argument("--rfam-seed", type=Path,
                        help="Import an already downloaded Rfam 11.0 seed")
    parser.add_argument(
        "--paper-selection", type=Path,
        help=("JSON file containing the historical DAFS ten-sequence selection; "
              "without it, generate a deterministic one-sample-per-family reconstruction"))
    parser.add_argument("--skip-linearturbofold", action="store_true",
                        help="Do not download the LinearTurboFold nonviral datasets")
    args = parser.parse_args()

    murlet = obtain(SOURCES["murlet"], args.murlet_archive)
    rfam = obtain(SOURCES["rfam_11_seed"], args.rfam_seed)
    manifest = prepare_murlet(murlet)
    rfam_manifests = prepare_rfam_paper(rfam, args.paper_selection)
    paper_directory = (
        "rfam11-paper-exact" if args.paper_selection else
        "rfam11-paper-reconstructed")
    linearturbofold = None if args.skip_linearturbofold else prepare_linearturbofold_nonviral()
    provenance = {
        "schema_version": 2,
        "sources": {
            name: {**{key: value for key, value in source.items() if key != "path"},
                   "local_path": str(Path(source["path"]).relative_to(ROOT)),
                   "verified": True}
            for name, source in SOURCES.items()
        },
        "prepared": {
            "murlet_dataset_count": len(manifest["datasets"]),
            "murlet_manifest": "murlet/manifest.json",
            "rfam_11_seed": str(rfam.relative_to(ROOT)),
            "dafs_paper_manifest": f"{paper_directory}/manifest.json",
            "dafs_paper_selection_status": rfam_manifests["paper"]["selection_status"],
            "dafs_pkfree_exact_selection_available": bool(
                rfam_manifests["paper"]["exact_selection_available"]),
            "dafs_paper_generated_alignment_counts": rfam_manifests["paper"][
                "generated_alignment_counts"],
            "rfam11_long_manifest": "rfam11-long/manifest.json",
            "rfam11_long_dataset_count": len(rfam_manifests["long"]["datasets"]),
            "linearturbofold_nonviral_manifest": (
                "linearturbofold-nonviral/manifest.json" if linearturbofold else None),
            "linearturbofold_nonviral_dataset_count": (
                len(linearturbofold["datasets"]) if linearturbofold else 0),
        },
        "linearturbofold": LINEARTURBOFOLD,
        "dafs_paper": DAFS_PAPER,
    }
    (ROOT / "sources.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Murlet datasets: {len(manifest['datasets'])}")
    print(f"Manifest: {ROOT / 'murlet' / 'manifest.json'}")
    print(f"Rfam 11.0 seed: {rfam}")
    print(f"DAFS paper manifest: {ROOT / paper_directory / 'manifest.json'}")
    print(f"Rfam long datasets: {len(rfam_manifests['long']['datasets'])}")
    if linearturbofold:
        print(f"LinearTurboFold nonviral datasets: {len(linearturbofold['datasets'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
