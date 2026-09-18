#!/usr/bin/env python3
"""Download and prepare the public datasets used by the DAFS paper."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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
    for fragment in text.replace("\r\n", "\n").split("//"):
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

    raw_root = RAW / "linearturbofold-nonviral"
    target = ROOT / "linearturbofold-nonviral"
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
    parser.add_argument("--skip-linearturbofold", action="store_true",
                        help="Do not download the LinearTurboFold nonviral datasets")
    args = parser.parse_args()

    murlet = obtain(SOURCES["murlet"], args.murlet_archive)
    rfam = obtain(SOURCES["rfam_11_seed"], args.rfam_seed)
    manifest = prepare_murlet(murlet)
    linearturbofold = None if args.skip_linearturbofold else prepare_linearturbofold_nonviral()
    provenance = {
        "schema_version": 1,
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
            "dafs_pkfree_exact_selection_available": False,
            "linearturbofold_nonviral_manifest": (
                "linearturbofold-nonviral/manifest.json" if linearturbofold else None),
            "linearturbofold_nonviral_dataset_count": (
                len(linearturbofold["datasets"]) if linearturbofold else 0),
        },
        "linearturbofold": LINEARTURBOFOLD,
    }
    (ROOT / "sources.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Murlet datasets: {len(manifest['datasets'])}")
    print(f"Manifest: {ROOT / 'murlet' / 'manifest.json'}")
    print(f"Rfam 11.0 seed: {rfam}")
    if linearturbofold:
        print(f"LinearTurboFold nonviral datasets: {len(linearturbofold['datasets'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
