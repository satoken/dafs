#!/usr/bin/env python3
"""Exercise the CLI switch in progressive, final, and constrained decoding."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def check_stacked_aux(binary, root, with_ilp):
    """A below-threshold common pair must support its profitable neighbor."""
    fasta = root / "stacked.fa"
    fasta.write_text(">seq1\nAAAAAAAAAA\n>seq2\nAAAAAAAAAA\n")
    folding = root / "stacked.bp"
    folding.write_text("".join(f"> {i}\n1 10:0.95\n2 9:0.1\n10\n" for i in [1, 2]))
    alignment = root / "identity.mp"
    alignment.write_text("> 1 2\n" + "".join(f"{i} {i}:1\n" for i in range(1, 11)))
    for iterations in [100] + ([0] if with_ilp else []):
        for allowed in [False, True]:
            metrics = root / f"stacked-{iterations}-{allowed}.jsonl"
            command = [binary, "--fold-aux", str(folding), "--align-aux", str(alignment),
                       "--fold-th=0.2", "--fold-th1=0.2", "--align-th=0.01",
                       "--weight=1", "--ribosum-weight=0", "--fourway-pct=0",
                       "--fold-pct=0", "--align-pct=0", f"--max-iter={iterations}",
                       "--metrics-jsonl", str(metrics),
                       *(["--allow-isolated"] if allowed else []), str(fasta)]
            result = subprocess.run(command, text=True, capture_output=True,
                                    timeout=20, check=True)
            lines = result.stdout.splitlines()
            structure = lines[lines.index(">SS_cons") + 1]
            expected = "(........)" if allowed else "((......))"
            if structure != expected:
                raise AssertionError(f"stacked/{iterations}/{allowed}: {structure} != {expected}")
            events = [json.loads(line) for line in metrics.read_text().splitlines()]
            score = next(e for e in events if e["event"] == "run_summary")["progressive_score"]
            expected_score = 11.4 if allowed else 11.2
            if abs(score - expected_score) > 1e-4:
                raise AssertionError(f"stacked/{iterations}/{allowed}: score {score} != {expected_score}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--with-ilp", action="store_true")
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    cases = [("IPknot", model, []) for model in
             ["lpc", "lpv"] + (["CONTRAfold"] if args.with_ilp else [])]
    cases += [("Nussinov", model, []) for model in ["lpc", "lpv", "CONTRAfold"]]
    cases += [("Nussinov", "lpc", ["--bp-update1"]),
              ("Nussinov", "lpc", ["--dense-lagrangian"]),
              ("Nussinov", "lpv", ["--alifold", "--bp-update1"]),
              ("Nussinov", "lpc", ["--final-ribosum-weight=0.2"])]
    if args.with_ilp:
        cases += [("Nussinov", "lpc", ["--max-iter=0"])]
    with tempfile.TemporaryDirectory(prefix="dafs-allow-isolated-") as directory:
        root = Path(directory)
        input_path = root / "isolated.fa"
        # With a minimum hairpin length of three, only the terminal GC pair
        # can form. It cannot satisfy either backend's stacking constraints.
        input_path.write_text(">seq1\nGAAAAC\n>seq2\nGAAAAC\n")
        for case_index, (decoder, model, extra) in enumerate(cases):
            label = f"{decoder}/{model}/{extra}"
            outputs = {}
            scores = {}
            for mode, flags in (("default", []),
                                ("allowed", ["--allow-isolated"]),
                                ("false", ["--allow-isolated=false"])):
                metrics = root / f"{case_index}-{mode}.jsonl"
                command = [binary, "-a", "LinearAlign", "-s", model,
                           *(["--ipknot"] if decoder == "IPknot" else []),
                           "-g", "999999,999999",
                           "-G", "999999,999999", "--max-iter=1",
                           "--ribosum-weight=0", "--metrics-jsonl", str(metrics),
                           *extra, *flags, str(input_path)]
                result = subprocess.run(command, text=True, capture_output=True,
                                        timeout=20, check=True)
                lines = result.stdout.splitlines()
                structure = lines[lines.index(">SS_cons") + 1]
                expected = "(....)" if mode == "allowed" else "......"
                if structure != expected:
                    raise AssertionError(f"{label}/{mode}: {structure!r} != {expected!r}")
                events = [json.loads(line) for line in metrics.read_text().splitlines()]
                start = next(event for event in events if event["event"] == "run_start")
                summary = next(event for event in events if event["event"] == "run_summary")
                scores[mode] = summary["progressive_score"]
                if start["allow_isolated"] != (mode == "allowed"):
                    raise AssertionError(f"{label}/{mode}: switch missing from metrics")
                iterations = [event for event in events if event["event"] == "dd_iteration"]
                if "--max-iter=0" in extra:
                    if iterations or start["structure_decoder"] != "exact-coupled-IP":
                        raise AssertionError(f"{label}/{mode}: exact coupled IP not exercised")
                else:
                    if not iterations:
                        raise AssertionError(f"{label}/{mode}: progressive decoding was not exercised")
                    first_score = iterations[0]["x_beam"]
                    # The single-row profile can fall below the posterior
                    # cutoff even though the final two-row profile pairs.
                    if ((mode != "allowed" and first_score > 0) or
                        (mode == "allowed" and "--alifold" not in extra and first_score <= 0)):
                        raise AssertionError(f"{label}/{mode}: progressive score {first_score}")
                outputs[mode] = result.stdout
            if outputs["default"] != outputs["false"]:
                raise AssertionError(f"{label}: explicit false changes default prediction")
            if "--max-iter=0" in extra and scores["allowed"] <= scores["default"]:
                raise AssertionError(f"{label}: isolated-pair constraint missing from coupled IP")
            print(f"{label}: default, allowed, explicit false passed")
        check_stacked_aux(binary, root, args.with_ilp)


if __name__ == "__main__":
    main()
