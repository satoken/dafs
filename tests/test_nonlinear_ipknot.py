"""CONTRAfold and auxiliary BPPs use layered DD, including solver-free builds."""

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


binary = str(Path(sys.argv[1]).resolve())


def run(root, name, options, fasta):
    metrics = root / f"{name}.jsonl"
    result = subprocess.run(
        [binary, "-a", "CONTRAlign", "-s", "CONTRAfold", "--no-alifold",
         "--refinement=0", "--max-iter=12", "--ribosum-weight=0",
         *options, "--metrics-jsonl", str(metrics), str(fasta)],
        text=True, capture_output=True, check=True, timeout=20)
    events = [json.loads(line) for line in metrics.read_text().splitlines()]
    start = next(e for e in events if e["event"] == "run_start")
    assert start["structure_decoder"] == "IPknot-DD", start
    assert start["pseudoknot_backend"] == "DD", start
    iterations = [e for e in events if e["event"] == "dd_iteration"]
    assert iterations, result.stderr
    for event in iterations:
        for side in ["x", "y"]:
            bound, score = event[f"{side}_bound"], event[f"{side}_beam"]
            assert math.isfinite(bound) and bound + 1e-4 >= score, event
        assert math.isfinite(event["certified_track_lagrangian"]), event
    for event in events:
        if event["event"] == "dd_summary":
            assert event["stop_reason"] != "agreement", event
    lines = result.stdout.splitlines()
    structure = lines[lines.index(">SS_cons") + 1]
    print(f"{name}: {structure}", flush=True)
    return structure, events


with tempfile.TemporaryDirectory(prefix="dafs-nonlinear-ipknot-") as directory:
    root = Path(directory)
    fasta = root / "input.fa"
    fasta.write_text(">a\nGGGAAACCC\n>b\nGGGAAUCCC\n")
    for selector in ["--ipknot", "--fold-decoder=IPknot"]:
        for beam in [1, 100]:
            structure, events = run(
                root, f"{selector}-{beam}",
                [selector, "-g99,99", "-G99,99", "-w20",
                 f"--fold-dd-beam={beam}", f"--fold-final-beam={beam}"], fasta)
            assert "(" in structure, structure
            start = next(e for e in events if e["event"] == "run_start")
            assert start["folding_dd_beam"] == beam, start
            assert start["folding_final_beam"] == beam, start

    # Two crossing isolated pairs require distinct pseudoknot levels. Their
    # objective is known independently: 14 alignment edges plus four pair
    # terms across the two identical profiles. This also exercises the
    # previously ILP-only auxiliary-input path and outer certificate handling.
    fasta.write_text(">a\nAAAAAAAAAAAAAA\n>b\nAAAAAAAAAAAAAA\n")
    folding = root / "crossing.bp"
    folding.write_text("".join(f"> {i}\n1 10:1\n5 14:1\n14\n" for i in [1, 2]))
    alignment = root / "identity.mp"
    alignment.write_text("> 1 2\n" + "".join(f"{i} {i}:1\n" for i in range(1, 15)))
    for allowed in [True, False]:
        structure, events = run(
            root, f"crossing-{allowed}",
            ["--fold-decoder=IPknot", "--fold-aux", str(folding),
             "--align-aux", str(alignment), "--fold-th=0.2,0.2",
             "--fold-th1=0.2,0.2", "--align-th=0.01", "--weight=1",
             "--fourway-pct=0", "--fold-pct=0", "--align-pct=0",
             *([] if allowed else ["--disallow-isolated"])], fasta)
        expected = "(...[....)...]" if allowed else "." * 14
        assert structure == expected, (structure, expected)
        optimum = 14 * 0.99 + (4 * 0.8 if allowed else 0)
        iterations = [e for e in events if e["event"] == "dd_iteration"]
        assert all(e["certified_lagrangian"] + 1e-4 >= optimum for e in iterations)
        assert all(e["certified_track_lagrangian"] + 1e-4 >= optimum for e in iterations)
