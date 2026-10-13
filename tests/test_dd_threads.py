"""DD subproblem workers preserve predictions, gradients and certificates."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
openmp = sys.argv[2] == '1'
ilp = sys.argv[3] == '1'
sanitizer_smoke = '--sanitizer-smoke' in sys.argv[4:]
env = dict(os.environ, OMP_NUM_THREADS='3', OMP_DYNAMIC='FALSE')


def deterministic(events):
    return [
        {k: v for k, v in event.items()
         if 'seconds' not in k and k != 'dd_workers'}
        for event in events
        if event['event'] in ('dd_iteration', 'dd_summary', 'run_summary')
    ]


with tempfile.TemporaryDirectory(prefix='dafs-dd-threads-') as directory:
    root = Path(directory)
    fasta = root / 'input.fa'
    fasta.write_text(
        '>a\nGGGAAACCC\n>b\nGGGAAUCCC\n>c\nGGGAAACCCC\n>d\nGGGUUUCCC\n'
        if sanitizer_smoke else
        '>a\nGGGAAACCCGGGUUUCCCGGGAAACCC\n'
        '>b\nGGGAAUCCCGGGUUUCCCGGGAAACCC\n'
        '>c\nGGGAAACCCGGGAAACCCGGGAAUCCC\n'
        '>d\nGGGAAACCCGGGUUCCCGGGAAACCC\n')
    common = ['--no-alifold', f'--max-iter={6 if sanitizer_smoke else 12}', '--seed=42']
    thread_counts = ([1, 3] if sanitizer_smoke else [1, 2, 3, 8, 0]) if openmp else [1, 0]
    cbp_modes = [['--dynamic-cbp']] if sanitizer_smoke else [[], ['--dynamic-cbp']]
    cases = [
        ('exact', 'CONTRAfold', 'CONTRAlign', 'Nussinov', []),
        ('beam', 'lpc', 'LinearAlign', 'Nussinov', []),
        ('exact-linear-align', 'CONTRAfold', 'LinearAlign', 'Nussinov', []),
        ('beam-exact-align', 'lpc', 'CONTRAlign', 'Nussinov', []),
        ('beam-dense', 'lpc', 'LinearAlign', 'Nussinov', ['--dense-lagrangian']),
        ('recovery', 'lpc', 'LinearAlign', 'Nussinov',
         ['--dd-recovery-interval=1', '--dd-tight-alignment',
          '--dd-column-bound', '--dd-outward-lb', '--dd-unpruned-bound']),
        ('ipknot', 'lpc', 'LinearAlign', 'IPknot', []),
    ]
    if ilp:
        cases.append(('ipknot-ilp', 'CONTRAfold', 'CONTRAlign', 'IPknot', []))
    for mode, model, align, decoder, extra in cases:
        print(f'Checking {mode}', flush=True)
        for isolated in ['--allow-isolated', '--disallow-isolated']:
            for dynamic in cbp_modes:
                baseline = None
                for threads in thread_counts:
                    metrics = root / 'metrics.jsonl'
                    command = [str(binary), *common, *extra, *dynamic,
                               '-s', model, '-a', align, '--fold-decoder', decoder,
                               isolated, f'--dd-threads={threads}',
                               '--metrics-jsonl', str(metrics), str(fasta)]
                    result = subprocess.run(command, env=env, capture_output=True,
                                            text=True, check=True, timeout=30)
                    events = [json.loads(line) for line in metrics.read_text().splitlines()]
                    start = next(e for e in events if e['event'] == 'run_start')
                    assert start['dd_threads'] == threads
                    summaries = [e for e in events if e['event'] == 'dd_summary']
                    assert len(summaries) == 3
                    expected = min(3, threads or 3) if openmp else 1
                    assert all(e['dd_workers'] == expected for e in summaries)
                    actual = result.stdout, deterministic(events)
                    if baseline is None:
                        baseline = actual
                    else:
                        assert actual == baseline, (mode, isolated, dynamic, threads)

    # Within-DP requests take precedence, including the runtime setting 0.
    for dp_threads in ([0, 2] if openmp else [0]):
        metrics = root / 'precedence.jsonl'
        command = [str(binary), *common, '-s', 'lpc', '-a', 'LinearAlign',
                   f'--fold-dp-threads={dp_threads}', '--metrics-jsonl',
                   str(metrics), str(fasta)]
        serial = subprocess.run(command, env=env, capture_output=True,
                                text=True, check=True, timeout=30)
        baseline = deterministic([json.loads(line) for line in metrics.read_text().splitlines()])
        parallel = subprocess.run([*command, f'--dd-threads={3 if openmp else 0}'],
                                  env=env, capture_output=True, text=True,
                                  check=True, timeout=30)
        events = [json.loads(line) for line in metrics.read_text().splitlines()]
        assert parallel.stdout == serial.stdout and deterministic(events) == baseline
        assert all(e['dd_workers'] == 1 for e in events if e['event'] == 'dd_summary')

    # A runtime request can serialize, and OpenMP may supply fewer workers
    # than explicitly requested. Metrics report the team that actually ran.
    if openmp:
        for overrides in [dict(OMP_NUM_THREADS='1'), dict(OMP_THREAD_LIMIT='1')]:
            metrics = root / 'runtime.jsonl'
            result = subprocess.run(
                [str(binary), *common, '-s', 'lpc', '-a', 'LinearAlign',
                 '--dd-threads=0', '--metrics-jsonl', str(metrics), str(fasta)],
                env=dict(env, **overrides), capture_output=True, text=True,
                check=True, timeout=30)
            events = [json.loads(line) for line in metrics.read_text().splitlines()]
            assert all(e['dd_workers'] == 1 for e in events if e['event'] == 'dd_summary')

    for threads, message in [(-1, 'nonnegative')] + ([] if openmp else [(2, 'OpenMP')]):
        result = subprocess.run([str(binary), *common, f'--dd-threads={threads}',
                                 str(fasta)], env=env, capture_output=True,
                                text=True, timeout=30)
        assert result.returncode != 0 and message in result.stdout + result.stderr

print('DD subproblem thread CLI parity passed')
