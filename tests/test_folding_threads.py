"""CLI wiring and deterministic DD/final results for all folding DP backends."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
input_path = Path(sys.argv[2]).resolve()
openmp = sys.argv[3] == '1'
env = dict(os.environ, OMP_NUM_THREADS='4', OMP_DYNAMIC='FALSE')
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    common = ['--no-alifold', '--dynamic-cbp', '--max-iter=8', '--seed=42']
    for mode, model, align, decoder in [
        ('exact', 'CONTRAfold', 'CONTRAlign', 'Nussinov'),
        ('beam', 'lpc', 'LinearAlign', 'Nussinov'),
        ('ipknot', 'lpc', 'LinearAlign', 'IPknot'),
    ]:
        for isolated in ['--allow-isolated', '--disallow-isolated']:
            baseline = None
            for threads in ([1, 2, 4, 0] if openmp else [1, 0]):
                metrics = root / f'{mode}-{isolated}-{threads}.jsonl'
                extra = ['--dd-recovery-interval=1'] if mode == 'beam' else []
                if decoder == 'Nussinov': extra.append('--final-ribosum-weight=.1')
                command = [str(binary), *common, *extra, '-s', model, '-a', align,
                           '--fold-decoder', decoder, isolated,
                           f'--fold-dp-threads={threads}', '--metrics-jsonl',
                           str(metrics), str(input_path)]
                result = subprocess.run(command, env=env, capture_output=True,
                                        text=True, check=True, timeout=30)
                events = [json.loads(line) for line in metrics.read_text().splitlines()]
                start = next(e for e in events if e['event'] == 'run_start')
                assert start['folding_dp_threads'] == threads
                assert start['nussinov_backend'] == {
                    'exact': 'DAFS', 'beam': 'IPknot', 'ipknot': 'none'}[mode]
                deterministic = (result.stdout, [
                    {k: e[k] for k in ('event', 'progressive_score', 'ub', 'lb',
                                       'gap', 'iterations', 'stop_reason') if k in e}
                    for e in events if e['event'] in ('dd_summary', 'run_summary')])
                if baseline is None:
                    baseline = deterministic
                else:
                    assert baseline == deterministic, (mode, isolated, threads)
    for threads, message in [(-1, 'nonnegative')] + ([] if openmp else [(2, 'OpenMP')]):
        result = subprocess.run([str(binary), *common, f'--fold-dp-threads={threads}',
                                 str(input_path)], env=env, capture_output=True,
                                text=True, timeout=30)
        assert result.returncode != 0 and message in result.stdout + result.stderr
print('Folding thread CLI parity passed')
