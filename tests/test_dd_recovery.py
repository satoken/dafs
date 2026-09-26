#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'benchmarks'))
from dd_exact_diagnostic import exact_primal, relaxation_lp, exact_dual_at, analyze

def main(binary):
    glpk=Path('/home/linuxbrew/.linuxbrew/lib/libglpk.so')
    # Crossing pairs are legal in the left relaxation, not in the primal.
    p={'length_x':5,'length_y':5,'x_weight':1.,'y_weight':1.,
       'fold_threshold':0.,'align_threshold':0.,
       'px':[[0,3,.8],[1,4,.8]],'py':[[0,3,.8],[1,4,.8]],
       'pz':[], 'cbp':[[0,3,0,3,0.],[1,4,1,4,0.]]}
    optimal,count=exact_primal(p)
    assert abs(optimal-1.6)<1e-10 and count>1
    q={k:[] for k in ['qx','qy','qz']}
    assert abs(exact_dual_at(p,q,'')-1.6)<1e-10
    if glpk.exists():
        assert abs(relaxation_lp(p,glpk)-3.2)<1e-8
        assert abs(relaxation_lp(p,glpk,full=True)-1.6)<1e-8
    command=[str(binary),'-a','LinearAlign','-s','lpc','--alifold','--dynamic-cbp',
             '--ribosum-weight=0.3','--max-iter=100',str(ROOT/'tests/data/tiny.fa')]
    with tempfile.TemporaryDirectory(prefix='dafs-recovery-test-') as directory:
        directory=Path(directory)
        def run(flags,name):
            metric=directory/(name+'.jsonl')
            completed=subprocess.run(command+flags+['--metrics-jsonl',str(metric)],text=True,capture_output=True,check=True)
            return completed.stdout,[json.loads(l) for l in metric.read_text().splitlines()],metric
        baseline,events,_=run([],'baseline')
        explicit,explicit_events,_=run(['--dd-recovery-interval=0','--dd-beam-eta=0.5'],'explicit')
        assert baseline==explicit
        diagnostic,_,path=run(['--dd-diagnostics'],'diagnostic')
        assert diagnostic==baseline
        if glpk.exists():analyze(path,glpk)
        for flags in [['--dd-recovery-interval=1'],['--dd-beam-projected-norm'],
                      ['--dd-beam-eta=0.75'],['--dd-recovery-interval=10','--dd-beam-projected-norm']]:
            _,es,_=run(flags,'candidate')
            for e in es:
                if e['event'] in ['dd_iteration','dd_summary']:
                    assert e['gap'] >= -1e-4*max(1,abs(e['lb']))
        for flags in [['--dd-recovery-interval=-1'],['--dd-beam-eta=0'],
                      ['--dd-beam-eta=nan'],['--dd-projected-norm','--dd-beam-projected-norm'],
                      ['--dd-diagnostics'],['--dd-recovery-interval=1','-s','CONTRAfold']]:
            r=subprocess.run(command+flags,capture_output=True,text=True)
            assert r.returncode!=0,flags
    print('Recovery CLI, disabled compatibility, feasibility, diagnostic oracle tests passed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('binary',type=Path);a=p.parse_args();main(a.binary.resolve())
