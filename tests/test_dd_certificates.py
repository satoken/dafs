import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    fasta=root/'mixed.fa'
    fasta.write_text('>one\nGCGCGCGCG\n>two\nCGCGCGCGC\n')
    common=['-a','LinearAlign-ProbConsRNA','-s','lpv','--alifold','--dynamic-cbp',
            '--ribosum-weight=0.3','--final-ribosum-weight=0','--weight=2',
            '--max-iter=500','--seed=42','--refinement=0']
    def run(name,extra=()):
        metric=root/(name+'.jsonl')
        result=subprocess.run([str(binary),*common,*extra,'--metrics-jsonl',str(metric),str(fasta)],
                              capture_output=True,text=True,check=True)
        events=[json.loads(line) for line in metric.read_text().splitlines()]
        return result.stdout,next(e for e in events if e['event']=='dd_summary')
    control,c=run('normal')
    off,o=run('off',['--dd-block-bound=0'])
    assert control==off and c['iterations']==o['iterations']
    for name,args in [('unpruned',['--dd-unpruned-bound']),('block',['--dd-block-bound=32']),
                      ('both',['--dd-unpruned-bound','--dd-block-bound=32'])]:
        prediction,s=run(name,args)
        assert prediction==control
        assert s['iterations']<c['iterations'] and s['stop_reason']=='certified_gap'
        assert 0<=s['gap']<=1e-4*max(1,abs(s['lb']))
    for args in [['--dd-block-bound=-1'],['--dd-block-bound=65'],
                 ['--dd-block-bound=8','--dense-lagrangian'],
                 ['--dd-unpruned-bound','--align-th=-0.1'],
                 ['--dd-unpruned-bound','-a','ProbCons']]:
        result=subprocess.run([str(binary),*common,*args,str(fasta)],capture_output=True,text=True)
        assert result.returncode!=0,args
    fasta.write_text('>one\nGCGAAACGC\n>two\nGGGAAACCC\n')
    common=['-a','LinearAlign','-s','lpc','--alifold','--dynamic-cbp',
            '--ribosum-weight=0','--final-ribosum-weight=0','--align-th=0.005',
            '--max-iter=500','--seed=42','--refinement=0']
    control,c=run('rounding-normal')
    prediction,s=run('rounding-outward',['--dd-outward-lb'])
    assert c['gap']<0 and c['iterations']==500
    assert prediction==control and s['stop_reason']=='certified_gap' and s['iterations']<500
    assert 0<=s['gap']<=1e-4*max(1,abs(s['lb']))
