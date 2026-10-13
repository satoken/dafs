#!/usr/bin/env python3
"""Run tiny-RNA fixed-merge diagnostics separately from Murlet evaluation."""
import argparse
import json
import subprocess
from pathlib import Path
from benchmark_dd_recovery import OUT, PREVIOUS, BASES, save
from dd_exact_diagnostic import analyze, relaxation_lp

CASES=[('stem_variant','GGGAAACCC','GGGAAUCCC'),
       ('stem_length','GGGAAACCC','GGAAACCC'),
       ('stem_shift','GCGAAACGC','GGGAAACCC'),
       ('opposite','CCCAAAGGG','GGGAAACCC'),
       ('mixed','GCGCGCGCG','CGCGCGCGC'),
       ('unstructured','AAAAAAAAA','UUUUUUUUU')]

def main():
    out=OUT/'exact-diagnostics'
    if out.exists():raise RuntimeError('Refusing to overwrite diagnostics')
    out.mkdir()
    original=json.loads((PREVIOUS/'alifold-config.json').read_text())
    templates={c['base_condition']:c for c in original['conditions'] if c['mode']=='normal' and c['max_iter']==100}
    results=[]
    for base in BASES:
        for name,x,y in CASES:
            fasta=out/(name+'.fa')
            if not fasta.exists():fasta.write_text(f'>one\n{x}\n>two\n{y}\n')
            for mode,budget in [('normal',100),('normal',500),('recovery',100)]:
                prefix=out/f'{base}-{name}-{mode}-{budget}'
                args=[f'--max-iter={budget}' if a.startswith('--max-iter=') else a for a in templates[base]['args']]
                if mode=='recovery':args.append('--dd-recovery-interval=10')
                command=[str(OUT/'dafs-validated'),*args,'--dd-diagnostics','--metrics-jsonl',str(prefix)+'.jsonl',str(fasta)]
                result=subprocess.run(command,capture_output=True,text=True,check=True)
                prefix.with_suffix('.prediction.fa').write_text(result.stdout)
                diagnostic=analyze(str(prefix)+'.jsonl',Path('/home/linuxbrew/.linuxbrew/lib/libglpk.so'))
                results.append({'base':base,'case':name,'mode':mode,'budget':budget,
                                'command':command,'diagnostic':diagnostic})
    save(out/'report.json',results)
    print(json.dumps({'diagnostic_runs':len(results),'all_exact_bounds_checked':True}))

def strong():
    report=json.loads((OUT/'exact-diagnostics/report.json').read_text())
    results=[]
    for r in report:
        if r['mode']!='normal' or r['budget']!=500:continue
        metric=Path(r['command'][r['command'].index('--metrics-jsonl')+1])
        events=[json.loads(line) for line in metric.read_text().splitlines()]
        p=next(e for e in events if e['event']=='dd_diagnostic_problem')
        d=r['diagnostic'][0]
        lp=relaxation_lp(p,Path('/home/linuxbrew/.linuxbrew/lib/libglpk.so'),full=True)
        assert d['exact_primal']-1e-5<=lp<=d['left_row_lp']+1e-5
        results.append({'base':r['base'],'case':r['case'],'primal':d['exact_primal'],
            'full_dd_lp':lp,'left_row_lp':d['left_row_lp'],'final_lb':d['final_lb'],
            'final_ub':d['final_ub'],'iterations':d['final_iterations']})
    save(OUT/'exact-diagnostics/strong-lp.json',results)
    print(json.dumps({'strong_lp_problems':len(results)}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--strong-only',action='store_true');a=p.parse_args()
    strong() if a.strong_only else main()
