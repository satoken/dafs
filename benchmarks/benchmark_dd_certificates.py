#!/usr/bin/env python3
"""Fixed-beam certificate ablation, using identical frozen code for controls."""
import argparse
import copy
import json
import shutil
import subprocess
import os
import random
import resource
import time
from statistics import fmean
from pathlib import Path
import benchmark_dd_recovery as previous
from benchmark_linearalifold_dd import ROOT, BASES, read, save, digest
from run_dd_diagnostics import CASES
from dd_exact_diagnostic import analyze

OUT=ROOT/'benchmarks/results-dd-certificates-20260920'
VARIANTS={'normal':[], 'unpruned':['--dd-unpruned-bound'],
          'block8':['--dd-block-bound=8'], 'block32':['--dd-block-bound=32'],
          'hybrid32':['--dd-unpruned-bound','--dd-block-bound=32']}
previous.OUT=OUT
previous.VARIANTS=VARIANTS

def prepare():
    if OUT.exists(): raise RuntimeError('Refusing overwrite')
    OUT.mkdir()
    for name in ('dafs-frozen','dafs-validated'):
        shutil.copy2(ROOT/'build-arm64/src/dafs',OUT/name)
    paths=[p for p in (ROOT/'src').rglob('*') if p.suffix in ('.cpp','.h','.hpp')]
    save(OUT/'source-hashes.json',{str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)})
    c=previous.config('pilot',[(b,v,n) for b in BASES for v in VARIANTS for n in (100,500)])
    print(c['experiment'])

def diagnostics():
    out=OUT/'diagnostics'
    if out.exists(): raise RuntimeError('Refusing overwrite')
    out.mkdir()
    original=read(previous.PREVIOUS/'alifold-config.json')
    templates={c['base_condition']:c for c in original['conditions'] if c['mode']=='normal' and c['max_iter']==100}
    results=[]
    for base in BASES:
        for name,x,y in CASES:
            fasta=out/(name+'.fa')
            if not fasta.exists():fasta.write_text(f'>one\n{x}\n>two\n{y}\n')
            for variant in VARIANTS:
                prefix=out/f'{base}-{name}-{variant}'
                args=[f'--max-iter=500' if a.startswith('--max-iter=') else a for a in templates[base]['args']]
                command=[str(OUT/'dafs-frozen'),*args,*VARIANTS[variant],
                         '--dd-diagnostics','--metrics-jsonl',str(prefix)+'.jsonl',str(fasta)]
                result=subprocess.run(command,capture_output=True,text=True,check=True)
                prefix.with_suffix('.prediction.fa').write_text(result.stdout)
                diagnostic=analyze(str(prefix)+'.jsonl',Path('/home/linuxbrew/.linuxbrew/lib/libglpk.so'))
                events=[json.loads(line) for line in Path(str(prefix)+'.jsonl').read_text().splitlines()]
                iterations=[e for e in events if e['event']=='dd_iteration']
                summary=next(e for e in events if e['event']=='dd_summary')
                results.append(dict(base=base,case=name,variant=variant,command=command,
                    diagnostic=diagnostic,summary=summary,terminal=iterations[-1],
                    negative_gap_iterations=sum(e['gap']<0 for e in iterations)))
    save(out/'report.json',results)
    print(json.dumps({'checked_runs':len(results)}))

def full():
    # Predeclared: compare strongest tested fixed-width certificate and the
    # essentially free no-pruning certificate at 500; also test hybrid at 100.
    c=previous.config('full',[(b,v,n) for b in BASES
        for v,n in [('normal',100),('normal',500),('unpruned',500),('hybrid32',100),('hybrid32',500)]])
    print(c['experiment'])

def timing(short=False):
    out=OUT/('serial-timing-short' if short else 'serial-timing')
    if out.exists():raise RuntimeError('Refusing overwrite')
    out.mkdir()
    config=read(OUT/'full-config.json')
    candidates=[c for c in config['conditions'] if c['max_iter']==500]
    datasets=(['murlet-dataset1-RF00031-2'] if short else
              ['murlet-dataset1-RF00001-2','murlet-dataset1-RF00059-2'])
    tasks=[(c,d,rep) for rep in range(3) for d in datasets for c in candidates]
    random.Random(42).shuffle(tasks)
    rows=[]
    env=dict(os.environ,OMP_NUM_THREADS='1')
    for c,d,rep in tasks:
        prefix=out/f"{c['id']}-{d}-{rep}"
        command=[config['binary'],*c['args'],'--metrics-jsonl',str(prefix)+'.jsonl',
                 str(ROOT/'benchmarks/data/murlet/inputs'/f'{d}.fa')]
        before=resource.getrusage(resource.RUSAGE_CHILDREN)
        start=time.perf_counter()
        result=subprocess.run(command,capture_output=True,text=True,check=True,env=env)
        wall=time.perf_counter()-start
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        events=[json.loads(line) for line in Path(str(prefix)+'.jsonl').read_text().splitlines()]
        summary=next(e for e in events if e['event']=='run_summary')
        prefix.with_suffix('.fa').write_text(result.stdout)
        rows.append(dict(base=c['base_condition'],variant=c['variant'],dataset=d,rep=rep,
            wall=wall,cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            dd_seconds=summary['dd_seconds'],iterations=summary['dd_iterations'],command=command))
    save(out/'report.json',rows)
    print(json.dumps({'serial_runs':len(rows)}))

def audit():
    import csv
    config=read(OUT/'full-config.json')
    rows=list(csv.DictReader((OUT/'full-cases.csv').open()))
    index={(r['condition'],r['variant'],int(r['budget']),r['dataset']):r for r in rows}
    comparisons=[]
    for base in BASES:
        for variant,budget,refbudget in [('unpruned',500,500),('hybrid32',500,500),('hybrid32',100,100),('hybrid32',100,500)]:
            candidates=[r for r in rows if r['condition']==base and r['variant']==variant and int(r['budget'])==budget]
            assert len(candidates)==85
            references=[index[base,'normal',refbudget,r['dataset']] for r in candidates]
            delta={k:fmean(float(c[k])-float(r[k]) for c,r in zip(candidates,references))
                   for k in ['mcc','sci','sps','dd_iterations','relative_gap','cpu_seconds','wall_seconds']}
            changed=sum(any(float(c[k])!=float(r[k]) for k in ['mcc','sci','sps']) for c,r in zip(candidates,references))
            comparisons.append(dict(base=base,variant=variant,budget=budget,reference_budget=refbudget,
                pairs=85,changed_quality_cases=changed,delta=delta))
    artifacts=[]
    conditions={c['id']:c for c in config['conditions']}
    for p in (OUT/'full/runs').glob('*/*/rep-*/result.json'):
        r=read(p);c=conditions[r['condition']]
        if c['variant']=='normal':continue
        ref=OUT/'full/runs'/f"{c['base_condition']}_normal_{c['max_iter']}"/r['dataset']/'rep-01/result.json'
        rr=read(ref)
        if digest(r['artifacts']['prediction'])!=digest(rr['artifacts']['prediction']):
            artifacts.append(dict(condition=c['id'],dataset=r['dataset']))
    report=dict(summary=read(OUT/'full-summary.json'),comparisons=comparisons,
                changed_prediction_files=artifacts,binary_sha256=digest(config['binary']))
    save(OUT/'report.json',report)
    print(json.dumps({'comparisons':len(comparisons),'changed_predictions':len(artifacts)}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','diagnostics','collect-pilot','full','collect-full','timing','timing-short','audit'])
    p.add_argument('--numeric',action='store_true');a=p.parse_args()
    if a.numeric:
        if a.action in ('full','collect-full','timing','timing-short','audit'):
            p.error('numeric follow-up currently supports pilot and diagnostics only')
        OUT=ROOT/'benchmarks/results-dd-certificates-numeric-20260920'
        VARIANTS={'normal':[], 'outward':['--dd-outward-lb'],
                  'unpruned_outward':['--dd-outward-lb','--dd-unpruned-bound'],
                  'hybrid32_outward':['--dd-outward-lb','--dd-unpruned-bound','--dd-block-bound=32']}
        previous.OUT=OUT;previous.VARIANTS=VARIANTS
    if a.action.startswith('collect-'):previous.collect(a.action[8:])
    elif a.action=='timing-short':timing(short=True)
    else:globals()[a.action]()
