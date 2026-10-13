#!/usr/bin/env python3
"""Alifold DD recovery ablations; immutable binary/configs and paired results."""
from __future__ import annotations
import argparse
import copy
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path
from statistics import fmean
from benchmark_linearalifold_dd import BASES, ROOT, digest, read, save
from benchmark_dd_budget import PILOT_DATASET_IDS, write_csv
from benchmark_luna_optimization import trace_metrics

OUT = ROOT / 'benchmarks/results-dd-recovery-20260920'
PREVIOUS = ROOT / 'benchmarks/results-linearalifold-dd-convergence-20260919'
VARIANTS = {
    'normal': [],
    'recovery': ['--dd-recovery-interval=10'],
    'beam_projected': ['--dd-beam-projected-norm'],
    'beam_eta025': ['--dd-beam-eta=0.25'],
    'beam_eta075': ['--dd-beam-eta=0.75'],
    'beam_eta100': ['--dd-beam-eta=1'],
    'recovery_projected': ['--dd-recovery-interval=10', '--dd-beam-projected-norm'],
}

def config(scope, selections):
    original = read(PREVIOUS/'alifold-config.json')
    result = {k:copy.deepcopy(original[k]) for k in
              ('timeout_seconds','repetitions','random_seed','environment','sci','dataset_manifests')}
    pilot=scope.endswith('pilot')
    binary=OUT/('dafs-frozen' if scope=='pilot' else 'dafs-validated')
    result.update(binary=str(binary), output_dir=str(OUT/scope), conditions=[],
                  dataset_filter={'ids':list(PILOT_DATASET_IDS)} if pilot else {'collections':['dataset1']})
    templates = {c['base_condition']:c for c in original['conditions'] if c['mode']=='normal' and c['max_iter']==100}
    for base, variant, budget in selections:
        c = copy.deepcopy(templates[base])
        c['args'] = [f'--max-iter={budget}' if a.startswith('--max-iter=') else a for a in c['args']]
        c['args'] += (['--dd-recovery-interval=1'] if variant=='recovery_dense' else VARIANTS[variant])
        c.update(id=f'{base}_{variant}_{budget}', variant=variant, max_iter=budget)
        c.pop('mode',None)
        result['conditions'].append(c)
    result['experiment'] = {'binary_sha256':digest(binary), 'scope':scope,
                            'beam':100, 'certified_polyak_scale':0.5,
                            'expected_runs':len(selections)*(12 if pilot else 85)}
    path = OUT/f'{scope}-config.json'
    if path.exists(): raise RuntimeError(f'Refusing overwrite: {path}')
    save(path,result)
    return result

def prepare():
    if OUT.exists(): raise RuntimeError(f'Refusing overwrite: {OUT}')
    OUT.mkdir(parents=True)
    shutil.copy2(ROOT/'build-arm64/src/dafs',OUT/'dafs-frozen')
    paths = [p for p in (ROOT/'src').rglob('*') if p.suffix in ('.cpp','.h','.hpp')]
    save(OUT/'source-hashes.json',{str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)})
    c = config('pilot',[(b,v,n) for b in BASES for v in VARIANTS
                        for n in ((100,200,300,500) if v=='normal' else (100,200))])
    print(c['experiment'])

def collect(scope):
    c = read(OUT/f'{scope}-config.json')
    assert read(OUT/scope/'manifest.json')['binary_sha256']==c['experiment']['binary_sha256']==digest(c['binary'])
    conditions = {x['id']:x for x in c['conditions']}
    rows=[]
    for path in sorted((OUT/scope/'runs').glob('*/*/rep-*/result.json')):
        r=read(path); condition=conditions[r['condition']]
        assert r['status']=='success' and r['metrics_validation']['valid'] and r['quality_complete'],path
        summaries=[]; last={}; run={}
        with (path.parent/'metrics.jsonl').open() as stream:
            for line in stream:
                if '"event":"dd_iteration"' in line:
                    # Keep last line per merge; parse only 9 terminal events.
                    key=re.search(r'"merge_id":(\d+)',line).group(1)
                    last[key]=line
                elif '"event":"dd_summary"' in line: summaries.append(json.loads(line))
                elif '"event":"run_summary"' in line:run=json.loads(line)
        metrics=trace_metrics([*summaries,*[json.loads(v) for v in last.values()],run],r)
        rows.append({'condition':condition['base_condition'],'variant':condition['variant'],
            'budget':condition['max_iter'],'dataset':r['dataset'],**metrics,
            'relative_gap':fmean(s['gap']/max(1,abs(s['lb'])) for s in summaries),
            'recovery_attempts':sum(s.get('recovery_attempts',0) for s in summaries),
            'recovery_improvements':sum(s.get('recovery_improvements',0) for s in summaries),
            'recovery_seconds':sum(s.get('recovery_seconds',0) for s in summaries),
            'certified_stops':sum(s['stop_reason']=='certified_gap' for s in summaries),
            'limit_stops':sum(s['stop_reason']=='max_iterations' for s in summaries)})
    assert len(rows)==c['experiment']['expected_runs'],len(rows)
    assert len({(r['condition'],r['variant'],r['budget'],r['dataset']) for r in rows})==len(rows)
    groups=defaultdict(list)
    for r in rows:groups[r['condition'],r['variant'],r['budget']].append(r)
    summary=[]
    for (base,v,n),rs in sorted(groups.items()):
        assert len(rs)==(12 if scope.endswith('pilot') else 85)
        s={'condition':base,'variant':v,'budget':n,'n':len(rs)}
        for k in ['mcc','sci','sps','wall_seconds','cpu_seconds','dd_wall_seconds','dd_iterations',
                  'dd_violations_total','certified_violations_total','relative_gap',
                  'recovery_attempts','recovery_improvements','recovery_seconds','certified_stops','limit_stops']:
            s[k]=fmean(r[k] for r in rs)
        summary.append(s)
    save(OUT/f'{scope}-summary.json',summary)
    write_csv(OUT/f'{scope}-summary.csv',summary)
    write_csv(OUT/f'{scope}-cases.csv',rows)
    print(json.dumps({'scope':scope,'valid_runs':len(rows),'groups':len(summary)}))
    return summary

def final_audit():
    import csv
    from benchmark_linearalifold_dd import paired_summary
    summary=read(OUT/'full-summary.json')+read(OUT/'control-full-summary.json')
    rows=list(csv.DictReader((OUT/'full-cases.csv').open()))+list(csv.DictReader((OUT/'control-full-cases.csv').open()))
    for r in rows:
        for k,v in list(r.items()):
            if k in ('condition','variant','dataset','dd_stop_reasons'):continue
            r[k]=float(v) if v else None
        r['budget']=int(r['budget']);r['repetition']=1
        # collect() has asserted these for every raw result before writing CSV.
        r.update(status='success',metrics_valid=True,quality_complete=True)
    groups=defaultdict(list)
    for r in rows:groups[r['condition'],r['variant'],r['budget']].append(r)
    comparisons=[]
    for (b,v,n),rs in groups.items():
        if v=='normal':continue
        comparisons.append(paired_summary(rs,groups[b,'normal',500],f'{b}:{v}:{n}:vs-normal500'))
        if n==100:
            comparisons.append(paired_summary(rs,groups[b,'normal',100],f'{b}:{v}:{n}:vs-normal100'))
        elif n==300:
            comparisons.append(paired_summary(rs,groups[b,'normal',300],f'{b}:{v}:{n}:vs-normal300'))
    old={(r['condition'],int(r['budget']),r['dataset']):r for r in
         csv.DictReader((PREVIOUS/'case-convergence.csv').open())
         if r['scope']=='alifold' and r['mode']=='normal'}
    drift=[]
    for r in rows:
        if r['variant']!='normal':continue
        o=old[r['condition'],r['budget'],r['dataset']]
        delta={k:r[k]-float(o[k]) for k in ['mcc','sci','sps']}
        if any(abs(v)>1e-12 for v in delta.values()):
            drift.append({'base':r['condition'],'budget':r['budget'],'dataset':r['dataset'],'delta':delta})
    changed=[name for name,sha in read(OUT/'validated-source-hashes.json').items()
             if digest(ROOT/name)!=sha]
    assert not changed,changed
    dense=read(OUT/'dense-pilot-summary.json')
    dense_count=sum(r['n'] for r in dense)
    assert dense_count==96
    audit={'pilot_runs':768,'full_runs':len(rows),'dense_pilot_runs':dense_count,
           'valid_runs':768+len(rows)+dense_count,'dense_pilot_summary':dense,
           'diagnostic_runs':len(read(OUT/'exact-diagnostics/report.json')),
           'summary':summary,'comparisons':comparisons,'normal_quality_drift_from_old_binary':drift,
           'validated_binary_sha256':digest(OUT/'dafs-validated'),
           'source_changes_since_full_freeze':changed,
           'analysis_hashes':{str(p.relative_to(ROOT)):digest(p) for p in
               [Path(__file__),ROOT/'benchmarks/dd_exact_diagnostic.py',ROOT/'benchmarks/run_dd_diagnostics.py']}}
    save(OUT/'report.json',audit)
    print(json.dumps({'valid_runs':audit['valid_runs'],'comparisons':len(comparisons),'old_binary_drift_cases':len(drift)}))

def select_full():
    if (OUT/'dafs-validated').exists():raise RuntimeError('Refusing to overwrite validated binary')
    shutil.copy2(ROOT/'build-arm64/src/dafs',OUT/'dafs-validated')
    paths=[p for p in (ROOT/'src').rglob('*') if p.suffix in ('.cpp','.h','.hpp')]
    save(OUT/'validated-source-hashes.json',{str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)})
    summary=read(OUT/'pilot-summary.json'); selections=[]; rationale=[]
    for b in BASES:
        reference=next(r for r in summary if r['condition']==b and r['variant']=='normal' and r['budget']==500)
        candidates=[r for r in summary if r['condition']==b and r['variant']!='normal' and r['budget']==200]
        # Predeclared compromise criterion: minimize largest mean quality
        # deficit vs normal500. Not a statistical equivalence test.
        def key(r):
            deficits=[max(0,reference[k]-r[k]) for k in ['mcc','sci','sps']]
            return max(deficits),sum(deficits),r['wall_seconds'],r['variant']
        best=min(candidates,key=key)
        rationale.append({'condition':b,'variant':best['variant'],'criterion':key(best)})
        selections.extend((b,'normal',n) for n in (100,500))
        selections.extend((b,best['variant'],n) for n in (100,200,300))
        # Always retain recovery-only ablation at100, regardless of selection.
        if best['variant']!='recovery': selections.append((b,'recovery',100))
    save(OUT/'selection.json',rationale)
    c=config('full',selections)
    print(json.dumps({'selection':rationale,**c['experiment']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','pilot','select','full','audit','dense-prepare','dense-pilot','control-prepare','control-full'])
    a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='select':select_full()
    elif a.command=='audit':final_audit()
    elif a.command=='dense-prepare':
        c=config('dense-pilot',[(b,'recovery_dense',n) for b in BASES for n in (100,200)])
        print(c['experiment'])
    elif a.command=='control-prepare':
        c=config('control-full',[(b,'normal',300) for b in BASES])
        print(c['experiment'])
    else:collect(a.command)
