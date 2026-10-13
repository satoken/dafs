#!/usr/bin/env python3
"""Strict rational audit of tiny fixed DD objectives, never a production oracle."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
from benchmark_dd_certificates import ROOT
from benchmark_linearalifold_dd import read,save

def exact_primal(p):
    n,m=p['length_x'],p['length_y']
    if max(n,m)>9:raise ValueError('tiny offline diagnostic only')
    px={(i,j):F(v) for i,j,v in p['px']};py={(i,j):F(v) for i,j,v in p['py']}
    pz={(i,j):F(v) for i,j,v in p['pz']}
    x={e:F(p['x_weight'])*(v-F(p['fold_threshold'])) for e,v in px.items()}
    y={e:F(p['y_weight'])*(v-F(p['fold_threshold'])) for e,v in py.items()}
    z={(i,k):pz.get((i,k),F(0))-F(p['align_threshold']) for i in range(n) for k in range(m)}
    w={(i,j,k,l):F(v) for i,j,k,l,v in p['cbp']}
    useful={e for e,v in z.items() if v>0}
    for i,j,k,l in w:useful.update(((i,k),(j,l)))
    forward={i:sorted(k for a,k in useful if a==i) for i in range(n)}
    mapping={};best=F(0)
    def visit(i,last,alignment):
        nonlocal best
        if i==n:
            scores={(a,b):x[a,b]+y[k,l]+v for (a,b,k,l),v in w.items()
                    if mapping.get(a)==k and mapping.get(b)==l}
            @lru_cache(None)
            def fold(a,b):
                if a>=b:return F(0)
                value=fold(a+1,b)
                for k in range(a+1,b+1):
                    if (a,k) in scores:value=max(value,scores[a,k]+fold(a+1,k-1)+fold(k+1,b))
                return value
            best=max(best,alignment+fold(0,n-1));return
        visit(i+1,last,alignment)
        for k in forward[i]:
            if k>last:
                mapping[i]=k;visit(i+1,k,alignment+z[i,k]);del mapping[i]
    visit(0,-1,F(0))
    return best

def main():
    root=ROOT/'benchmarks/results-dd-certificates-numeric-20260920'
    report=read(root/'diagnostics/report.json')
    cache={};checks=[]
    for r in report:
        key=r['base'],r['case']
        if key not in cache:
            path=Path(r['command'][r['command'].index('--metrics-jsonl')+1])
            p=next(json.loads(s) for s in path.read_text().splitlines()
                   if '"event":"dd_diagnostic_problem"' in s)
            cache[key]=exact_primal(p)
        f=cache[key];s=r['summary']
        checks.append(dict(base=r['base'],case=r['case'],variant=r['variant'],
            exact_primal=float(f),rational_primal=str(f),
            lb_valid=F(s['lb'])<=f,ub_valid=F(s['best_ub'])>=f,
            lb_excess=float(F(s['lb'])-f),ub_shortfall=float(f-F(s['best_ub']))))
    # The new downward LB must be conservative with ZERO tolerance.
    assert all(c['lb_valid'] for c in checks if c['variant']!='normal')
    save(root/'rational-audit.json',checks)
    print(json.dumps({'fixed_problems':len(cache),'checked_final_bounds':len(checks),
                      'outward_lb_strictly_valid':True}))

if __name__=='__main__':main()
