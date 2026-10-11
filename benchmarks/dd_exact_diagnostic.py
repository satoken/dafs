#!/usr/bin/env python3
"""Offline tiny fixed-merge exact primal / left-row LP relaxation diagnostics.

Exponential enumeration is capped at length9, never used by DAFS production.
GLPK simplex solves the independent certificate's LP, not the integer model.
"""
import argparse
import ctypes as C
import json
from functools import lru_cache
from pathlib import Path

def matrix(values): return {(i,j):float(v) for i,j,v in values}

def coefficients(p):
    x={e:p['x_weight']*(v-p['fold_threshold']) for e,v in matrix(p['px']).items() if e[1]>=e[0]+2}
    y={e:p['y_weight']*(v-p['fold_threshold']) for e,v in matrix(p['py']).items() if e[1]>=e[0]+2}
    probs=matrix(p['pz'])
    z={(i,k):probs.get((i,k),0)-p['align_threshold'] for i in range(p['length_x']) for k in range(p['length_y'])}
    w={(i,j,k,l):v for i,j,k,l,v in p['cbp']}
    return x,y,z,w

def has_stacks(pairs):
    pairs=set(pairs)
    return all((i-1,j+1) in pairs or (i+1,j-1) in pairs for i,j in pairs)

def structures(scores,n,allow_isolated=True,mapping=None):
    @lru_cache(None)
    def configs(i,j):
        if i>=j:return ((),)
        out=list(configs(i+1,j))
        for k in range(i+2,j+1):
            if (i,k) not in scores:continue
            out.extend(((i,k),)+a+b for a in configs(i+1,k-1) for b in configs(k+1,j))
        return tuple(out)
    return {c:sum(scores[e] for e in c) for c in configs(0,n-1)
            if allow_isolated or (has_stacks(c) and
                (mapping is None or has_stacks((mapping[i],mapping[j]) for i,j in c)))}

def fold(n, scores, allow_isolated=True, mapping=None):
    if not allow_isolated:
        return max(structures(scores,n,False,mapping).values())
    @lru_cache(None)
    def dp(i,j):
        if i>=j:return 0.0
        best=dp(i+1,j)
        for k in range(i+1,j+1):
            if (i,k) in scores:
                best=max(best,scores[i,k]+dp(i+1,k-1)+dp(k+1,j))
        return best
    return dp(0,n-1)

def exact_primal(p):
    n,m=p['length_x'],p['length_y']
    if max(n,m)>9:raise ValueError('offline exact enumeration capped at length9')
    x,y,z,w=coefficients(p)
    useful={(i,k) for (i,k),v in z.items() if v>0}
    for i,j,k,l in w:useful.update(((i,k),(j,l)))
    forward={i:sorted(k for a,k in useful if a==i) for i in range(n)}
    mapping={};best=0.0;count=0
    def visit(i,last,alignment):
        nonlocal best,count
        if i==n:
            count+=1
            scores={(a,b):x[a,b]+y[k,l]+v for (a,b,k,l),v in w.items()
                    if mapping.get(a)==k and mapping.get(b)==l}
            best=max(best,alignment+fold(n,scores,p.get('allow_isolated',True),mapping));return
        visit(i+1,last,alignment)
        for k in forward[i]:
            if k<=last:continue
            mapping[i]=k;visit(i+1,k,alignment+z[i,k]);del mapping[i]
    visit(0,-1,0.0)
    return best,count

def relaxation_lp(p, library, full=False):
    if max(p['length_x'],p['length_y'])>9:
        raise ValueError('offline LP diagnostics capped at length9')
    lib=C.CDLL(str(library))
    def fn(name, restype, *args):
        f=getattr(lib,name);f.restype=restype;f.argtypes=list(args);return f
    ptr=C.c_void_p;I=C.c_int;D=C.c_double
    create=fn('glp_create_prob',ptr);delete=fn('glp_delete_prob',None,ptr)
    fn('glp_term_out',I,I)(0)
    addcols=fn('glp_add_cols',I,ptr,I);addrows=fn('glp_add_rows',I,ptr,I)
    colbound=fn('glp_set_col_bnds',None,ptr,I,I,D,D)
    rowbound=fn('glp_set_row_bnds',None,ptr,I,I,D,D)
    obj=fn('glp_set_obj_coef',None,ptr,I,D)
    load=fn('glp_load_matrix',None,ptr,I,C.POINTER(I),C.POINTER(I),C.POINTER(D))
    simplex=fn('glp_simplex',I,ptr,ptr)
    status=fn('glp_get_status',I,ptr);value=fn('glp_get_obj_val',D,ptr)
    x,y,z,w=coefficients(p)
    blocks={'x':x,'y':y,'z':z,'w':w}
    if full:
        # Convex hulls of exact noncrossing structures and exact monotone
        # alignments: optimum of the original DD relaxation on this tiny
        # support, as distinct from the cheaper left/row relaxation.
        useful={e for e,v in z.items() if v>0}
        for i,j,k,l in w:useful.update(((i,k),(j,l)))
        @lru_cache(None)
        def aligns(i,last):
            if i==p['length_x']:return ((),)
            out=list(aligns(i+1,last))
            for k in range(last+1,p['length_y']):
                if (i,k) in useful:
                    out.extend(((i,k),)+rest for rest in aligns(i+1,k))
            return tuple(out)
        blocks={'x':structures(x,p['length_x'],p.get('allow_isolated',True)),
                'y':structures(y,p['length_y'],p.get('allow_isolated',True)),
                'z':{c:sum(z[e] for e in c) for c in aligns(0,-1)},'w':w}
    index={(b,e):i+1 for i,(b,e) in enumerate((b,e) for b,d in blocks.items() for e in d)}
    constraints=[]
    if full:
        for b in ['x','y','z']:
            constraints.append((5,1.,1.,{index[b,c]:1. for c in blocks[b]}))
    else:
        for b,size in [('x',p['length_x']),('y',p['length_y']),('z',p['length_x'])]:
            for i in range(size):
                constraints.append((3,0.,1.,{index[b,e]:1. for e in blocks[b] if e[0]==i}))
    for b in ['x','y','z']:
        for e in {'x':x,'y':y,'z':z}[b]:
            row=({index[b,c]:1. for c in blocks[b] if e in c} if full else {index[b,e]:1.})
            for a,c,k,l in w:
                if (b=='x' and e==(a,c)) or (b=='y' and e==(k,l)) or (b=='z' and e in ((a,k),(c,l))):
                    row[index['w',(a,c,k,l)]]=-1.
            constraints.append((2 if b=='z' else 5,0.,0.,row))
    problem=create()
    try:
        fn('glp_set_obj_dir',None,ptr,I)(problem,2)
        addcols(problem,len(index));addrows(problem,len(constraints))
        for (b,e),i in index.items():colbound(problem,i,4,0.,1.);obj(problem,i,blocks[b][e])
        entries=[]
        for i,(kind,lo,hi,row) in enumerate(constraints,1):
            rowbound(problem,i,kind,lo,hi)
            entries.extend((i,j,v) for j,v in row.items())
        ia=(I*(len(entries)+1))(0,*[a for a,b,v in entries])
        ja=(I*(len(entries)+1))(0,*[b for a,b,v in entries])
        ar=(D*(len(entries)+1))(0.,*[v for a,b,v in entries])
        load(problem,len(entries),ia,ja,ar)
        if simplex(problem,None)!=0 or status(problem)!=5:raise RuntimeError('GLPK failed to find LP optimum')
        return value(problem)
    finally:delete(problem)

def exact_dual_at(p,state,prefix):
    x,y,z,w=coefficients(p)
    qx,qy,qz=[matrix(state[prefix+k]) for k in ['qx','qy','qz']]
    sx={e:v-qx.get(e,0.) for e,v in x.items()};sy={e:v-qy.get(e,0.) for e,v in y.items()}
    sz={e:v+qz.get(e,0.) for e,v in z.items()}
    @lru_cache(None)
    def align(i,k):
        if i==p['length_x'] or k==p['length_y']:return 0.
        return max(align(i+1,k),align(i,k+1),sz[i,k]+align(i+1,k+1))
    sw=sum(max(0.,v+qx.get((i,j),0.)+qy.get((k,l),0.)-qz.get((i,k),0.)-qz.get((j,l),0.))
           for (i,j,k,l),v in w.items())
    allowed=p.get('allow_isolated',True)
    return fold(p['length_x'],sx,allowed)+fold(p['length_y'],sy,allowed)+align(0,0)+sw

def analyze(path, library):
    events=[json.loads(line) for line in Path(path).read_text().splitlines()]
    results=[]
    for p in [e for e in events if e['event']=='dd_diagnostic_problem']:
        optimal,count=exact_primal(p);relaxed=relaxation_lp(p,library)
        assert relaxed>=optimal-1e-5
        states=[]
        for s in events:
            if s['event']!='dd_diagnostic_state' or s['merge_id']!=p['merge_id']:continue
            exact_beam=exact_dual_at(p,s,'beam_');exact_cert=exact_dual_at(p,s,'cert_')
            assert min(s['best_ub'],exact_beam,exact_cert)>=optimal-1e-4
            assert s['lb']<=optimal+1e-4
            states.append({'iteration':s['iteration'],'lb':s['lb'],'ub':s['best_ub'],
                'primal_shortfall':optimal-s['lb'],'left_row_lp_gap':relaxed-optimal,
                # Prefix UB also includes beam certificates, and can be lower
                # than this particular LP optimum. Do not call negative values
                # optimization error; report them as signed differences.
                'ub_minus_left_row_lp':s['best_ub']-relaxed,
                'exact_dual_at_beam_q':exact_beam,'exact_dual_at_cert_q':exact_cert})
        final=next(e for e in events if e['event']=='dd_summary' and e['merge_id']==p['merge_id'])
        assert final['lb']<=optimal+1e-4 and final['best_ub']>=optimal-1e-4
        results.append({'merge':p['merge_id'],'exact_primal':optimal,'left_row_lp':relaxed,
                        'alignments_enumerated':count,'checkpoints':states,
                        'final_iterations':final['iterations'],'stop_reason':final['stop_reason'],
                        'final_lb':final['lb'],'final_ub':final['best_ub'],
                        'primal_shortfall':optimal-final['lb'],
                        'lp_relaxation_gap':relaxed-optimal,
                        'final_ub_minus_lp':final['best_ub']-relaxed})
    if not results:raise ValueError('no eligible diagnostic merges')
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('metrics',type=Path)
    ap.add_argument('--glpk',type=Path,default=Path('/home/linuxbrew/.linuxbrew/lib/libglpk.so'))
    ap.add_argument('--output',type=Path);args=ap.parse_args()
    result=analyze(args.metrics,args.glpk)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text)
