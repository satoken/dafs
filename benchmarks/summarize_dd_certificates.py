#!/usr/bin/env python3
"""Export tables from validated certificate experiments without hand transcription."""
import json
from collections import defaultdict
from statistics import fmean,median
from benchmark_linearalifold_dd import ROOT,BASES,read,save,digest

def main():
    out=ROOT/'benchmarks/results-dd-certificates-20260920'
    numeric=ROOT/'benchmarks/results-dd-certificates-numeric-20260920'
    full=read(out/'full-summary.json');follow=read(numeric/'pilot-summary.json')
    report=read(out/'report.json');timing=read(out/'serial-timing/report.json')
    short_timing=read(out/'serial-timing-short/report.json')
    labels=dict(zip(BASES,['A','B','C','D']))
    lines=['# DD証明上界改善：集計表','',
           'A=lpc×LinearAlign、B=lpc×LinearAlign-CONTRAlign、C=lpv×LinearAlign、D=lpv×LinearAlign-ProbConsRNA。',
           'normal=従来、unpruned=枝刈りなし上界、hybrid32=枝刈りなし+block32。全てalifoldあり・beam100。','',
           '## Murlet dataset1 全85ケース','',
           '各行85ケースの平均。DD反復は1ケースの9 mergeの合計。違反数は各mergeの最終反復におけるbeam trackの制約違反数の合計。gapは各mergeの (UB−LB)/max(1,|LB|) の平均。停止数は1ケース平均（最大9）。bulk wallには併走の影響がある。','',
           '|構成|方式|max-iter|MCC|SCI|SPS|wall秒|CPU秒|DD秒|DD反復|違反数|gap %|証明停止数|',
           '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in full:
        assert r['n']==85
        values=[labels[r['condition']],r['variant'],str(r['budget'])]
        values += [f"{r[k]:.5f}" for k in ('mcc','sci','sps')]
        values += [f"{r[k]:.3f}" for k in ('wall_seconds','cpu_seconds','dd_wall_seconds')]
        values += [f"{r[k]:.1f}" for k in ('dd_iterations','dd_violations_total')]
        values += [f"{100*r['relative_gap']:.3f}",f"{r['certified_stops']:.3f}"]
        lines.append('|'+ '|'.join(values)+'|')
    lines += ['', '## 同じケースへの対応差','',
              '|構成|方式|max-iter|比較先normal|ΔMCC|ΔSCI|ΔSPS|ΔDD反復|品質指標が変わったケース数 /85|',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in report['comparisons']:
        d=r['delta']
        lines.append(f"|{labels[r['base']]}|{r['variant']}|{r['budget']}|{r['reference_budget']}|{d['mcc']:+.5f}|{d['sci']:+.5f}|{d['sps']:+.5f}|{d['dd_iterations']:+.1f}|{r['changed_quality_cases']}|")
    lines += ['', '## 逐次再計測（72実行）','',
              '他の今回のベンチマーク終了後、2ケース（RF00001-2 / RF00059-2）×4構成×3方式×3反復をseed42でシャッフルして逐次実行。SCIの外部評価時間は含めない。各ケースで3回の中央値を取り、その2ケース平均を示す。85ケース全体の速度比ではない。','',
              '|構成|方式|wall秒|CPU秒|DD秒|DD反復|対normal wall比|',
              '|---|---|---:|---:|---:|---:|---:|']
    groups=defaultdict(list)
    for r in timing:groups[r['base'],r['variant'],r['dataset']].append(r)
    aggregated={}
    for b in BASES:
        for v in ('normal','unpruned','hybrid32'):
            rs=[rs for (base,variant,_),rs in groups.items() if base==b and variant==v]
            assert len(rs)==2 and all(len(s)==3 for s in rs)
            aggregated[b,v]={k:fmean(median(x[k] for x in s) for s in rs)
                              for k in ('wall','cpu','dd_seconds','iterations')}
            s=aggregated[b,v];ratio=s['wall']/aggregated[b,'normal']['wall']
            lines.append(f"|{labels[b]}|{v}|{s['wall']:.3f}|{s['cpu']:.3f}|{s['dd_seconds']:.3f}|{s['iterations']:.1f}|{ratio:.3f}|")
    lines += ['', '### 短い配列の追加逐次計測（36実行）','',
              'RF00031-2（入力長61〜73）：枝刈りなし判定の効く長さ域を追加で確認。各方式3回の中央値。上の2ケースと混ぜずに示す。','',
              '|構成|方式|wall秒|DD秒|DD反復|対normal wall比|',
              '|---|---|---:|---:|---:|---:|']
    for b in BASES:
        control=median(r['wall'] for r in short_timing if r['base']==b and r['variant']=='normal')
        for v in ('normal','unpruned','hybrid32'):
            rs=[r for r in short_timing if r['base']==b and r['variant']==v]
            assert len(rs)==3
            wall=median(r['wall'] for r in rs)
            dd=median(r['dd_seconds'] for r in rs);it=median(r['iterations'] for r in rs)
            lines.append(f"|{labels[b]}|{v}|{wall:.3f}|{dd:.3f}|{it:.1f}|{wall/control:.3f}|")
    lines += ['', '## 下向き再採点の予備評価（12ケース、max-iter=500）','',
              '別の凍結バイナリ内で比較。outwardは下界再採点、unpruned_outward/hybrid32_outwardは上界改善との併用。下界変更がPolyak更新へ伝わるため、上界のみの評価とは分ける。','',
              '|構成|方式|MCC|SCI|SPS|DD反復|gap %|',
              '|---|---|---:|---:|---:|---:|---:|']
    for r in follow:
        if r['budget']!=500:continue
        assert r['n']==12
        lines.append(f"|{labels[r['condition']]}|{r['variant']}|{r['mcc']:.5f}|{r['sci']:.5f}|{r['sps']:.5f}|{r['dd_iterations']:.1f}|{100*r['relative_gap']:.3f}|")
    lines+=['','## 検査数','',
            f"Murlet：予備480 + 全体{sum(r['n'] for r in full)} + 下界予備{sum(r['n'] for r in follow)} = {480+sum(r['n'] for r in full)+sum(r['n'] for r in follow)}実行。全てsuccess、metrics valid、MCC/SCI/SPS評価完備をcollect処理で検査。",'',
            '短い固定問題：120+96=216実行。うち下界追加の24問題・96最終上下界は有理数による整数最適値とも比較。別途逐次時間計測72+36=108実行。','']
    path=ROOT/'docs/dafs_dd_certificates_tables_2026-09-20_ja.md'
    path.write_text('\n'.join(lines))
    changed=[p for p,s in read(numeric/'source-hashes.json').items() if digest(ROOT/p)!=s]
    assert not changed,changed
    save(out/'final-audit.json',dict(murlet_valid_runs=480+sum(r['n'] for r in full)+sum(r['n'] for r in follow),
        diagnostic_runs=216,serial_runs=len(timing)+len(short_timing),numeric_binary_sha256=digest(numeric/'dafs-frozen'),
        upper_binary_sha256=digest(out/'dafs-frozen'),source_changes_since_numeric_freeze=changed,
        tables_sha256=digest(path),serial_summary={f'{b}:{v}':s for (b,v),s in aggregated.items()}))
    print(path)

if __name__=='__main__':main()
