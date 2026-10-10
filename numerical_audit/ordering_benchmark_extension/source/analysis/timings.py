#!/usr/bin/env python3
"""Report recorded classical work without conflating measurement intervals."""
import argparse,csv,hashlib,json,statistics
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 summary=json.loads((a.output/'summary.json').read_text());cases=summary['cases'];rows=[];v2=[];inputs={}
 for case in cases:
  f=a.root/'optimized/designs'/(case+'.design.json');d=json.loads(f.read_text());inputs[str(f.relative_to(a.root))]=sha(f)
  for name,r in d['optimizer_reports'].items():
   rows.append({'case':case,'ordering':name,'units':r['units'],'nonzero_pairs':r['nonzero_pair_count'],'pair_setup_seconds':r['pair_setup_seconds'],'search_seconds':r['search_seconds'],'evaluations':r['evaluations'],'rescored':r['rescored'],'accepted_moves_both_starts':sum(len(s['history']) for s in r['starts']),'selected_start_accepted_moves':len(r['accepted_moves']),'cache_group':'pauli' if 'pauli' in name else 'conserving'})
  f=a.root/'published/v2/designs'/(case+'.json');d=json.loads(f.read_text());inputs[str(f.relative_to(a.root))]=sha(f)
  meta=d['methods']['tranter_v2_l1_first_order']['metadata'];v2.append({'case':case,'latest_native_seconds':meta.get('latest_native_seconds'),'metadata':{k:v for k,v in meta.items() if k not in ('history','order','native_history','eq4_order','chronological_order')}})
 scopes={}
 for name in ('fixed','optimized','published'):
  f=a.root/name/'timing_scope.json';scopes[name]=json.loads(f.read_text());inputs[str(f.relative_to(a.root))]=sha(f)
 with (a.output/'optimizer_runtime.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 stats={}
 for name in sorted({r['ordering'] for r in rows}):
  rr=[r for r in rows if r['ordering']==name];ss=[r['search_seconds'] for r in rr];stats[name]={'min_search_seconds':min(ss),'median_search_seconds':statistics.median(ss),'max_search_seconds':max(ss),'total_examined_positions':sum(r['evaluations'] for r in rr)}
 result={'status':'complete','interpretation':'descriptive timings under concurrent research jobs; not isolated hardware speedup benchmarks','optimizer_rows':rows,'optimizer_summary':stats,'v2_preparation':v2,'timing_scopes':scopes,'inputs_sha256':inputs,'source_sha256':sha(Path(__file__))}
 (a.output/'timings.json').write_text(json.dumps(result,indent=2)+'\n')
 lines=['# Recorded runtimes and classical work','','These are observed times from concurrent single-threaded research jobs, not an isolated speed comparison. See environment.json for hardware and software. All raw evolution rows retain their own seconds field.','','## Optimizer setup and search','','Pair-cache setup and search are separate. Conserving JD and J0 share one cache; the two Pauli arms share another. Repeated setup values in the CSV describe this shared cache and must not be summed as independent measurements. Search includes between-sweep checkpoint overhead. Evaluation matching limits examined insertion positions, not setup cost or wall time.','','| Case | Conserving JD search (s) | Pauli JD search (s) | Evaluation-matched Pauli search (s) | Conserving / Pauli positions |','|---|---:|---:|---:|---:|']
 for case in cases:
  rr={r['ordering']:r for r in rows if r['case']==case};c=rr['class_refined_response'];q=rr['class_pauli_JD'];b=rr['class_pauli_JD_budget'];lines.append(f"| {case} | {c['search_seconds']:.3f} | {q['search_seconds']:.3f} | {b['search_seconds']:.3f} | {c['evaluations']:,} / {q['evaluations']:,} |")
 lines+=['','Full setup, both-start accepted moves, exact-positive rescores and objective counts are in optimizer_runtime.csv. Native Tranter V2 exact-search timings and reuse metadata are in timings.json. One V2 order supports two formula arms and its preparation is counted once. Exact positive-scaling reuse is not a fresh search.','','## Propagation interval differences','','- Fixed and published first-order row.seconds time evolution only, excluding metrics and output files. Compilation before evolution is excluded; first-call JIT inside evolution can be included.','- Published V2 Strang times merged-pulse assembly and evolution, excluding metrics and output files.','- Optimized row.seconds additionally includes pulse construction, compressed state output and the design hash; metrics and the state hash are excluded.','- These intervals must not be compared directly as algorithm timings. End-to-end case preparation has its own boundaries and may include verified structure reuse.','', 'Raw seconds permit workload diagnosis. Accuracy-versus-rotation and accuracy-versus-evaluation comparisons remain the controlled comparisons in this campaign. A wall-time speedup claim would need fresh, repeated, isolated runs with matching measurement boundaries.']
 (a.output/'TIMINGS.md').write_text('\n'.join(lines)+'\n');print(json.dumps(stats))
if __name__=='__main__':main()
