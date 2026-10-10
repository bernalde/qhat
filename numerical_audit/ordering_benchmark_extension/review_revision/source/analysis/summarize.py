#!/usr/bin/env python3
"""Paired, floor-aware comparisons generated solely from frozen raw rows."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,hashlib,json,math,statistics
from pathlib import Path
import numpy as np

FLOOR=1e-12
CASES=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]+[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
CAND='class_refined_response'
LABELS={'jw_signed':'Pauli signed','jw_magnitude':'Pauli magnitude','jw_lex':'Pauli lexicographic','jw_color':'Pauli coloring','parent_signed':'Parent signed','parent_magnitude':'Parent magnitude','parent_color':'Parent coloring','owned_signed':'Owned signed','owned_parent_magnitude':'Owned parent magnitude','owned_l1':'Owned block l1','owned_reverse':'Owned reverse signed','owned_color':'Owned coloring','class_lex':'Complete class lexicographic','class_magnitude':'Complete class magnitude','class_color':'Complete class coloring','class_suffix_l1':'Complete class suffix l1','class_bch_l1':'Complete class BCH l1','refined_magnitude':'Refined class magnitude','tranter_deplete_groups':'Tranter depleteGroups','tranter_equalise_groups':'Tranter equaliseGroups','hastings_interleaved_lex':'Hastings literal, lexical completion','bch1_full_l1_insertion':'BCH1 full-l1 insertion adaptation','hastings_interleaved_consolidated':'Hastings consolidated, lexical completion','tomesh_beam32':'Tomesh-inspired beam32','tranter_commutator':'Tranter commutator','tranter_reverse_commutator':'Tranter reverseCommutator',CAND:'Conserving response JD','class_refined_J0':'Conserving unweighted J0','class_pauli_JD':'Pauli response JD','class_pauli_JD_budget':'Pauli response JD, matched evaluations','strang_class_refined_response':'Strang, conserving JD order','strang_jw_magnitude':'Strang, Pauli magnitude','strang_jw_signed':'Strang, Pauli signed','strang_class_lex':'Strang, complete class lexicographic','tranter_v2_l1_first_order':'Tranter V2-l1 order, first-order use','tranter_v2_l1_strang':'Strang, Tranter V2-l1 insertion'}

def require(v,m):
    if not v:raise ValueError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def comparison(candidate,baseline):
    require(len(candidate)==len(baseline),'Unequal paired sample lengths')
    require(all(c['case']==b['case'] for c,b in zip(candidate,baseline)),'Misaligned paired cases')
    triples=[(c['case'],float(c['one_minus_overlap']),float(b['one_minus_overlap'])) for c,b in zip(candidate,baseline)]
    deltas=[x-y for _,x,y in triples]
    require(all(math.isfinite(x) and math.isfinite(y) for _,x,y in triples),'Nonfinite comparison')
    ratios=[(case,x/y) for case,x,y in triples if y>FLOOR]
    return {'wins':sum(d < -FLOOR for d in deltas),'ties':sum(abs(d)<=FLOOR for d in deltas),'losses':sum(d > FLOOR for d in deltas),'median_candidate_over_comparator':statistics.median(x for _,x in ratios) if ratios else None,'ratio_cases':len(ratios),'omitted_denominator_floor_cases':[case for case,_,y in triples if y<=FLOOR],'per_case':[{'case':case,'candidate_E':x,'comparator_E':y,'ratio':x/y if y>FLOOR else None} for case,x,y in triples]}

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--preview-without-v2',action='store_true',help='Interim layout preview only; cannot be packaged as complete');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    campaigns=[a.root/'fixed',a.root/'optimized',a.root/'published/extra',a.root/'published/v2'];expected=[2996,728,392,196];rows=[];inputs={}
    if a.preview_without_v2:campaigns=campaigns[:-1];expected=expected[:-1]
    for root,count in zip(campaigns,expected):
        if (root/'rows').exists():paths=sorted((root/'rows').rglob('*.json'))
        else:paths=sorted((root/'results').glob('*.b*.json'))
        require(len(paths)==count,f'Incomplete {root}: {len(paths)}/{count}')
        for f in paths:
            row=json.loads(f.read_text());row['campaign']=str(root.relative_to(a.root));row['row_path']=str(f.relative_to(a.root));row['timing_scope']='pulse_evolution_state_io_design_hash' if root.name=='optimized' else ('pulse_assembly_and_evolution' if row['formula_order']==2 else 'evolution_only');rows.append(row);inputs[row['row_path']]=sha(f)
    indexed={}
    for row in rows:
        token=(row['case'],row['ordering'],int(row['budget_units']));require(token not in indexed,'Duplicate numerical row '+str(token));indexed[token]=row
        require(row['status'] in ('complete','infeasible'),'Unexpected result status')
        require(row['pauli_rotations']<=row['rotation_budget'],'Over-budget row')
        if row['status']=='complete':require(row['steps']>0 and math.isfinite(float(row['one_minus_overlap'])),'Invalid completed row')
    names=sorted({r['ordering'] for r in rows if not r['ordering'].startswith('random_')});candidate=[indexed[(c,CAND,100)] for c in CASES];comparisons={}
    require(all(r['status']=='complete' for r in candidate),'Incomplete candidate')
    for name in names:
        if name==CAND:continue
        base=[indexed[(c,name,100)] for c in CASES];require(all(x['status']=='complete' for x in base),'Missing headline completed row')
        require(all(c['case']==b['case'] and c['T']==b['T'] and c['rotation_budget']==b['rotation_budget'] for c,b in zip(candidate,base)),'Unmatched comparison target/time/budget')
        comparisons[name]={'all':comparison(candidate,base),'development':comparison(candidate[:8],base[:8]),'prospective':comparison(candidate[8:],base[8:]),'formula_order':base[0]['formula_order']}
    strang={}
    for budget in (2,5,10,100):
        cc=[indexed[(c,'strang_class_refined_response',budget)] for c in CASES]
        strang[str(budget)]={}
        for name in names:
            if name=='strang_class_refined_response' or indexed[(CASES[0],name,budget)]['formula_order']!=2:continue
            bb=[indexed[(c,name,budget)] for c in CASES]
            require(all(r['status']=='complete' for r in cc+bb),'Incomplete symmetric comparison')
            strang[str(budget)][name]={'all':comparison(cc,bb),'development':comparison(cc[:8],bb[:8]),'prospective':comparison(cc[8:],bb[8:])}
    random={}
    for family in ('owned','refined','pauli'):
        medians=[];case_stats=[]
        for case in CASES:
            group=[indexed[(case,f'random_{family}_{i:02d}',100)] for i in range(20)];require(all(r['status']=='complete' for r in group),'Incomplete random distribution');values=[x['one_minus_overlap'] for x in group];q=np.quantile(values,[0,.25,.5,.75,1])
            medians.append({'case':case,'one_minus_overlap':float(q[2])});cand=indexed[(case,CAND,100)]['one_minus_overlap'];case_stats.append({'case':case,'min':float(q[0]),'q25':float(q[1]),'median':float(q[2]),'q75':float(q[3]),'max':float(q[4]),'sample_count':20,'candidate_better_than_samples':sum(cand-v < -FLOOR for v in values)})
        random[family]={'comparison_to_median':comparison(candidate,medians),'cases':case_stats}
    thresholds={}
    for eps in (1e-4,1e-6,1e-8):
        methods={}
        for name in names:
            items=[]
            for case in CASES:
                valid=[r for (c,n,b),r in indexed.items() if c==case and n==name and r['status']=='complete' and r['one_minus_overlap']<=eps]
                best=min(valid,key=lambda r:r['pauli_rotations']) if valid else None
                items.append({'case':case,'rotations':best['pauli_rotations'] if best else None,'budget_units':best['budget_units'] if best else None,'E':best['one_minus_overlap'] if best else None})
            methods[name]=items
        thresholds[str(eps)]=methods
    leakage={name:{'max':max(indexed[(c,name,100)]['particle_number_leakage'] for c in CASES),'per_case':[indexed[(c,name,100)]['particle_number_leakage'] for c in CASES]} for name in names}
    prep=[]
    for case in CASES:
        path=a.root/'optimized/designs'/f'{case}.design.json';d=json.loads(path.read_text())
        prep.append({'case':case,'U':len(d['coefficients']),'n_qubits':d['n_qubits'],'dimension':len(d['affine_basis_indices']),'reports':{n:{k:v for k,v in rep.items() if k in ('evaluations','rescored','search_seconds','pair_setup_seconds','units','nonzero_pair_count','evaluation_budget','J')} for n,rep in d['optimizer_reports'].items()},'JD_J0_same_order':d['optimizer_reports'][CAND]['order']==d['optimizer_reports']['class_refined_J0']['order']})
    result={'status':'interim_pending_v2' if a.preview_without_v2 else 'complete','rows':len(rows),'completed_evolutions':sum(x['status']=='complete' for x in rows),'infeasible':sum(x['status']=='infeasible' for x in rows),'methods_per_case':len(names)+60,'deterministic_methods':names,'cases':CASES,'analysis_floor':FLOOR,'comparison_budget':'100 times the number U of retained nonidentity Pauli words','comparisons':comparisons,'strang_comparisons':strang,'random':random,'least_tested_work':thresholds,'number_leakage_b100':leakage,'optimizer_preparation':prep,'source_sha256':sha(Path(__file__)),'input_sha256':inputs}
    (a.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Expanded ordering benchmark results','',('INTERIM: V2 comparison omitted while its exact search finishes. Not a complete campaign summary.' if a.preview_without_v2 else 'Complete frozen comparison.'),'',f"{len(rows)} planned rows; {result['completed_evolutions']} completed evolutions; {result['infeasible']} infeasible budgets. All 14 targets are fixed; four development molecular families at two geometries and three prospective families at two coefficient multipliers are not 14 independent molecular samples.",'','Candidate: complete nonzero flip classes with individually movable diagonal Paulis, determinant-response JD insertion. Comparisons below use budget 100U, T=1 and the same HF state and target. W/T/L refer to the candidate; absolute tie tolerance 1e-12. Ratios omit comparator denominators<=1e-12.','','| Comparator | Dev W/T/L | Prospective W/T/L | Overall W/T/L | Median E candidate/comparator | Ratio cases |','|---|---:|---:|---:|---:|---:|']
    for name,info in comparisons.items():
        def wtl(k):
            r=info[k];return f"{r['wins']}/{r['ties']}/{r['losses']}"
        rr=info['all'];ratio=rr['median_candidate_over_comparator'];lines.append(f"| {LABELS.get(name,name)} | {wtl('development')} | {wtl('prospective')} | {wtl('all')} | {ratio:.6g} | {rr['ratio_cases']} |" if ratio is not None else f"| {LABELS.get(name,name)} | {wtl('development')} | {wtl('prospective')} | {wtl('all')} | omitted (floor) | 0 |")
    lines+=['','## Symmetric use of the frozen conserving order','','At budget 10U (less affected by the error floor), W/T/L now favor the symmetric conserving-order formula. This is an empirical transfer of a first-order surrogate, not a second-order objective certificate.','','| Symmetric comparator | W/T/L | Median E conserving/comparator | Ratio cases |','|---|---:|---:|---:|']
    for name,info in strang['10'].items():
        r=info['all'];ratio=r['median_candidate_over_comparator'];ratio_text=f'{ratio:.6g}' if ratio is not None else 'omitted';lines.append(f"| {LABELS.get(name,name)} | {r['wins']}/{r['ties']}/{r['losses']} | {ratio_text} | {r['ratio_cases']} |")
    lines+=['','Fixed-random distributions are archived separately for every family/case; no order was selected using its measured error. Strang is a different formula order; its budgets account for exact commuting boundary fusion. Evaluation-matched Pauli search does not match classical wall time. Full methods, input hashes, threshold costs and leakage are in summary.json.']
    (a.output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    fields=['case','panel','ordering','formula_order','budget_units','rotation_budget','steps','pauli_rotations','status','one_minus_overlap','infidelity','particle_number_leakage','spin_sector_leakage','seconds','timing_scope','campaign','row_path']
    with (a.output/'all_rows.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps({k:result[k] for k in ('rows','completed_evolutions','infeasible','methods_per_case')}))
if __name__=='__main__':main()
