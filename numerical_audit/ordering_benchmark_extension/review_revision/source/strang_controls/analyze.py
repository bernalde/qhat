#!/usr/bin/env python3
"""Analyze all predeclared orientations without selecting on observed error."""
import os,resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,json,statistics,sys
from pathlib import Path
import run as r
ROOT=Path(__file__).resolve().parent
TOL=1e-12

def outcomes(pairs):
    counts={'win':0,'tie':0,'loss':0};ratios=[]
    for row in pairs:
        a,b=row['candidate_error'],row['comparator_error'];counts['win' if a<b-TOL else 'loss' if a>b+TOL else 'tie']+=1
        if b>TOL:ratios.append(a/b)
    return {'pairs':len(pairs),'outcomes':counts,'ratio_cases':len(ratios),'median_error_ratio':statistics.median(ratios) if ratios else None,'case_rows':pairs}

def masks(word):
    x=z=0
    for q,p in word:
        if p in 'XY':x|=1<<q
        if p in 'YZ':z|=1<<q
    return x,z

def endpoints(schedule):
    mm=[masks(word) for word,_ in schedule]
    def commutes(i,j):return ((mm[i][0]&mm[j][1])^(mm[i][1]&mm[j][0])).bit_count()%2==0
    a=b=0;U=len(mm)
    while a<U and all(commutes(a,j) for j in range(a)):a+=1
    while b<U and all(commutes(U-1-b,U-1-j) for j in range(b)):b+=1
    return a,b

def analyze():
    rows=json.loads((ROOT/'results.json').read_text());r.require(len(rows)==504 and len({(x['case'],x['ordering'],x['budget_units']) for x in rows})==504,'Campaign incomplete/duplicates');names=sorted({x['ordering'] for x in rows});by={(x['case'],x['ordering'],x['budget_units']):x for x in rows};plans=[]
    for case in r.CASES:
        d=r.verified_design(case)
        for name,method in d['methods'].items():
            U=len(method['schedule']);a,b=endpoints(method['schedule']);lengths=method['unit_lengths'];r.require(a>=lengths[0] and b>=lengths[-1],'Bad endpoint commutation')
            plans.append({'case':case,'ordering':name,'U':U,'prefix':a,'suffix':b,'overlap':max(0,a+b-U),'unit_first':lengths[0],'unit_last':lengths[-1]})
            for budget in r.BUDGETS:
                row=by[case,name,budget];slope=2*U-lengths[0]-lengths[-1];r.require(row['steps']==(budget*U-lengths[0])//slope,'Nonmaximal r');r.require(row['pauli_rotations']==row['steps']*slope+lengths[0],'Actual count mismatch');row['uniform_fused_rotations']=row['steps']*(2*U-a-b)+a;r.require(row['uniform_fused_rotations']<=row['pauli_rotations'],'Fusion increased cost')
    comparisons={};thresholds={}
    for orientation in ('forward','reverse'):
        candidate='class_refined_response__'+orientation;others=[n for n in names if n!=candidate];comparisons[candidate]={}
        for b in r.BUDGETS:
            report={}
            for name in others:
                pairs=[]
                for case in r.CASES:
                    a=by[case,candidate,b];bb=by[case,name,b];pairs.append({'case':case,'candidate_error':a['one_minus_overlap'],'comparator_error':bb['one_minus_overlap'],'candidate_rotations':a['pauli_rotations'],'comparator_rotations':bb['pauli_rotations'],'candidate_steps':a['steps'],'comparator_steps':bb['steps']})
                report[name]=outcomes(pairs)
            comparisons[candidate][str(b)]=report
        thresholds[candidate]={}
        for eps in (1e-4,1e-6,1e-8):
            report={}
            for name in others:
                counts={'lower_tested_level':0,'same_tested_level':0,'higher_tested_level':0,'candidate_only_success':0,'comparator_only_success':0,'neither_success':0};pairs=[]
                for case in r.CASES:
                    ar=[by[case,candidate,b] for b in r.BUDGETS if by[case,candidate,b]['one_minus_overlap']<=eps];br=[by[case,name,b] for b in r.BUDGETS if by[case,name,b]['one_minus_overlap']<=eps]
                    if not ar or not br:
                        counts['candidate_only_success' if ar else 'comparator_only_success' if br else 'neither_success']+=1;continue
                    a=ar[0];b=br[0];counts['lower_tested_level' if a['budget_units']<b['budget_units'] else 'higher_tested_level' if a['budget_units']>b['budget_units'] else 'same_tested_level']+=1;pairs.append({'case':case,'candidate_level':a['budget_units'],'comparator_level':b['budget_units'],'candidate_rotations':a['pauli_rotations'],'comparator_rotations':b['pauli_rotations'],'candidate_uniform_rotations':a['uniform_fused_rotations'],'comparator_uniform_rotations':b['uniform_fused_rotations'],'emitted_ratio_comparator_over_candidate':b['pauli_rotations']/a['pauli_rotations'],'uniform_ratio_comparator_over_candidate':b['uniform_fused_rotations']/a['uniform_fused_rotations']})
                report[name]={'outcomes':counts,'pairs':pairs,'warning':'Same tested budget level does not resolve minimum work, even when endpoint fusion yields slightly different exact emitted costs.'}
            thresholds[candidate][str(eps)]=report
    orientation={}
    for key in ('class_refined_response','class_pauli_JD','class_pauli_JD_budget'):
        orientation[key]={}
        for b in r.BUDGETS:orientation[key][str(b)]=outcomes([{'case':case,'candidate_error':by[case,key+'__forward',b]['one_minus_overlap'],'comparator_error':by[case,key+'__reverse',b]['one_minus_overlap']} for case in r.CASES])
    robust={}
    for b in r.BUDGETS:
        robust[str(b)]={}
        for key in ('class_pauli_JD','class_pauli_JD_budget'):
            counts={'all_orientation_pairs_conserving_better':0,'all_orientation_pairs_pauli_better':0,'orientation_dependent_or_tied':0};pairs=[]
            for case in r.CASES:
                cons=[by[case,'class_refined_response__'+o,b]['one_minus_overlap'] for o in ('forward','reverse')];pauli=[by[case,key+'__'+o,b]['one_minus_overlap'] for o in ('forward','reverse')];label='all_orientation_pairs_conserving_better' if max(cons)<min(pauli)-TOL else 'all_orientation_pairs_pauli_better' if max(pauli)<min(cons)-TOL else 'orientation_dependent_or_tied';counts[label]+=1;pairs.append({'case':case,'conserving_forward_reverse':cons,'pauli_forward_reverse':pauli,'outcome':label})
            robust[str(b)][key]={'outcomes':counts,'pairs':pairs,'scope':'Robustness diagnostic across both fixed orientations, not performance of an error-selected schedule.'}
    result={'rows':504,'new_evolutions':420,'reused_baselines':84,'cases':14,'budgets':r.BUDGETS,'tie_tolerance':TOL,'ratio_floor':TOL,'comparisons':comparisons,'thresholds':thresholds,'orientation_comparisons_forward_vs_reverse':orientation,'orientation_robustness':robust,'maximum_particle_number_leakage':{name:max(x['particle_number_leakage'] for x in rows if x['ordering']==name) for name in names},'maximum_norm_drift':max(abs(x['approx_norm']-1) for x in rows),'uniform_fusion_plans':plans,'source_sha256':r.sha(__file__),'inputs_sha256':{'results.json':r.sha(ROOT/'results.json'),'ALL_PLANS_FROZEN.json':r.sha(ROOT/'ALL_PLANS_FROZEN.json')}}
    r.write(ROOT/'summary.json',result)
    fields=sorted(set().union(*(x.keys() for x in rows)))
    with (ROOT/'all_rows.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    lines=['# Strang transfer and orientation controls','','All 14 previously frozen targets and JD orders are retained. Five new arms give 420 evolutions; 84 forward conserving states are reused with exact provenance. No search or data-dependent selection is performed. The first chronological unit is outermost. Pauli arms are singleton products; conserving reverse reverses units while retaining each commuting unit\'s word order. All pairings use the same H, HF state, T=1 and nominal Pauli-rotation budget.','','## Accuracy at 10U','','Absolute ties use 1e-12; median ratios include comparator errors above 1e-12. A smaller ratio favors conserving Strang.','','| Conserving orientation | Comparator | W/T/L | Median E ratio | Ratio cases |','|---|---|---:|---:|---:|']
    for candidate,bb in comparisons.items():
        for name,record in bb['10'].items():
            q=record['outcomes'];v=record['median_error_ratio'];lines.append(f"| {candidate} | {name} | {q['win']}/{q['tie']}/{q['loss']} | {v:.8g} | {record['ratio_cases']} |")
    lines+=['','## Orientation robustness at 10U','','| Pauli search | Conserving better for all four orientation pairs | Pauli better for all four | Orientation-dependent or tied |','|---|---:|---:|---:|']
    for key,v in robust['10'].items():q=v['outcomes'];lines.append(f"| {key} | {q['all_orientation_pairs_conserving_better']} | {q['all_orientation_pairs_pauli_better']} | {q['orientation_dependent_or_tied']} |")
    lines+=['','## First successful tested budget level','','Each entry is conserving lower / same / higher tested level. Same-level exact rotation differences are not evidence of a resolved minimum-work ranking. Success is checked separately; counts below use joint successes. No interpolation is used.','','| Conserving orientation | E target | Comparator | Lower/same/higher | Candidate-only / comparator-only / neither success |','|---|---|---|---:|---:|']
    for candidate,ee in thresholds.items():
        for eps,nn in ee.items():
            for name,item in nn.items():q=item['outcomes'];lines.append(f"| {candidate} | {eps} | {name} | {q['lower_tested_level']}/{q['same_tested_level']}/{q['higher_tested_level']} | {q['candidate_only_success']}/{q['comparator_only_success']}/{q['neither_success']} |")
    lines+=['','The uniform-fusion sensitivity retains each measured state and r, changing only exact word-rotation accounting via C_uniform(r)=r(2U-a-b)+a, where a,b are maximal pairwise-commuting prefix/suffix lengths. These are costs for the existing circuits, not simulations with additional steps. Detailed costs are in all_rows.csv and summary.json.','','## Limits','','Both orientations are reported individually. Reversal is equally good under first-order JD but its Strang error can differ; no error-informed orientation selector is presented. The matched-evaluation Pauli arm accepted only one or two moves per start in the original search, so it is a deliberately truncated best-improvement control, not evidence of general search-efficiency superiority. The larger-search arm is the substantive partition comparison. This controls search prescription and existing chosen orders, not separate reoptimization of a second-order objective. The panel comprises seven molecular families and correlated geometry/coefficient variants; no general arbitrary-state, operator-norm, or compiled hardware claim follows.','']
    (ROOT/'REPORT.md').write_text('\n'.join(lines));print(json.dumps({'rows':504,'at_10U':{k:{n:{kk:vv for kk,vv in x.items() if kk!='case_rows'} for n,x in v['10'].items()} for k,v in comparisons.items()}},indent=2))
if __name__=='__main__':analyze()
