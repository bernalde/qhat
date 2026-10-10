"""Paired headline comparisons against the previously frozen response candidate."""
from pathlib import Path
import json,statistics,sys,zipfile
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent));from common import *

def main():
    prior=json.loads((MECH/'diagonal_refinement/results/results.json').read_text())
    with zipfile.ZipFile(HELDOUT) as z:prior+=json.loads(z.read('results/results.json'))
    candidate={r['case']:r for r in prior if r['ordering']=='class_refined_response' and r['steps']==100 and r['case'] in CASES}
    require(len(candidate)==14,'Missing archived candidate');summary={};values=[]
    for campaign in ('extra','v2'):
        result=HERE/campaign/'results.json'
        if not result.exists():continue
        rows=json.loads(result.read_text());methods=sorted({r['ordering'] for r in rows})
        for name in methods:
            stats={}
            for panel,cases in [('development',DEV_CASES),('prospective_followup',PROS_CASES),('all',CASES)]:
                sample=[r for r in rows if r['case'] in cases and r['ordering']==name and r['budget_units']==100 and r['status']=='complete'];wins=ties=losses=0;ratios=[]
                for row in sample:
                    c=candidate[row['case']];require(c['pauli_rotations']==row['rotation_budget'],'Candidate budget mismatch');e=c['one_minus_overlap'];base=row['one_minus_overlap'];delta=e-base
                    if abs(delta)<=1e-12:ties+=1
                    elif delta<0:wins+=1
                    else:losses+=1
                    if base>1e-12:ratios.append(e/base)
                    if panel=='all':values.append({'case':row['case'],'comparator':name,'candidate_E':e,'comparator_E':base,'candidate_rotations':c['pauli_rotations'],'comparator_rotations':row['pauli_rotations'],'available_budget':row['rotation_budget'],'comparator_steps':row['steps'],'comparator_number_leakage':row['particle_number_leakage']})
                stats[panel]={'candidate_win_tie_loss':[wins,ties,losses],'median_candidate_over_comparator_E':statistics.median(ratios) if ratios else None,'ratios_with_denominator_above_1e_minus12':len(ratios),'max_comparator_number_leakage':max((r['particle_number_leakage'] for r in sample),default=0.)}
            summary[name]=stats
    write_atomic(HERE/'published_comparison_summary.json',{'candidate_provenance':'previously frozen exact same target/state/candidate at r100; cross-check new optimized campaign before manuscript use','tie_tolerance':1e-12,'statistics':summary,'paired_rows':values});print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
