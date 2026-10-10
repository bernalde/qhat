#!/usr/bin/env python3
"""Independent final fixed-row arithmetic and read-only manuscript correspondence."""
import hashlib,json,math,re,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT.parent/'paper'
FIX=ROOT/'fixed'
FLOOR=1e-12
CASES=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]+[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
def require(v,m):
    if not v:raise ValueError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compare(aa,bb):
    require(len(aa)==len(bb)==14,'Need 14 paired cases')
    require(all(a['case']==b['case'] for a,b in zip(aa,bb)),'Case mismatch')
    delta=[a['one_minus_overlap']-b['one_minus_overlap'] for a,b in zip(aa,bb)]
    ratios=[a['one_minus_overlap']/b['one_minus_overlap'] for a,b in zip(aa,bb) if b['one_minus_overlap']>FLOOR]
    return dict(wins=sum(v < -FLOOR for v in delta),ties=sum(abs(v)<=FLOOR for v in delta),losses=sum(v>FLOOR for v in delta),median_ratio=statistics.median(ratios) if ratios else None,ratio_count=len(ratios))
paths=sorted((FIX/'results').glob('*.b*.json'))
require(len(paths)==2996,'Wrong fixed row count')
rows=[json.loads(p.read_text()) for p in paths]
idx={(r['case'],r['ordering'],r['budget_units']):r for r in rows}
require(len(idx)==len(rows),'Duplicate rows')
case_counts={}
for case in CASES:
    rr=[r for r in rows if r['case']==case]
    require(len(rr)==214,'Wrong per-case count '+case)
    require(sum(r['status']=='complete' for r in rr)==206,'Wrong complete count '+case)
    for r in rr:
        require(r['T']==1 and r['formula_order']==1,'Fixed formula/time mismatch')
        require(r['pauli_rotations']<=r['rotation_budget'],'Budget exceeded')
        require(r['status'] in ('complete','infeasible'),'Invalid row status')
        if r['status']=='complete':
            require(r['steps']>0 and 0<=r['one_minus_overlap']<=1 and math.isfinite(r['one_minus_overlap']),'Invalid error')
            state=FIX/'results'/r['state_file']
            require(state.is_file(),'Missing state')
        else:require(r['steps']==0,'Infeasible row has nonzero steps')
    case_counts[case]={'planned':214,'complete':206,'infeasible':8}
candidate=[json.loads((ROOT/'optimized/rows'/c/'class_refined_response.b100.json').read_text()) for c in CASES]
require(all(r['status']=='complete' for r in candidate),'Incomplete candidate')
methods=sorted({r['ordering'] for r in rows if not r['ordering'].startswith('random_')})
require(len(methods)==22,'Wrong deterministic method count')
comp={}
for method in methods:
    bb=[idx[(c,method,100)] for c in CASES]
    require(all(a['T']==b['T'] and a['rotation_budget']==b['rotation_budget'] for a,b in zip(candidate,bb)),'Unmatched time/budget')
    comp[method]=compare(candidate,bb)
random={}
for family in ('owned','refined','pauli'):
    median=[];wins=[]
    for c,a in zip(CASES,candidate):
        vals=[idx[(c,f'random_{family}_{i:02d}',100)]['one_minus_overlap'] for i in range(20)]
        median.append({'case':c,'one_minus_overlap':statistics.median(vals)})
        wins.append({'case':c,'candidate_better_than_samples':sum(a['one_minus_overlap']-v < -FLOOR for v in vals),'sample_count':20,'median_random_E':statistics.median(vals),'min_random_E':min(vals),'max_random_E':max(vals)})
    random[family]={'comparison_to_median':compare(candidate,median),'per_case':wins}
result={'status':'pass','row_count':len(rows),'complete':2884,'infeasible':112,'cases':case_counts,'fixed_comparisons':comp,'random_comparisons':random,'source_sha256':sha(Path(__file__)),'fixed_results_sha256':sha(FIX/'results.json'),'error_floor':FLOOR,'state_check':'Every completed row has its saved state file; hash/metric propagation replay is a separate root audit.'}
(FIX/'final_numerical_review.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Independent fixed-campaign final review','','All 2,996 planned rows are present: 2,884 completed evolutions and 112 explicitly infeasible budgets. Each of 14 targets has 206 completed and eight infeasible rows. Checked unique method/case/budget keys, valid status, time and formula order, nonnegative finite errors, cost within the common budget, and saved-state presence. State hashes and independent propagation are covered by the root audit; this script does not claim to redo them.','','At budget $100U$, the fresh conserving response candidate is paired with the same target, HF initial state, evolution time and rotation budget. The comparison below uses the absolute $10^{-12}$ tie floor. These 14 targets are seven molecular families at paired geometries or coefficient scales, not 14 independent molecules.','','| Fixed first-order comparator | Candidate W/T/L | Median error ratio |','|---|---:|---:|']
for m,v in comp.items():lines.append(f"| {m} | {v['wins']}/{v['ties']}/{v['losses']} | {v['median_ratio']:.8g} |")
lines+=['','The candidate wins against all 22 deterministic fixed first-order controls in all 14 targets. This is a statement about this specified panel, not a general ordering theorem. Classical preprocessing cost and second-order formulas require the separate optimized/published comparisons.','','| Fixed-random family | Candidate versus sample median W/T/L | Median ratio | Strictly better than individual samples |','|---|---:|---:|---:|']
for fam,v in random.items():
    d=v['comparison_to_median'];nw=sum(t['candidate_better_than_samples'] for t in v['per_case'])
    lines.append(f"| {fam} | {d['wins']}/{d['ties']}/{d['losses']} | {d['median_ratio']:.8g} | {nw}/280 |")
lines+=['','Random permutations use 20 predeclared seeds per family and are reused across every time step. No measured-error selection is made. Stronger controls matter: the full-budget unconstrained Pauli response optimizer wins 10 of 14 against the conserving first-order candidate, and symmetric formulas often do better at the same rotation budget.','','## Historical baseline correction','','Of 96 historical schedule checks, 64 reproduce literal historical flat-sequence hashes and 32 reproduce the same ordered commuting blocks with different internal word order. Those 32 ideal products agree exactly for every time and step count. Saved historical versus fresh states differ by at most 1.1998993524231055e-14 in vector norm and 6.0507798587110945e-18 in the reported error. See HISTORICAL_BASELINE_ALIGNMENT.md and historical_internal_order_check.json. No fixed measured schedule was changed after this distinction was identified.','','All frozen sources, 14 design hashes and eight QHAT dependency hashes remain unchanged at campaign completion (final_execution_check.json). Row timing measures evolution only; use timing_scope.json before comparing times across campaigns.']
require(all(v['wins']==14 for v in comp.values()),'The prose claim all22 wins is not supported')
(FIX/'FINAL_REVIEW.md').write_text('\n'.join(lines)+'\n')
lean=PAPER/'theory/ordering_mechanism/lean/FermionicTrotter/DiagonalResponseBridge.lean'
guide=PAPER/'theory/ordering_mechanism/lean_proof_guide.md'
tex=PAPER/'theory/diagonal_response_bridge_proofs.tex'
source=lean.read_text();appendix=tex.read_text();names=re.findall(r'^theorem (\w+)',source,re.M)
paragraphs=re.findall(r'\\paragraph\{(\d+)\. \\texttt\{([^}]+)\}\.\}',appendix)
covered=[n.replace('\\_','_') for _,n in paragraphs]
require(len(names)==19 and names==covered,'Appendix does not cover 19 source theorems in order')
proof_audit={'status':'pass','review_kind':'read-only human mathematical/source correspondence plus exact theorem-name coverage; no Lean recompilation','covered_theorems':names,'inputs':{str(p.relative_to(PAPER)):sha(p) for p in (lean,guide,tex,PAPER/'theory/ordering_mechanism/lean/FermionicTrotter/Algebra.lean')},'scope_checks':['Arbitrary real matrices; no imported Hermitian, number conservation, nonzero coefficient or body-rank hypothesis.','Zero orbitals allowed; arbitrary label type, finite list and repeated labels/masks retained.','Commutator orientation AB-BA and head-tail BCH recursion agree with Algebra.lean.','The diagonal pair-amplitude identity has no support hypothesis; diagonal commutator entry does.','Gauge acts on actual matrices and extracted amplitudes, not an assumed free scalar model.','The appendix explicitly excludes sinc integration, the full variance identity, complete JD, ranking and insertion-search formalization.'],'issues':[]}
(FIX/'lean_bridge_correspondence_audit.json').write_text(json.dumps(proof_audit,indent=2)+'\n')
(FIX/'LEAN_BRIDGE_CORRESPONDENCE.md').write_text('''# Read-only Lean bridge appendix correspondence\n\nPassed: all 19 public theorem names appear exactly once, in source order, in the new manuscript appendix. The mathematical steps follow the frozen proof guide and actual Lean source. Source hashes and the complete statement list are recorded in lean_bridge_correspondence_audit.json. No Lean recompilation was performed or claimed in this transcription check.\n\nThe support-to-weighted-shift derivation, literal matrix multiplication, pair commutator, molecular-mask specialization, actual orbital-sign conjugation, diagonal energy-gap identity and list BCH accumulation all preserve their Lean hypotheses and orientation. In particular, the scalar diagonal identity requires no support condition, whereas its commutator-entry specialization does; the appendix correctly distinguishes these statements. Repeated labels and repeated masks are allowed, and the label type need not be finite.\n\nThe proof of the BCH-entry identity expands exactly the imported head-tail recursion `comm head tail.sum + bch tail`. The prose invokes the earlier entry specialization for readability where Lean simplifies the same shift conditional directly; this is an equivalent mathematical step, not an alternative hypothesis or result.\n\nThe appendix correctly labels these as supporting formalizations rather than 19 new contributions and explicitly leaves the sinc integral, variance identity, complete response objective, ranking condition and numerical search outside the Lean bridge. No mathematical or source-correspondence defect was found. The external Ryabinkin locator and the prior 42-constant audit claim were not independently re-reviewed in this source-transcription check.\n''')
print(json.dumps({'fixed_status':'pass','fixed_rows':len(rows),'comparators':len(comp),'random':{k:v['comparison_to_median'] for k,v in random.items()},'lean_theorems_covered':len(names)},indent=2))
