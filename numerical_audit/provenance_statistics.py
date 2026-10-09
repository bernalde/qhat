#!/usr/bin/env python3
"""Read-only reconstruction of Albert's archived QHAT paper statistics.

Reads an explicitly selected checkout; writes a new JSON report outside it.
The bootstrap is an exploratory dependence sensitivity, not population inference.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,pathlib,subprocess
import numpy as np
from scipy import stats

def load_csv(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def pearson(x,y):
    r=stats.pearsonr(x,y)
    return {'r':float(r.statistic),'naive_iid_p_value':float(r.pvalue),'n_rows':len(x)}
def bootstrap(x,y,groups,seed=20261009,samples=10000,center=False):
    rng=np.random.default_rng(seed);names=sorted(set(groups));indices=[np.flatnonzero(np.array(groups)==g) for g in names];rs=[]
    for _ in range(samples):
        ix=np.concatenate([indices[i] for i in rng.integers(len(names),size=len(names))]);a=x[ix];b=y[ix]
        if center:a=a-a.mean();b=b-b.mean()
        den=np.linalg.norm(a)*np.linalg.norm(b)
        if den:rs.append(float(a@b/den))
    return {'clusters':len(names),'resamples':len(rs),'seed':seed,'percentile_95_interval':np.quantile(rs,[.025,.975]).tolist(),'interpretation':'exploratory dependence sensitivity on selected cases; not a population confidence claim'}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('study',type=pathlib.Path);parser.add_argument('output',type=pathlib.Path);args=parser.parse_args();base=args.study
    inputs=['analysis/jw_magnitude_descending_case_comparison.csv','analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv','analysis/cancellation_hypothesis_validation/full_ablation_results.csv','analysis/cancellation_hypothesis_validation/full_analysis/validation_statistics.csv']
    joined=load_csv(base/inputs[0]);published=load_csv(base/inputs[1]);raw=load_csv(base/inputs[2]);families=collections.defaultdict(dict)
    for row in raw:families[row['case_id']][row['schedule']]=row
    invariant=['number_of_pauli_terms','number_of_fermionic_terms','evolution_time','trotter_steps','coefficient_tolerance','bch_pair_abs_sum','nominal_exponential_count','exact_sector_dimension','number_of_unowned_pauli_terms']
    checks=[];logs=[];contrasts=[];molecules=[];centered=[]
    for case,schedules in sorted(families.items()):
        ref=schedules['fermionic_signed_reference'];signed=schedules['jw_signed_baseline'];mag=schedules['jw_magnitude_baseline'];g=list(schedules.values())
        checks.append({'case':case,'varying_invariant_fields':{k:[row.get(k) for row in g] for k in invariant if len({row.get(k) for row in g})!=1},'bch_pair_abs_sum':float(ref['bch_pair_abs_sum'])})
        a,b,c,d=[math.log10(float(row[field])) for row,field in [(mag,'bch_cancellation_ratio'),(ref,'bch_cancellation_ratio'),(mag,'one_minus_overlap'),(ref,'one_minus_overlap')]]
        logs.append([a,b,c,d]);molecules.append(ref['molecule'])
        x=math.log10(float(signed['bch_cancellation_ratio'])/float(mag['bch_cancellation_ratio']));y=math.log10(float(signed['one_minus_overlap'])/float(mag['one_minus_overlap']));contrasts.append([x,y])
        alt=np.array([[math.log10(float(r['bch_cancellation_ratio'])/float(ref['bch_cancellation_ratio'])),math.log10(float(r['one_minus_overlap'])/float(ref['one_minus_overlap']))] for r in [signed,mag]])
        centered.extend(alt-alt.mean(axis=0))
    logs=np.array(logs);contrasts=np.array(contrasts);centered=np.array(centered);x=logs[:,0]-logs[:,1];y=logs[:,2]-logs[:,3]
    identity=np.stack((contrasts/2,-contrasts/2),axis=1).reshape(-1,2)
    covariance=np.cov(logs.T,ddof=1);parts={'Cov(log Rmag, log Emag)':float(covariance[0,2]),'-Cov(log Rmag, log Eferm)':float(-covariance[0,3]),'-Cov(log Rferm, log Emag)':float(-covariance[1,2]),'Cov(log Rferm, log Eferm)':float(covariance[1,3])}
    cache={};joinchecks=[]
    for row in joined:
        entries=[]
        for tag,order in [('jw_signed','signed_coefficient_lexicographic'),('jw_magnitude','jw_magnitude_descending_lexicographic'),('fermionic_signed','fermionic_signed_coefficient_lexicographic')]:
            name=row[tag+'_source_csv']
            if not name:continue
            if name not in cache:cache[name]=load_csv(base/name)
            found=[r for r in cache[name] if r.get('tensor_path')==row['tensor_path'] and r.get('ordering')==order and r.get('status')=='success' and int(r['trotter_steps'])==int(row['trotter_steps']) and float(r['evolution_time'])==float(row['evolution_time']) and float(r.get('coefficient_tolerance') or '1e-12')==float(row['coefficient_tolerance'])]
            if len(found)!=1:raise RuntimeError(('ambiguous join',row['case_id'],tag,len(found)))
            entries.extend(found)
        unequal={k:sorted({e[k] for e in entries}) for k in ('number_of_pauli_terms','n_qubits','nominal_exponential_count','exact_sector_dimension') if len({e[k] for e in entries})!=1}
        if unequal:joinchecks.append({'case':row['case_id'],'unequal':unequal})
    out={'archive_commit':subprocess.check_output(['git','-C',str(base),'rev-parse','HEAD'],text=True).strip(),'analysis_script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'input_sha256':{p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in inputs},
      'historical_aggregate':{'rows':len(joined),'fixed_magnitude_outcomes':dict(collections.Counter(r['outcome'] for r in published)),'best_of_two_jw_outcomes':dict(collections.Counter(r['fermionic_vs_strongest_jw_winner'] for r in joined)),'source_csv_count':len(cache),'source_csv_sha256':{p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in cache},'selected_raw_schedule_rows':sum(2+bool(r['fermionic_signed_source_csv']) for r in joined),'invariant_field_mismatches':joinchecks,'absent_identifiers':['tensor content hash','Hamiltonian coefficient hash','source commit','initial state hash']},
      'full_panel':{'cases':len(families),'rows':len(raw),'molecule_families':dict(collections.Counter(molecules)),'invariant_checks':checks,'fallback_values':dict(collections.Counter(r['number_of_unowned_pauli_terms'] for r in raw)),'minimum_reported_error':min(float(r['one_minus_overlap']) for r in raw),'direct_cross_case_ratio_correlation':pearson(x,y),'direct_cross_case_family_bootstrap':bootstrap(x,y,molecules,center=True),'within_case_centered_correlation':pearson(*centered.T),'reference_cancellation_max_identity_residual':float(np.max(abs(centered-identity))),'within_case_unreplicated_jw_pair_contrast_correlation':pearson(*contrasts.T),'jw_pair_contrast_sign_agreements':int(np.sum(contrasts[:,0]*contrasts[:,1]>0)),'within_case_centered_case_bootstrap':bootstrap(contrasts[:,0],contrasts[:,1],list(range(len(contrasts)))),'within_case_centered_family_bootstrap':bootstrap(contrasts[:,0],contrasts[:,1],molecules),'direct_ratio_covariance':float(np.cov(x,y,ddof=1)[0,1]),'direct_ratio_covariance_decomposition':parts,'interpretation':'Repeated molecular families, selected favorable/control panel, and correlated per-case contrasts prevent treating all plotted rows as independent prospective tests. Common reference does not alone prove spurious correlation.'}}
    later=[]
    for campaign in ('codex_mechanism_experiment/runs/mechanism-t06qcgac','codex_external_validation/runs/transfer-vtrgy3ds'):
        path=base/campaign/'results.json';rows=json.loads(path.read_text());coherence=[r for r in rows if isinstance(r.get('coherence_factor'),(int,float))]
        later.append({'campaign':campaign,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(rows),'rows_with_defect_decomposition':len(coherence),'maximum_C_over_steps':max(r['coherence_factor']/r['steps'] for r in coherence),'rows_with_C_below_one':sum(r['coherence_factor']<1 for r in coherence),'max_I_minus_D_times_C':max(abs(r['infidelity']-r['projected_local_power']*r['coherence_factor']) for r in coherence),'interpretation':'Internal consistency of recorded outputs, not a new evolution rerun.'})
    out['later_campaign_record_checks']=later
    with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k in ['archive_commit','historical_aggregate']},indent=2));print(json.dumps({k:v for k,v in out['full_panel'].items() if k not in ['invariant_checks']},indent=2))
if __name__=='__main__':main()
