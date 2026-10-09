#!/usr/bin/env python3
"""Post-run descriptive summary; it does not select schedules or states."""
import argparse,json
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);args=ap.parse_args()
    rows=[];checks=[];spin=[];coherent=[]
    for directory in sorted(args.root.glob('*q')):
        if not (directory/'results.json').exists():continue
        result=json.loads((directory/'results.json').read_text());design=json.loads((directory/'frozen_state_design.json').read_text())
        rows+=result['rows']+result['archived_HF_comparator_rows']
        checks.append({'case':design['case'],'seconds':result['elapsed_seconds'],'peak_rss_KiB':result['peak_rss_KiB'],
                       **json.loads((directory/'reference_checks.json').read_text())})
        spin.append({'case':design['case'],'states':design['state_info'],'spectrum':design['spectrum']})
        with np.load(directory/'states.npz',allow_pickle=False) as z:
            g=z['initial_sector_ground'];e=z['initial_sector_next_distinct_energy'];ref=z['exact_equal_coherent'];ar=np.vdot(g,ref);br=np.vdot(e,ref)
            for method in ('fermionic_signed','jw_signed','jw_magnitude'):
                v=z[method+'_equal_coherent'];v=v/np.linalg.norm(v);a=np.vdot(g,v);b=np.vdot(e,v)
                coherent.append({'case':design['case'],'ordering':method,
                     'relative_eigenpair_phase_error_radians':float(np.angle((b/a)/(br/ar))),
                     'probability_outside_initial_eigenpair':float(np.linalg.norm(v-a*g-b*e)**2),
                     'ground_minus_next_population':float(abs(a)**2-abs(b)**2),
                     'interpretation':'Post-run descriptive projection diagnostic, not used to select states or orderings. Relative phase is measured against the full-H exact reference.'})
    summary={'new_evolution_rows':len([r for r in rows if r['state']!='HF_archived']),
             'archived_HF_comparator_rows':len([r for r in rows if r['state']=='HF_archived']),
             'rows':rows,'checks':checks,'initial_state_diagnostics':spin,'coherent_diagnostics':coherent}
    for filename,obj in [('summary.json',summary),('coherent_diagnostics.json',coherent)]:
        with (args.root/filename).open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
    lines=['# H₂O larger-state control','',
           'This exploratory control addresses state dependence of two previously successful HF cases. It introduces no new ordering rule. Both are reconstructed HGBS-5 active spaces of one molecule, in the original frozen orbital gauge. They are not independent molecular samples.', '',
           'The target coefficients and three chronological schedules are exactly those of QHAT revision `0be7d4857537668eec957dabaa8438ea5c0d33f2`. Every new evolution uses T = 1 and r = 100 in the full Hilbert space, with the same omitted scalar. No particle- or spin-sector projection is applied during propagation.', '',
           'The state-preparation protocol was frozen before approximate errors: the first canonical vector of the complete ground-energy cluster, the first canonical vector of the next distinct-energy cluster in the full Nα = Nβ = 3 sector, and their equal coherent superposition. Initial arrays and hashes were written before approximate propagation. This differs deliberately from taking the second listed eigenvector when a ground cluster is degenerate.', '',
           '| Qubits | State | Owned E | JW signed E | JW magnitude E | JW signed / owned | JW magnitude / owned |',
           '|---:|---|---:|---:|---:|---:|---:|']
    for case in sorted({r['case'] for r in rows}):
        n=sum(map(int,case.split('_as-')[1].split('-')))
        for name in ('HF_archived','sector_ground','sector_next_distinct_energy','equal_coherent'):
            d={r['ordering']:r for r in rows if r['case']==case and r['state']==name}
            a=d['fermionic_signed']['one_minus_overlap'];b=d['jw_signed']['one_minus_overlap'];c=d['jw_magnitude']['one_minus_overlap']
            lines.append(f'| {n} | {name} | {a:.6e} | {b:.6e} | {c:.6e} | {b/a:.4g} | {c/a:.4g} |')
    lines+=['', 'E = 1 − |〈normalized exact, normalized approximate〉|, evaluated through a stable phase-aligned distance. A ratio above one favors owned ordering. HF rows are rechecked archived states, not new propagations.', '',
       'Both selected ground states are singlets and almost HF determinants; the next distinct-energy states are Ms = 0 triplets with negligible HF overlap. The coherent states have 〈S²〉 = 1 and Var(S²) = 1. They are spin mixtures, not excited singlets. Full-H eigenvector residuals and small off-sector target entries are measured and retained.', '',
       'The independent checks include written-letter sparse Hamiltonian assembly; tensor-axis Hamiltonian action and a sampled twelve-rotation product checked against QHAT; dense/ARPACK eigenvalue agreement; full sparse exponential references compared with split-time action and eigenphase predictions. The standalone saved-state checker imports no QHAT and recomputes hashes, stable errors, and number/spin leakage. It is not a second independent full product-formula propagation.', '',
       'Energy expectations are normalized and reported as drift from the exact evolved state, in Hartree. These are an additional observable, not interchangeable with overlap error. `coherent_diagnostics.json` reports relative phase and population outside the initial eigenpair; this is descriptive post-run analysis.', '',
       'The 14-qubit case was run first. Its 33-second runtime justified extending the already frozen protocol to 16 qubits. That extension was authorized independently of observed errors. Each process had a 6 GiB resident-memory cap and two threads.', '',
       'Run `recheck_saved.py --root <this directory> --output <new JSON path>` with the pinned environment; it supports Python `-O`. Scientific replay of `run_state_control.py` currently expects the original pinned QHAT checkout at `/home/bernalde/repos/qhat`; changing that checkout requires an isolated copy at the pinned revision or a separately documented path adaptation. Source snapshots and exact coefficient/schedule designs are included.', '',
       'These controls do not resolve orbital-sign sensitivity, geometry transfer, comparative gate cost, or performance on unseen molecules. In particular, success in the original orbital gauge is not evidence of invariance of a signed coefficient ordering under physically equivalent orbital phases.']
    with (args.root/'REPORT.md').open('x') as f:f.write('\n'.join(lines)+'\n')
    print('\n'.join(lines[:19]))

if __name__=='__main__':main()
