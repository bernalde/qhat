#!/usr/bin/env python3
"""Generate manuscript tables and figures from the complete paired summary."""
import os
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k]="1"
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(__file__).resolve().parent/'mpl-cache'))
import argparse,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from summarize import LABELS,CASES,CAND

def tex(s):return s.replace('_',r'\_').replace('&',r'\&')
def wtl(d):return f"{d['wins']}/{d['ties']}/{d['losses']}"
def main():
    p=argparse.ArgumentParser();p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();s=json.loads((a.analysis/'summary.json').read_text());a.output.mkdir(parents=True,exist_ok=True)
    comp=s['comparisons'];names=list(comp)
    selected=['jw_signed','jw_magnitude','owned_signed','refined_magnitude','hastings_interleaved_consolidated','tranter_deplete_groups','tranter_equalise_groups','bch1_full_l1_insertion','class_refined_J0','class_pauli_JD','class_pauli_JD_budget','strang_jw_magnitude','strang_class_refined_response']
    if 'tranter_v2_l1_first_order' in names:selected.insert(8,'tranter_v2_l1_first_order')
    if 'tranter_v2_l1_strang' in names:selected.append('tranter_v2_l1_strang')
    selected=[n for n in selected if comp[n]['formula_order']==1]
    def make_table(keys,caption,label):
        if s['status']!='complete':caption='Interim layout preview; the V2 comparator is still pending. '+caption
        lines=[r'\begin{table}[htbp]',r'\centering',r'\small',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\begin{tabular}{p{.43\textwidth}ccc}',r'\toprule',r'Comparator & Dev. W/T/L & Follow-up W/T/L & Ratio ($n$)\\',r'\midrule']
        for n in keys:
            c=comp[n];m=c['all']['median_candidate_over_comparator'];ratio='floor (0)' if m is None else f"{m:.3g} ({c['all']['ratio_cases']})"
            if c['formula_order']==2 and c['all']['ratio_cases']<14:ratio += r' $^{\dagger}$'
            lines.append(tex(LABELS.get(n,n))+' & '+wtl(c['development'])+' & '+wtl(c['prospective'])+' & '+ratio+r'\\')
        lines.extend([r'\bottomrule',r'\end{tabular}',r'\end{table}']);return '\n'.join(lines)+'\n'
    (a.output/'headline_table.tex').write_text(make_table(selected,'Conserving response optimizer versus selected controls at budget $100U_c$. W/T/L counts favor the conserving first-order optimizer. Reported ratios are medians of $E_{\\rm conserving}/E_{\\rm comparator}$; smaller favors it. Denominators at or below $10^{-12}$ are omitted; $n$ is the number of included ratios. All displayed controls here are first order. Follow-up means Li$_2$/F$_2$/NH$_3$ on previously revealed targets, not a new holdout.','tab:expanded-headline'))
    chunks=[names[i:i+18] for i in range(0,len(names),18)]
    (a.output/'all_comparisons.tex').write_text('\n'.join(make_table(chunk,'All deterministic paired comparisons at budget $100U_c$; conventions match Table~\\ref{tab:expanded-headline}. Strang ratios marked $^{\\dagger}$ are floor-limited, using only the indicated $n$ cases; their omitted cases are those with the most accurate comparators.','tab:expanded-all-'+str(i)) for i,chunk in enumerate(chunks)))
    # State errors and least sampled work. Both numerical definitions are in
    # the figure caption/main text; points below floor are visibly capped.
    import csv
    rows=list(csv.DictReader((a.analysis/'all_rows.csv').open()));rmap={(r['case'],r['ordering'],int(r['budget_units'])):r for r in rows}
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,1,figsize=(6.5,6.8),gridspec_kw={'height_ratios':[1.15,1]},layout='constrained')
    styles=[(CAND,'Conserving $J_D$','#0072B2','o'),('class_pauli_JD','Pauli $J_D$','#CC79A7','s'),('class_pauli_JD_budget','Pauli $J_D$, capped evaluations','#E69F00','^'),('strang_jw_magnitude','Strang, Pauli magnitude','#009E73','D'),('strang_class_refined_response','Strang, conserving order','#555555','v')]
    x=np.arange(14)
    for at,(name,label,color,marker) in enumerate(styles):
        y=np.array([float(rmap[(case,name,100)]['one_minus_overlap']) for case in CASES]);axes[0].scatter(x+(at-2)*.11,np.maximum(y,1e-12),s=24,c=color,marker=marker,label=label)
    labels=[c.replace('_s','\ng').replace('_coeff','\nc').replace('1p','1.').replace('f2','F2').replace('li2','Li2').replace('nh3','NH3') for c in CASES]
    axes[0].set_xticks(x,labels,fontsize=7);axes[0].set_yscale('log');axes[0].axvline(7.5,color='.8',lw=.8);axes[0].axhline(1e-12,color='.6',lw=.8,ls=':');axes[0].set_ylabel('$E$ at rotation budget $100U_c$');axes[0].set_title('(a) Equal target, state, time and rotation budget',loc='left',pad=48,fontsize=9);axes[0].legend(ncol=2,fontsize=7,loc='lower left',bbox_to_anchor=(0,1.02),borderaxespad=0);axes[0].grid(axis='y',alpha=.2)
    choices=['jw_signed','jw_magnitude','hastings_interleaved_consolidated','class_pauli_JD','class_pauli_JD_budget','strang_jw_magnitude','strang_class_refined_response'];threshold=s['least_tested_work']['1e-06'];candidate={r['case']:r['rotations'] for r in threshold[CAND]}
    ticks=[]
    for j,name in enumerate(choices):
        vals=[(i,r['rotations']/candidate[r['case']]) for i,r in enumerate(threshold[name]) if r['rotations'] is not None and candidate[r['case']] is not None]
        axes[1].scatter([j+(i-6.5)*.025 for i,v in vals],[v for i,v in vals],c=['#0072B2' if i<8 else '#D55E00' for i,v in vals],s=22,alpha=.85)
        short={'jw_signed':'Pauli\nsigned','jw_magnitude':'Pauli\nmagnitude','hastings_interleaved_consolidated':'Hastings\nconsolidated','class_pauli_JD':'Pauli $J_D$','class_pauli_JD_budget':'Pauli $J_D$\ncapped eval.','strang_jw_magnitude':'Strang\nPauli mag.','strang_class_refined_response':'Strang\nconserving'}
        ticks.append(short[name]+f'\n({len(vals)}/14)')
    axes[1].set_xticks(range(len(choices)),ticks,fontsize=7);axes[1].set_yscale('log',base=2);axes[1].axhline(1,color='.4',ls='--',lw=1);axes[1].set_ylabel('$C_{\\mu,c}(10^{-6})/C_{\\mathrm{cons},c}(10^{-6})$');axes[1].set_title('(b) Least tested rotation count reaching $E\\leq10^{-6}$',loc='left',fontsize=9);axes[1].grid(axis='y',alpha=.2)
    fig.savefig(a.output/'expanded_ordering_comparison.pdf');fig.savefig(a.output/'expanded_ordering_comparison.png',dpi=180);plt.close(fig)
    # Standalone full work curves, suitable for the research supplement.
    fig,axes=plt.subplots(7,2,figsize=(10,15),layout='constrained',sharex=True)
    for case,ax in zip(CASES,axes.flat):
        U=next(r['U'] for r in s['optimizer_preparation'] if r['case']==case)
        for name,label,color,marker in styles:
            rr=[r for r in rows if r['case']==case and r['ordering']==name and r['status']=='complete'];rr.sort(key=lambda r:int(r['pauli_rotations']))
            ax.plot([int(r['pauli_rotations'])/U for r in rr],[max(float(r['one_minus_overlap']),1e-12) for r in rr],marker=marker,color=color,label=label,ms=3,lw=1)
        ax.set(xscale='log',yscale='log',title=case.replace('_s',' geometry ').replace('_coeff',' coefficient ').replace('1p','1.').replace('f2','F2').replace('li2','Li2').replace('nh3','NH3'),ylabel='$E$');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=6);axes[-1,0].set_xlabel('Emitted rotations / $U_c$');axes[-1,1].set_xlabel('Emitted rotations / $U_c$');fig.savefig(a.output/'all_work_curves.pdf');plt.close(fig)
    print('Wrote tables and publication figures to',a.output)
if __name__=='__main__':main()
