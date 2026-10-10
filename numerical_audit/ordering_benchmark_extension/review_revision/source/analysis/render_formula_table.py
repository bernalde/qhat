#!/usr/bin/env python3
"""Render the floor-resolved formula-order table from review recomputation."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--checks',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();j=json.loads(a.checks.read_text());labels={'strang_class_lex':'Complete-class label','strang_class_refined_response':'Conserving $J_D$','strang_jw_magnitude':'Pauli magnitude','strang_jw_signed':'Pauli signed','tranter_v2_l1_strang':r'Tranter $V_2$, coefficient-$\ell^1$'}
s=[r'\begin{table}[htbp]',r'\centering\small',r'\caption{First-order conserving candidate versus symmetric controls at $10U_c$. W/T/L favor the first-order candidate; medians divide its error by the symmetric comparator error. All fourteen comparator errors exceed the denominator floor.}',r'\label{tab:review-formula-order}',r'\begin{tabular}{lcc}',r'\toprule',r'Symmetric comparator & W/T/L & Median ratio ($n=14$)\\',r'\midrule']
for k,v in j['first_order_vs_strang']['10'].items():
 if v['ratio_cases']!=14:raise ValueError('Table caption needs changed floor coverage')
 s.append(f"{labels[k]} & {v['wins']}/{v['ties']}/{v['losses']} & {v['median_ratio']:.3g}"+r'\\')
s += [r'\bottomrule',r'\end{tabular}',r'\end{table}'];a.output.write_text('\n'.join(s)+'\n')
