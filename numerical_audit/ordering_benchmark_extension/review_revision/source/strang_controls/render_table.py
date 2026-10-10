#!/usr/bin/env python3
import resource
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
import json
from pathlib import Path
root=Path(__file__).resolve().parent;s=json.loads((root/'summary.json').read_text())
labels={'class_pauli_JD__forward':'Larger Pauli search, forward','class_pauli_JD__reverse':'Larger Pauli search, reverse','class_pauli_JD_budget__forward':'Truncated Pauli search, forward','class_pauli_JD_budget__reverse':'Truncated Pauli search, reverse'}
lines=[r'\begin{table}[htbp]',r'\centering',r'\small',r'\begin{tabular}{lrrrr}',r'\toprule',r' & \multicolumn{2}{c}{Conserving forward} & \multicolumn{2}{c}{Conserving reverse} \\',r'Pauli Strang comparator & W/T/L & Median ratio & W/T/L & Median ratio \\',r'\midrule']
for name,label in labels.items():
 cols=[]
 for orient in ('forward','reverse'):
  x=s['comparisons']['class_refined_response__'+orient]['10'][name];q=x['outcomes'];cols+=[f"{q['win']}/{q['tie']}/{q['loss']}",f"{x['median_error_ratio']:.3g}"]
 lines.append(label+' & '+' & '.join(cols)+r' \\')
x=s['comparisons']['class_refined_response__forward']['10']['class_refined_response__reverse'];q=x['outcomes']
lines.append('Conserving reverse & '+f"{q['win']}/{q['tie']}/{q['loss']} & {x['median_error_ratio']:.3g} & --- & ---"+r' \\')
lines += [r'\bottomrule',r'\end{tabular}',r'\caption{Strang use of the frozen Hartree--Fock-response orders at rotation budget $10U_c$, across all 14 cases. W/T/L compares the conserving error with the comparator using absolute tolerance $10^{-12}$; a median error ratio below one favors the conserving order. Every reported ratio includes all 14 cases with comparator error above the floor. Forward uses the first unit of the frozen order outermost; reverse reverses unit order. No order is reoptimized for Strang or selected using these errors. The truncated Pauli search is the previous insertion-position-budget-matched best-improvement search, which accepted only one or two moves per start.}',r'\label{tab:strang-pauli-response-controls}',r'\end{table}','']
# Do not silently claim 14 ratio cases if the floor actually removed any.
if any(s['comparisons']['class_refined_response__'+o]['10'][n]['ratio_cases']!=14 for o in ('forward','reverse') for n in labels):
 lines=[line.replace('Every reported ratio includes all 14 cases with comparator error above the floor.','Ratios include only comparator errors above $10^{-12}$; the complete counts appear in the accompanying raw-data report.') for line in lines]
(root/'controls_table.tex').write_text('\n'.join(lines))
