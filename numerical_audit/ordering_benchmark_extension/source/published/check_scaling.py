"""Check exact positive uniform-scale reuse and cubic/quadratic diagnostics."""
import json,math
from pathlib import Path
import algorithms as a
HERE=Path(__file__).resolve().parent

def main():
    output=[]
    for family in ('li2','f2','nh3'):
        for campaign in ('extra','v2'):
            p0=HERE/campaign/'designs'/f'{family}_coeff1p00.json';p1=HERE/campaign/'designs'/f'{family}_coeff1p25.json'
            if not p0.exists() or not p1.exists():raise RuntimeError('Incomplete scale pair '+family+' '+campaign)
            d0=json.loads(p0.read_text());d1=json.loads(p1.read_text());t0=a.block(d0['coefficients']);t1=a.block(d1['coefficients']);n0=int(d0['denominator']);n1=int(d1['denominator'])
            if t0.keys()!=t1.keys() or any(4*t1[w]*n0!=5*t0[w]*n1 for w in t0):raise RuntimeError('Targets are not exactly5/4scaled')
            num_factor=next(iter(t1.values()))//next(iter(t0.values()))
            if any(t1[w]!=num_factor*t0[w] for w in t0):raise RuntimeError('Not a common integer numerator multiplier')
            for name,m0 in d0['methods'].items():
                m1=d1['methods'][name];s0=m0['schedule'];s1=m1['schedule']
                if len(s0)!=len(s1) or any(w!=v or int(c1)!=num_factor*int(c0) for (w,c0),(v,c1) in zip(s0,s1)):raise RuntimeError('Scale-independent ordering changed '+name)
                if name=='tomesh_beam32':
                    if m0['metadata']['clique_order']!=m1['metadata']['clique_order'] or int(m1['metadata']['adjacent_noncommuting_pair_weight_numerator'])!=num_factor**2*int(m0['metadata']['adjacent_noncommuting_pair_weight_numerator']):raise RuntimeError('Tomesh scaled diagnostic mismatch')
                if campaign=='v2':
                    h0=m0['metadata']['history'];h1=m1['metadata']['history']
                    if any(i!=j or p!=q or int(b)!=num_factor**3*int(a) for (i,p,a),(j,q,b) in zip(h0,h1)):raise RuntimeError('V2 cubic history mismatch')
            output.append({'family':family,'campaign':campaign,'methods':len(d0['methods']),'physical_coefficient_multiplier':[5,4],'integer_numerator_multiplier':num_factor,'all_orders_identical':True,'diagnostics_scale_exactly':True})
    (HERE/'scaling_checks.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output))
if __name__=='__main__':main()
