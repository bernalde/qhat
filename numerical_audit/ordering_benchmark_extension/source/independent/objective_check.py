#!/usr/bin/env python3
"""Reconstruct optimized BCH scores using cumulative operators grouped by flip.

No production pair table or incremental move code is imported. All cumulative
operators and destination amplitudes use integers. Diagonal transport weights
are evaluated only after accumulating the actual commutator action.
"""
import argparse,json,math,resource,time
from pathlib import Path
from core import require,block,physical_masks,sha,object_sha

def descriptors(raw,n):return [(physical_masks(w,n)[1],physical_masks(w,n)[2]//2,c) for w,c in block(raw).items()]
def val(terms,v):return sum(c*(-1 if ((z&v).bit_count()+y)&1 else 1) for z,y,c in terms)
def bch_action(units,order,n,hf):
    prefix={};prefix0={};out={}
    desc=[descriptors(b,n) for b in units]
    flips=[physical_masks(next(iter(block(b))),n)[0] for b in units]
    for j in order:
        current=desc[j];x=flips[j];f0=val(current,hf)
        for y,terms in prefix.items():
            if x==y:continue
            d=x^y
            fj_at_y=val(current,hf^y)
            prev_at_x=sum(c*(-1 if ((z&(hf^x)).bit_count()+yy)&1 else 1) for (z,yy),c in terms.items())
            v=f0*prev_at_x-prefix0[y]*fj_at_y
            out[d]=out.get(d,0)+v
        pp=prefix.setdefault(x,{})
        for z,yy,c in current:
            k=(z,yy);pp[k]=pp.get(k,0)+c
            if not pp[k]:del pp[k]
        prefix0[x]=prefix0.get(x,0)+f0
    return {k:v for k,v in out.items() if v}

def objective(amplitudes,diagonal,hf,den,weighted):
    e0=val(diagonal,hf)/den;parts=[]
    for dest,v in amplitudes.items():
        if not dest:continue
        gap=(val(diagonal,hf^dest)/den-e0)/2
        rho=(math.sin(gap)/gap)**2 if weighted and gap else 1.
        parts.append(.25*rho*(v/(den*den))**2)
    return math.fsum(parts)

def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='*');a=p.parse_args();resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));a.output.mkdir(parents=True,exist_ok=True)
    files=sorted((a.campaign/'designs').glob('*.design.json'))
    if a.cases:files=[p for p in files if p.name.removesuffix('.design.json') in a.cases]
    checks=[]
    for path in files:
        d=json.loads(path.read_text());case=d['case'];out=a.output/(case+'.json')
        if out.exists():
            old=json.loads(out.read_text());require(old['design_sha256']==sha(path) and old['source_sha256']==sha(__file__),'Stale objective check');checks.append(old);continue
        start=time.monotonic();n=d['n_qubits'];hf=int(d['affine_basis_indices'][0]);den=int(d['denominator']);diag=[]
        for word,v in block(d['coefficients']).items():
            x,z,y=physical_masks(word,n)
            if not x:diag.append((z,y//2,v))
        report=[]
        for name,metadata in d['optimizer_reports'].items():
            units=d['blocks']['pauli' if 'pauli' in name else 'refined'];order=metadata['order'];require(sorted(order)==list(range(len(units))),'Invalid optimizer permutation')
            flattened=[[w,str(v)] for i in order for w,v in units[i]]
            require(object_sha(flattened)==d['methods'][name]['schedule_sha256'],'Order/schedule mismatch')
            amps=bch_action(units,order,n,hf);J=objective(amps,diag,hf,den,name!='class_refined_J0');delta=abs(J-metadata['J']);require(delta<max(1e-20,1e-10*max(J,metadata['J'])),'Objective mismatch '+case+'/'+name)
            evaluations=sum(x['evaluations'] for x in metadata['starts']);require(evaluations==metadata['evaluations'],'Evaluation count mismatch')
            for s in metadata['starts']:
                require(len(s['history'])<=10,'Too many accepted insertions')
                for h in s['history']:require(h['new_J']<h['old_J']-max(1e-30,1e-12*h['old_J']),'Accepted nonimprovement')
            if metadata['evaluation_budget'] is not None:require(evaluations<=metadata['evaluation_budget'],'Exceeded search budget')
            report.append({'ordering':name,'independent_J':J,'recorded_J':metadata['J'],'absolute_difference':delta,'amplitude_sha256':object_sha({str(k):str(v) for k,v in sorted(amps.items())}),'evaluations':evaluations})
        require(d['optimizer_reports']['class_pauli_JD_budget']['evaluation_budget']==d['optimizer_reports']['class_refined_response']['evaluations'],'Unmatched objective-evaluation budget')
        result={'case':case,'status':'passed','design_sha256':sha(path),'source_sha256':sha(__file__),'checks':report,'seconds':time.monotonic()-start};out.write_text(json.dumps(result,indent=2)+'\n');checks.append(result);print('OBJECTIVES',case,[r['absolute_difference'] for r in report],result['seconds'],flush=True)
    (a.output/'summary.json').write_text(json.dumps({'status':'passed','cases':len(checks),'objectives':sum(len(x['checks']) for x in checks),'checks':checks},indent=2)+'\n')
if __name__=='__main__':main()
