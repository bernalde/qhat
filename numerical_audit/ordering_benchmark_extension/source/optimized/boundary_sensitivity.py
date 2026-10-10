"""Uniform commuting-endpoint compilation; postprocess existing states only."""
from pathlib import Path
import json,math,statistics,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c

def anti(a,b):
    x,z,_=c.masks(a);xx,zz,_=c.masks(b)
    return bool(((x&zz)^(z&xx)).bit_count()%2)
def endpoints(words):
    U=len(words);a=0;b=0
    while a<U and all(not anti(words[a],words[j]) for j in range(a)):a+=1
    while b<U and all(not anti(words[U-1-b],words[U-1-j]) for j in range(b)):b+=1
    return a,b,max(0,a+b-U)
def fused_count(U,a,b,r):return r*(2*U-a-b)+a

def fused_pulses(U,a,b,r):
    c.require(r>=1,'Nonpositive step count');overlap=max(0,a+b-U);out=[]
    if overlap:
        global_words=list(range(U-b,a));prefix=list(range(U-b));middle=[];suffix=list(range(a,U))
    else:global_words=[];prefix=list(range(a));middle=list(range(a,U-b));suffix=list(range(U-b,U))
    out.extend((i,2*r) for i in global_words)
    if len(global_words)==U:
        c.require(len(out)==fused_count(U,a,b,r) and all(mult==2*r for _,mult in out),'All-commuting pulse/count failure')
        return out
    out.extend((i,1) for i in prefix)
    for step in range(r):
        out.extend((i,1) for i in middle);out.extend((i,2) for i in suffix);out.extend((i,1) for i in reversed(middle));out.extend((i,1 if step==r-1 else 2) for i in prefix)
    totals=[0]*U
    for i,mult in out:totals[i]+=mult
    c.require(totals==[2*r]*U and len(out)==fused_count(U,a,b,r),'Boundary fusion pulse/count failure')
    return out

def dense_checks():
    pauli={'I':np.eye(2,dtype=complex),'X':np.array([[0,1],[1,0]],dtype=complex),'Y':np.array([[0,-1j],[1j,0]],dtype=complex),'Z':np.diag([1.,-1.]).astype(complex)}
    def matrix(w,n):
        tags=dict(w);out=np.array([[1.+0j]])
        for i in range(n):out=np.kron(out,pauli[tags.get(i,'I')])
        return out
    def product(pulses,matrices,coeffs,r):
        out=np.eye(len(matrices[0]),dtype=complex);I=np.eye(len(out),dtype=complex)
        for i,mult in pulses:
            angle=coeffs[i]*mult/(2*r);out=(math.cos(angle)*I-1j*math.sin(angle)*matrices[i])@out
        return out
    rng=np.random.default_rng(2026100917);toys=[('overlap_AB_BC_but_not_AC',2,[((0,'X'),),((1,'Z'),),((0,'Z'),)]),('all_commuting',2,[((0,'X'),(1,'X')),((0,'Y'),(1,'Y')),((0,'Z'),(1,'Z'))]),('nonoverlap',2,[((0,'X'),),((0,'Z'),),((1,'X'),),((1,'Z'),)])]
    from itertools import product as cartesian
    for it in range(80):
        n=2+it%2;pool=[tuple((i,p) for i,p in enumerate(tags) if p!='I') for tags in cartesian('IXYZ',repeat=n) if any(p!='I' for p in tags)]
        choose=rng.choice(len(pool),size=int(rng.integers(2,min(12,len(pool))+1)),replace=False);toys.append((f'random{it}',n,[pool[int(i)] for i in choose]))
    records=[]
    for name,n,words in toys:
        a,b,o=endpoints(words);U=len(words);coeffs=rng.uniform(-.7,.7,U);matrices=[matrix(w,n) for w in words]
        for r in (1,2,5):
            unfused=[(i,1) for _ in range(r) for i in list(range(U))+list(reversed(range(U)))];fused=fused_pulses(U,a,b,r)
            error=float(np.linalg.norm(product(unfused,matrices,coeffs,r)-product(fused,matrices,coeffs,r)))
            c.require(error<3e-12,'Dense boundary-fusion oracle failed')
            records.append({'case':name,'r':r,'U':U,'prefix':a,'suffix':b,'overlap':o,'fused_rotations':len(fused),'unitary_frobenius_difference':error})
    c.require(records[0]['prefix']==2 and records[0]['suffix']==2 and records[0]['overlap']==1,'Overlap witness not tested')
    return {'dense_products':len(records),'maximum_unitary_frobenius_difference':max(r['unitary_frobenius_difference'] for r in records),'records':records}

def main():
    begin=time.monotonic();oracle=dense_checks();plans=[];newrows=[];identities=[];pending=[];inputs={}
    ownrows=json.loads((HERE/'results.json').read_text());inputs['optimized/results.json']=c.sha(HERE/'results.json')
    for case in c.CASES:
        path=HERE/'designs'/f'{case}.design.json';d=json.loads(path.read_text());inputs['optimized/designs/'+path.name]=c.sha(path)
        refined=[(w,v) for b in d['blocks']['refined'] for w,v in c.block(b).items()];complete=c.decode_schedule(d['methods']['strang_class_lex']['schedule']);c.require(refined==complete,'Refined label and complete-class label lists differ')
        identities.append({'case':case,'literal_flattened_refined_label_equals_complete_class_label':True,'schedule_sha256':d['methods']['strang_class_lex']['schedule_sha256'],'meaning':'Their Strang unitary is identical for every same r; contiguous commuting diagonal factors combine exactly.'})
        sources=[('optimized',d,[r for r in ownrows if r['case']==case])]
        vp=HERE.parent/'published/v2/designs'/f'{case}.json'
        if vp.exists():
            vd=json.loads(vp.read_text());inputs['published/v2/designs/'+vp.name]=c.sha(vp);vrows=[]
            for p in sorted((HERE.parent/'published/v2/results').glob(case+'.*.b*.json')):
                row=json.loads(p.read_text())
                if row.get('status')=='complete':vrows.append(row);inputs['published/v2/results/'+p.name]=c.sha(p)
            sources.append(('published/v2',vd,vrows))
        else:pending.append(case)
        for campaign,source,rows in sources:
            for name,method in source['methods'].items():
                if method['formula_order']!=2:continue
                schedule=c.decode_schedule(method['schedule']);words=[w for w,v in schedule];U=len(words);a,b,o=endpoints(words);slope=2*U-a-b
                lengths=method.get('unit_lengths',[1]*U);L1,Llast=lengths[0],lengths[-1];oldr=(100*U-L1)//(2*U-L1-Llast);newr=(100*U-a)//slope if slope else 1
                c.require(a>=L1 and b>=Llast,'Uniform prefix shorter than declared unit')
                plans.append({'case':case,'campaign':campaign,'ordering':name,'U':U,'declared_first_words':L1,'declared_last_words':Llast,'maximal_commuting_prefix':a,'maximal_commuting_suffix':b,'overlap':o,'old_steps_at_budget100U':oldr,'possible_uniform_steps_at_budget100U':newr,'new_step_count_was_not_simulated_by_this_analysis':True,'schedule_sha256':method['schedule_sha256']})
                for row in rows:
                    if row['ordering']!=name or row.get('status')!='complete':continue
                    r=row['steps'];count=fused_count(U,a,b,r);c.require(count<=row['pauli_rotations'],'Uniform compilation increased cost')
                    newrows.append({'case':case,'campaign':campaign,'ordering':name,'original_budget_units':row['budget_units'],'steps':r,'one_minus_overlap':row['one_minus_overlap'],'original_emitted_rotations':row['pauli_rotations'],'uniform_endpoint_rotations':count,'state_file':row.get('state_file'),'state_sha256':row.get('state_sha256'),'schedule_sha256':method['schedule_sha256'],'inference':'Existing state/error reused under proved exact unitary equivalence; no molecular state was propagated for this cost sensitivity.'})
    thresholds={};baseline='strang_class_refined_response';names=sorted(set(p['ordering'] for p in plans))
    for eps in (1e-4,1e-6,1e-8):
        work={}
        for case in c.CASES:
            for name in names:
                eligible=[r['uniform_endpoint_rotations'] for r in newrows if r['case']==case and r['ordering']==name and r['one_minus_overlap']<=eps]
                work[case,name]=min(eligible) if eligible else None
        report={}
        for name in names:
            if name==baseline:continue
            paired=[];counts={'conserving_win':0,'tie':0,'conserving_loss':0}
            for case in c.CASES:
                aa=work[case,baseline];bb=work[case,name]
                if aa is None or bb is None:continue
                counts['tie' if aa==bb else 'conserving_win' if aa<bb else 'conserving_loss']+=1;paired.append({'case':case,'conserving_rotations':aa,'comparator_rotations':bb,'comparator_over_conserving':bb/aa})
            report[name]={'paired_cases':len(paired),'outcomes':counts,'median_comparator_over_conserving':statistics.median(p['comparator_over_conserving'] for p in paired) if paired else None,'pairs':paired}
        thresholds[str(eps)]=report
    result={'scope':'Cost-only sensitivity on existing symmetric-product states, with uniform maximal commuting endpoint fusion. No new molecular propagations and no newly saturated-budget error claims.','plans':plans,'cost_rows':newrows,'threshold_work':thresholds,'refined_label_correspondence':identities,'dense_checks':oracle,'pending_v2_plans':pending,'input_sha256':inputs,'source_sha256':c.sha(__file__),'seconds':time.monotonic()-begin}
    c.write_atomic(HERE/'boundary_sensitivity.json',result)
    lines=['# Uniform boundary-fusion sensitivity','','This postprocessing retains every measured state, step count and error. It changes only the rotation accounting using a single compilation rule for all symmetric schedules. It does **not** create observations at the larger step counts that this compilation might permit under 100U.','','For U ordered words, let a be the maximal pairwise-commuting prefix length and b the corresponding suffix length. Overlap c=max(0,a+b-U) contains words commuting with the entire list: each overlap word is in both sets, whose union covers all words. Extract those c words into one full-time rotation each. The remaining prefix and suffix are disjoint and commute internally. Fuse the suffix at each palindrome center and the prefix across successive steps. The cost is r(2U-a-b)+a, including the overlap case. If all words commute, a=b=U and the entire product costs U rotations, independent of r.','','For a complete diagonal class followed by complete nonzero flip classes, refining the diagonal class into consecutive singleton words leaves the flattened list unchanged. Since those diagonal terms commute, their forward and reverse half-exponentials combine into the corresponding half-exponential of D. The complete-class label Strang formula and refined-label Strang formula therefore define exactly the same unitary at every same r. All 14 actual flattened lists agree exactly.','','Dense tests include X0,Z1,Z0, where AB and BC commute but AC do not, plus an all-commuting list and 80 random lists. All integer pulse durations are checked before comparing full unitaries.','',f"Checked {oracle['dense_products']} dense products; maximum Frobenius discrepancy {oracle['maximum_unitary_frobenius_difference']:.3g}.",'','## Feasible step-count sensitivity at 100U','','| Schedule | Plans | Cases permitting more steps | Largest old→possible r |','|---|---:|---:|---|']
    for name in names:
        pp=[p for p in plans if p['ordering']==name];changed=[p for p in pp if p['possible_uniform_steps_at_budget100U']>p['old_steps_at_budget100U']];big=max(pp,key=lambda p:p['possible_uniform_steps_at_budget100U']-p['old_steps_at_budget100U']);lines.append(f"| {name} | {len(pp)} | {len(changed)} | {big['case']}: {big['old_steps_at_budget100U']}→{big['possible_uniform_steps_at_budget100U']} |")
    lines+=['','These possible extra steps were not simulated by this sensitivity. The following table uses only already-tested r values with recompiled costs. A ratio above 1 means the comparator needs more tested rotations than conserving-order Strang.','','| Error target | Comparator | Paired cases | Conserving W/T/L | Median comparator/conserving work |','|---|---|---:|---|---:|']
    for eps,report in thresholds.items():
        for name,r in report.items():
            cc=r['outcomes'];ratio=r['median_comparator_over_conserving'];lines.append(f"| {eps} | {name} | {r['paired_cases']} | {cc['conserving_win']}/{cc['tie']}/{cc['conserving_loss']} | "+('unavailable' if ratio is None else f'{ratio:.6g}')+' |')
    lines+=['',f"Pending V2 plans: {', '.join(pending) if pending else 'none'}. V2 state coverage may still be incomplete until that campaign finishes; every table reports its paired count. Original emitted costs remain unchanged in the raw records. This is not a full commuting-gate optimizer or a hardware gate-count estimate.",'']
    (HERE/'BOUNDARY_SENSITIVITY.md').write_text('\n'.join(lines));print(json.dumps({'plans':len(plans),'cost_rows':len(newrows),'pending_v2':pending,'dense_max':oracle['maximum_unitary_frobenius_difference'],'thresholds':thresholds},indent=2),flush=True)
if __name__=='__main__':main()
