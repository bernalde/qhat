#!/usr/bin/env python3
"""Frozen deterministic and fixed-random molecular ordering controls.

Preparation is independent of circuit outcomes. Each design is atomic and all
14 designs are hashed before any evolution. Every run checkpoints separately.
"""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,hashlib,json,math,random,resource,sys,time
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from common import *
from qhat.numerical_audit.followup_classes import exact_blocks as eb
sys.path.insert(0,str(HERE.parent/'published'))
import algorithms as published


def dense(w,n):
    a=['I']*n
    for q,p in w:a[q]=p
    return ''.join(a)
def flatten(bs,order):return [(w,v) for i in order for w,v in bs[i].items()]
def freeze_method(schedule,metadata):
    raw=encode_schedule(schedule)
    return {'schedule':raw,'metadata':metadata,'terms':len(schedule),'schedule_sha256':object_sha(raw),'formula_order':1}
def structural_fingerprint(d):
    nums=[abs(int(v)) for b in d['blocks']['parent'] for _,v in b]+[abs(int(v)) for _,v in d['coefficients']]
    gcd=math.gcd(*nums)
    payload={'parent':[[[w,str(int(v)//gcd)] for w,v in b] for b in d['blocks']['parent']], 'target':[[w,str(int(v)//gcd)] for w,v in d['coefficients']]}
    return object_sha(payload)
def pauli_color(total,n):
    keys=sorted(total,key=lambda w:dense(w,n));adj=published.adjacency([published.masks(w) for w in keys]);colors={};groups=[];masks=[]
    for i in sorted(range(len(keys)),key=lambda i:(-adj[i].bit_count(),i)):
        c=0
        while c<len(masks) and masks[c]&adj[i]:c+=1
        if c==len(masks):masks.append(0);groups.append([])
        masks[c]|=1<<i;colors[i]=c
    for i in range(len(keys)):groups[colors[i]].append(i)
    return [(keys[i],total[keys[i]]) for group in groups for i in group],{'colors':len(groups),'algorithm':'largest-degree greedy; dense IXYZ ties; dense IXYZ order within color'}
def class_candidates(bs,comms):
    # Fixed first class i sends [G_i,G_j] to mutually distinct flip masks as j
    # varies. The suffix coefficient-l1 objective is therefore sum edge weights.
    m=len(bs);weights={(i,j):sum(map(abs,c.values())) for (i,j),c in comms.items()};degree=[0]*m;neighbors=[[] for _ in bs]
    for (i,j),v in weights.items():degree[i]+=v;degree[j]+=v;neighbors[i].append((j,v));neighbors[j].append((i,v))
    todo=set(range(m));suffix=[]
    while todo:
        i=min(todo,key=lambda j:(degree[j],j));todo.remove(i);suffix.append(i)
        for j,v in neighbors[i]:degree[j]-=v
    order=list(range(m));current={}
    for c in comms.values():eb.add_into(current,c)
    value=sum(map(abs,current.values()));initial=value;accepted=[]
    for it in range(2):
        for at in range(m-1):
            i,j=order[at:at+2];c=comms.get((min(i,j),max(i,j)),{});factor=-2 if i<j else 2
            trial=value+sum(abs(current.get(k,0)+factor*v)-abs(current.get(k,0)) for k,v in c.items())
            if trial<value:
                eb.add_into(current,c,factor);order[at:at+2]=j,i;accepted.append([it,at,i,j]);value=trial
    # Different implementation reconstructs the final polynomial.
    check={};pos={i:j for j,i in enumerate(order)}
    for (i,j),c in comms.items():eb.add_into(check,c,1 if pos[i]<pos[j] else -1)
    require(check==current and sum(map(abs,check.values()))==value,'BCH incremental update failed')
    return suffix,order,{'initial_l1_integer':str(initial),'final_l1_integer':str(value),'accepted_swaps':accepted}

def prepare_case(case):
    p=HERE/'designs'/f'{case}.json'
    if p.exists():return json.loads(p.read_text())
    started=time.monotonic();data=normalized_case(case);d=data['design'];bs=data['blocks'];total=data['total'];n=d['n_qubits'];methods={};orders={};diagnostics={}
    def add(name,schedule,**metadata):
        acc={}
        for w,v in schedule:eb.add_into(acc,{w:v})
        require(acc==total,'Wrong target in '+name)
        methods[name]=freeze_method(schedule,metadata)
    def addorder(name,kind,order,**metadata):
        require(sorted(order)==list(range(len(bs[kind]))),'Invalid permutation '+name)
        orders[name]=order
        schedule=decode_schedule(d['flat_schedules'][name]) if case in PROS_CASES and name=='owned_signed' else flatten(bs[kind],order)
        if case in PROS_CASES and name=='owned_signed':metadata['internal_word_order']='literal archived flat schedule within each commuting owned block'
        add(name,schedule,unit_partition=kind,unit_order=order,**metadata)
    for name in ('jw_signed','jw_magnitude'):add(name,[(w,total[w]) for w in data['jw_orders'][name]],tie_policy='literal archived floating-coefficient order; exact integer coefficients used in evolution')
    add('jw_lex',[(w,total[w]) for w in sorted(total,key=lambda w:dense(w,n))],tie_policy='dense IXYZ lexical')
    schedule,meta=pauli_color(total,n);add('jw_color',schedule,**meta)
    pm=sorted(range(len(bs['parent'])),key=lambda i:(-abs(int(d['parent_metadata'][i]['weight'])),i))
    addorder('parent_signed','parent',list(range(len(bs['parent']))))
    addorder('parent_magnitude','parent',pm,tie_policy='descending absolute canonical parent weight; original parent index')
    addorder('owned_signed','owned',list(range(len(bs['owned']))))
    owner={}
    for i,b in enumerate(bs['parent']):
        for w in b:owner.setdefault(w,i)
    owned_owners=[]
    for b in bs['owned']:
        oo={owner[w] for w in b};require(len(oo)==1,'Ambiguous owned block');owned_owners.append(next(iter(oo)))
    rank={i:j for j,i in enumerate(pm)}
    addorder('owned_parent_magnitude','owned',sorted(range(len(bs['owned'])),key=lambda i:(rank[owned_owners[i]],i)),ownership='frozen signed-parent ownership, only existing blocks reordered')
    addorder('owned_l1','owned',sorted(range(len(bs['owned'])),key=lambda i:(-sum(map(abs,bs['owned'][i].values())),i)))
    addorder('owned_reverse','owned',list(reversed(range(len(bs['owned'])))))
    addorder('class_lex','class',list(range(len(bs['class']))))
    addorder('class_magnitude','class',sorted(range(len(bs['class'])),key=lambda i:(-sum(map(abs,bs['class'][i].values())),d['class_flip_masks'][i])))
    addorder('refined_magnitude','refined',sorted(range(len(bs['refined'])),key=lambda i:(-sum(map(abs,bs['refined'][i].values())),data['refined_flip_masks'][i],i)))
    fingerprint=structural_fingerprint(d);cache=HERE/'structure_cache'/f'{fingerprint}.json'
    if cache.exists():
        cached=json.loads(cache.read_text());diagnostics['structural_cache']=fingerprint
        for kind in ('parent','owned','class'):addorder(kind+'_color',kind,cached['orders'][kind+'_color'],algorithm='exact operator graph; largest-degree greedy; original index ties')
        for name in ('class_suffix_l1','class_bch_l1'):addorder(name,'class',cached['orders'][name],algorithm=cached['algorithm'][name])
    else:
        comms={}
        for kind in ('parent','owned','class'):
            t=time.monotonic()
            if case in DEV_CASES and kind+'_color' in d['orders']:
                order=d['orders'][kind+'_color'];info=d['graphs'][kind];new_comms={}
                if kind=='class':
                    for i in range(len(bs[kind])):
                        for j in range(i+1,len(bs[kind])):
                            c=eb.commutator(bs[kind][i],bs[kind][j])
                            if c:new_comms[i,j]=c
            else:order,info,new_comms=eb.graph_and_color(bs[kind],keep_commutators=kind=='class')
            diagnostics[kind+'_graph_seconds']=time.monotonic()-t
            addorder(kind+'_color',kind,order,algorithm='exact operator graph; largest-degree greedy; original index ties')
            if kind=='class':comms=new_comms
            print('GRAPH',case,kind,round(time.monotonic()-t,3),flush=True)
        suffix,bch,bchinfo=class_candidates(bs['class'],comms)
        if case in DEV_CASES:require(suffix==d['orders']['class_suffix_l1'],'Suffix differs from archived development order')
        addorder('class_suffix_l1','class',suffix,algorithm='exact weighted minimum-degree deletion for distinct complete flip classes; no strict-decrease adjacent swap exists')
        addorder('class_bch_l1','class',bch,algorithm='two exact strict-improvement adjacent passes on global BCH coefficient-l1, starting class_lex')
        diagnostics['bch']=bchinfo
        if case in DEV_CASES:
            old=json.loads((QHAT/'numerical_audit/followup_classes/results/bch_supplement'/f'{case}.bch_design.json').read_text())
            require(bch==old['order'],'BCH order differs from archived development supplement')
        write_atomic(cache,{'fingerprint':fingerprint,'orders':{k:orders[k] for k in ('parent_color','owned_color','class_color','class_suffix_l1','class_bch_l1')},'algorithm':{k:methods[k]['metadata']['algorithm'] for k in ('class_suffix_l1','class_bch_l1')}})
    for kind,label in (('owned','owned'),('refined','refined'),('pauli','pauli')):
        units=bs[kind] if kind!='pauli' else [{w:total[w]} for w in sorted(total,key=lambda w:dense(w,n))]
        for rep in range(20):
            seed=2026100900+rep;order=list(range(len(units)));random.Random(seed).shuffle(order)
            add('random_'+label+'_'+str(rep).zfill(2),flatten(units,order),unit_partition=kind,unit_order=order,seed=seed,randomization='one fixed permutation repeated at every time step; no accuracy-based selection')
    require(len(methods)==78,'Expected 18 basic +60 random')
    t=time.monotonic()
    for name,value in published.build_published_schedules(d).items():add(name,value['schedule'],**value['metadata'])
    diagnostics['published_seconds']=time.monotonic()-t
    require(len(methods)==82,'Expected total82 methods')
    d['extension_methods']=list(methods);d['methods']=methods;d['extension_preparation']=diagnostics;d['extension_preparation']['seconds']=time.monotonic()-started;d['extension_panel']='development' if case in DEV_CASES else 'prospective_followup';d['extension_protocol']='fixed22 +20 each owned/refined/Pauli random; T1; matched emitted rotation budgets'
    write_atomic(p,d);print('PREPARED',case,len(methods),round(time.monotonic()-started,2),flush=True);return d

def freeze():
    designs={case:sha(HERE/'designs'/f'{case}.json') for case in CASES}
    source={str(p.relative_to(HERE.parent)):sha(p) for p in (Path(__file__),HERE.parent/'common.py',HERE.parent/'published/algorithms.py',HERE/'PROTOCOL.md')}
    record={'cases':CASES,'design_sha256':designs,'source_sha256':source,'expected_rows':len(CASES)*(22*len(BUDGETS)+60),'budget_units':BUDGETS,'random_replicates':20,'status':'all schedules frozen before first new error'}
    target=HERE/'freeze.json'
    if target.exists():require(json.loads(target.read_text())==record,'Freeze changed')
    else:write_atomic(target,record)
    return record

def run_case(case):
    import numpy as np
    from qhat.numerical_audit.run_frozen_audit import evolve,physical_metrics
    frozen=json.loads((HERE/'freeze.json').read_text());path=HERE/'designs'/f'{case}.json'
    require(sha(path)==frozen['design_sha256'][case],'Design mutated')
    for p,h in frozen['source_sha256'].items():require(sha(HERE.parent/p)==h,'Source mutated '+p)
    d=json.loads(path.read_text());raw,states=load_case(case);require(d['extension_input_provenance']==raw['extension_input_provenance'],'Input changed');compiled=compile_case(d,states)
    unit=len(d['coefficients']);out=HERE/'results';out.mkdir(exist_ok=True);rows=[]
    for name,method in d['methods'].items():
        schedule=decode_schedule(method['schedule']);terms=[compiled['term'](w,v) for w,v in schedule]
        for b in ((100,) if name.startswith('random_') else BUDGETS):
            rowpath=out/f'{case}.{name}.b{b}.json';statepath=out/f'{case}.{name}.b{b}.npz';steps=unit*b//len(terms)
            if rowpath.exists():
                row=json.loads(rowpath.read_text());require(row['extension_design_sha256']==frozen['design_sha256'][case] and row['schedule_sha256']==method['schedule_sha256'],'Checkpoint mismatch')
                if row['status']=='complete':require(statepath.exists() and sha(statepath)==row['state_sha256'],'State checkpoint corrupt')
                rows.append(row);continue
            row={'case':case,'panel':d['extension_panel'],'ordering':name,'T':1.,'formula_order':1,'budget_units':b,'rotation_budget':unit*b,'steps':steps,'role':'matched_rotation_budget','extension_design_sha256':frozen['design_sha256'][case],'schedule_sha256':method['schedule_sha256'],'pauli_rotations':steps*len(terms)}
            if steps==0:row['status']='infeasible';write_atomic(rowpath,row);rows.append(row);continue
            t=time.monotonic();actual=evolve(compiled['initial'],terms,1.,steps)
            row.update(status='complete',seconds=time.monotonic()-t,state_file=statepath.name,state_key='actual',**physical_metrics(compiled['exact'],actual,compiled['number_mask'],compiled['spin_mask']))
            tmp=statepath.with_suffix('.partial.npz');np.savez_compressed(tmp,actual=actual);os.replace(tmp,statepath);row['state_sha256']=sha(statepath);write_atomic(rowpath,row);rows.append(row)
            print('EVOLVE',case,name,b,steps,row['one_minus_overlap'],round(row['seconds'],3),flush=True)
    require(len(rows)==214,'Incomplete case rows');write_atomic(out/f'{case}.completion.json',{'rows':len(rows),'status':'complete'});return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare','freeze','run','collect'));ap.add_argument('--cases',nargs='*',default=CASES);a=ap.parse_args();resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    if a.action=='prepare':
        for case in a.cases:prepare_case(case)
    elif a.action=='freeze':freeze()
    elif a.action=='run':
        for case in a.cases:run_case(case)
    else:
        rows=[]
        for case in CASES:
            require((HERE/'results'/f'{case}.completion.json').exists(),'Incomplete case '+case)
            for p in sorted((HERE/'results').glob(case+'.*.b*.json')):rows.append(json.loads(p.read_text()))
        require(len(rows)==2996,'Incomplete aggregate');write_atomic(HERE/'results.json',rows);write_atomic(HERE/'completion.json',{'status':'complete','rows':len(rows),'complete':sum(r['status']=='complete' for r in rows),'infeasible':sum(r['status']=='infeasible' for r in rows)})
if __name__=='__main__':main()
