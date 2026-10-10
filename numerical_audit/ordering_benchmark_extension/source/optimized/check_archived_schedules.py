"""Check literal old schedule reproduction, including within-unit word order."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import common as c

def flatten(raw,order):
    bs=[c.block(b) for b in raw];return [(w,v) for i in order for w,v in bs[i].items()]
def main():
    checks=[]
    for case in c.CASES:
        p=HERE/'designs'/f'{case}.design.json';d=json.loads(p.read_text());n=c.normalized_case(case);raw=n['design'];expected={}
        if case in c.PROS_CASES:
            for old,new in (('class_refined_response','class_refined_response'),('jw_signed','strang_jw_signed'),('jw_magnitude','strang_jw_magnitude'),('class_lex','strang_class_lex')):expected[new]=c.decode_schedule(raw['flat_schedules'][old])
        else:
            old=json.loads((c.MECH/'diagonal_refinement/results'/f'{case}.design.json').read_text());expected['class_refined_response']=flatten(old['blocks']['class'],old['orders']['class_refined_response'])
            for name in ('class_refined_J0','class_pauli_JD'):
                old=json.loads((c.MECH/'objective_ablation/results'/name/f'{case}.design.json').read_text());expected[name]=flatten(old['blocks']['class'],old['orders'][name])
            for name in ('jw_signed','jw_magnitude'):expected['strang_'+name]=[(w,n['total'][w]) for w in n['jw_orders'][name]]
            expected['strang_class_lex']=flatten(raw['blocks']['class'],raw['orders']['class_lex'])
        for name,schedule in expected.items():
            actual=c.decode_schedule(d['methods'][name]['schedule']);c.require(actual==schedule,'Changed archived literal schedule: '+case+' '+name)
            checks.append({'case':case,'ordering':name,'terms':len(actual),'schedule_sha256':d['methods'][name]['schedule_sha256'],'design_sha256':c.sha(p),'exact_word_and_coefficient_sequence':True})
    c.require(len(checks)==72,'Wrong archived sequence-check count')
    out={'cases':c.CASES,'checks':checks,'count':len(checks),'source_sha256':c.sha(__file__)};c.write_atomic(HERE/'archived_schedule_checks.json',out);print('ARCHIVED_LITERAL_SCHEDULES',len(checks),flush=True)
if __name__=='__main__':main()
