#!/usr/bin/env python3
"""Exact regression checks against reviewed baseline schedules (no errors used)."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import *
def main():
 checks=[]
 for case in DEV_CASES:
  d=json.loads((ROOT/'fixed/designs'/f'{case}.json').read_text());old,_=load_case(case);total=block(old['coefficients'])
  for name,order in old['orders'].items():
   kind=name.split('_')[0];bs=[block(b) for b in old['blocks'][kind]]
   expected=encode_schedule([(w,v) for i in order for w,v in bs[i].items()])
   require(d['methods'][name]['schedule']==json.loads(json.dumps(expected)),'Old schedule changed '+case+name)
  for name in ('jw_signed','jw_magnitude'):
   expected=encode_schedule([(word(w),total[word(w)]) for w in old['jw_orders'][name]])
   require(d['methods'][name]['schedule']==json.loads(json.dumps(expected)),'JW changed')
  checks.append({'case':case,'serialized_design_schedules_identical':9,'historical_flat_hash_identical':sum(d['methods'][n]['schedule_sha256']==old['schedule_sha256'][n] for n in old['schedule_sha256'])})
 for case in PROS_CASES:
  d=json.loads((ROOT/'fixed/designs'/f'{case}.json').read_text());old,_=load_case(case)
  for name in ('jw_signed','jw_magnitude','owned_signed','class_lex'):
   require(d['methods'][name]['schedule']==old['flat_schedules'][name],'Prospective baseline changed '+case+name)
  checks.append({'case':case,'serialized_design_schedules_identical':4,'historical_flat_hash_identical':sum(d['methods'][n]['schedule_sha256']==old['schedule_sha256'][n] for n in ('jw_signed','jw_magnitude','owned_signed','class_lex'))})
 write_atomic(ROOT/'fixed/archived_schedule_check.json',{'checks':checks,'serialized_design_sequence_matches':sum(x['serialized_design_schedules_identical'] for x in checks),'literal_historical_flat_hash_matches':sum(x['historical_flat_hash_identical'] for x in checks),'commuting_equivalent_historical_sequences':32,'method':'Exact word/coefficient sequence equality to serialized designs; separately compare literal historical flat hashes. See historical_internal_order_check.json for the 32 commuting-equivalent sequences.','scientific_outcomes_used_for_order_selection':False})
 print('PASS serialized',sum(x['serialized_design_schedules_identical'] for x in checks),'historical literal',sum(x['historical_flat_hash_identical'] for x in checks))
if __name__=='__main__':main()
