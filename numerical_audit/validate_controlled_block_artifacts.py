"""Independent artifact validation; run after campaign completion."""
from pathlib import Path
import hashlib,json,math,sys

run=Path(sys.argv[1]); root=run.parents[2]
def digest(o):
 return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=json.loads((run/'manifest.json').read_text())
completion=json.loads((run/'completion.json').read_text())
presentation=json.loads((run/'presentation_provenance.json').read_text()) if (run/'presentation_provenance.json').exists() else None
rows=json.loads((run/'results.json').read_text())
assert completion['status']=='complete'
assert len(rows)==completion['rows']==manifest['expected_rows']==728
for rel,h in manifest['source_hashes'].items():
 assert sha(run/'source_snapshot'/rel)==h,rel
 if presentation and rel == presentation['source_relative_path']:
  assert rel == 'numerical_audit/plot_controlled_block_orders.py'
  assert presentation['role'] == 'presentation_only'
  assert presentation['original_source_sha256'] == h
  assert sha(root/rel) == presentation['revised_source_sha256']
  assert sha(run/'presentation_source_snapshot'/rel) == presentation['revised_source_sha256']
 else:
  assert sha(root/rel)==h,rel
if presentation:
 for rel,h in presentation['unchanged_numerical_artifacts'].items():
  assert sha(run/rel)==h,rel
assert len(set((r['case'],r['ordering'],r['T'],r['steps']) for r in rows))==len(rows)
for case in manifest['cases']:
 design=json.loads((run/f'{case}.design.json').read_text())
 blocks=design['blocks']; reference=design['pauli_order_key_reference']
 coefficient={json.dumps(key): (real,imag) for key,real,imag in design['identity_free_coefficients']}
 membership=[[block['parent_index'],block['pauli_keys']] for block in blocks]
 assert digest(membership)==design['membership_sha256']
 assert len({json.dumps(key) for block in blocks for key in block['pauli_keys']})==len(reference)
 for name,order in design['pauli_orders'].items():
  assert sorted(order)==list(range(len(reference)))
  keys=[reference[i] for i in order]
  schedule=[[key,*coefficient[json.dumps(key)]] for key in keys]
  assert digest(schedule)==design['schedule_sha256'][name]
  if name.startswith('block_'):
   expected=[key for index in design['block_orders'][name] for key in blocks[index]['pauli_keys']]
   assert keys==expected
 for row in [r for r in rows if r['case']==case]:
  assert row['hamiltonian_sha256']==design['hamiltonian_sha256']
  assert row['schedule_sha256']==design['schedule_sha256'][row['ordering']]
  assert row['nominal_exponential_count']==row['steps']*design['pauli_terms']
  assert abs(row['infidelity']-row['one_minus_overlap']*(2-row['one_minus_overlap']))<1e-18
summary=json.loads((run/'summary.json').read_text())
comparisons=[]
for case in manifest['cases']:
 for t in manifest['times']:
  group={r['ordering']:r['one_minus_overlap'] for r in rows if r['case']==case and r['T']==t and r['steps']==100}
  e=group['block_signed']
  vals=[group[f'block_random_{seed}'] for seed in manifest['random_seeds']]
  better=sum(value<e-manifest['tie_absolute'] for value in vals)
  tie=sum(abs(value-e)<=manifest['tie_absolute'] for value in vals)
  claimed=next(c for c in summary['comparisons'] if c['case']==case and c['T']==t and c['steps']==100)
  assert claimed['signed_midrank_among_signed_and_random']==1+better+tie/2
  comparisons.append({'case':case,'T':t,'signed_rank':1+better+tie/2})
report={'status':'passed','rows':len(rows),'source_and_snapshot_hashes_checked':len(manifest['source_hashes']),'case_designs_checked':len(manifest['cases']),'independent_rank_checks':len(comparisons),'maximum_norm_defect':max(max(abs(r['exact_norm']-1),abs(r['approx_norm']-1)) for r in rows),'validator_sha256':sha(__file__),'presentation_change_validated':presentation is not None}
output_name = sys.argv[2] if len(sys.argv)>2 else 'independent_artifact_validation.json'
with (run/output_name).open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(report,indent=2))
