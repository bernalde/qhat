#!/usr/bin/env python3
"""Package a complete, hash-bound extension and its exact QHAT snapshot."""
import argparse,hashlib,json,subprocess,zipfile
from pathlib import Path
from verify_archive import require,file_hash,verify

CASES=[m+'_s'+s for m in ('CH2O','H2S','HCN','SiH4') for s in ('1.00','1.25')]+[m+'_coeff'+s for m in ('f2','li2','nh3') for s in ('1p00','1p25')]
CAMPAIGNS=(('fixed','fixed',2996),('optimized','optimized',728),('published/extra','extra',392),('published/v2','v2',196))

def read_json(path):
    return json.loads(path.read_text())

def contained(root,name):
    root=root.resolve();path=(root/name).resolve()
    require(path.is_relative_to(root),'Path escapes archive root: '+str(name))
    return path

def check_hash(path,expected):
    require(path.is_file(),'Missing bound file: '+str(path))
    actual=file_hash(path)
    require(actual==expected,'Bound file changed: '+str(path))
    return actual

def check_bindings(root,bindings):
    for name,digest in bindings.items():check_hash(contained(root,name),digest)

def campaign_bindings(extension,campaign_name,audit_name,cases,expected_rows):
    """Revalidate recorded byte bindings only; never rerun propagation."""
    root=extension/campaign_name;audit=extension/'independent'/audit_name
    summary=read_json(audit/'summary.json')
    require(summary['status']=='passed' and summary['case_count']==len(cases) and summary['rows']==expected_rows,'Incomplete independent audit '+audit_name)
    require(len(summary['cases'])==len(cases) and set(summary['cases'])==set(cases),'Wrong independent case inventory '+audit_name)
    freezes=[p for p in (root/'freeze.json',root/'ALL_SCHEDULES_FROZEN.json') if p.is_file()]
    require(len(freezes)==1,'Missing or ambiguous campaign freeze '+campaign_name)
    check_hash(freezes[0],summary['freeze_sha256']);frozen=read_json(freezes[0])
    require(len(frozen['cases'])==len(cases) and set(frozen['cases'])==set(cases),'Wrong frozen case inventory')
    require(set(frozen['design_sha256'])==set(cases),'Wrong frozen design inventory')
    check_bindings(extension,frozen['source_sha256'])
    rows={};provenance={};completed=0
    for case in cases:
        designs=[p for p in (root/'designs'/(case+'.json'),root/'designs'/(case+'.design.json')) if p.is_file()]
        require(len(designs)==1,'Missing or ambiguous design '+case);dp=designs[0]
        design_hash=check_hash(dp,frozen['design_sha256'][case]);design=read_json(dp)
        provenance[case]=design['extension_input_provenance']
        record=read_json(audit/(case+'.json'))
        require(record['case']==case and record['status']=='passed' and record['design_sha256']==design_hash,'Stale design audit '+case)
        check_hash(extension/'independent/check_campaign.py',record['source_sha256'])
        check_hash(extension/'independent/core.py',record['core_sha256'])
        paths=list((root/'results').glob(case+'.*.b*.json'))+list((root/'rows'/case).glob('*.json'))
        current={}
        for path in paths:
            row=read_json(path);key=(row['ordering'],row['budget_units'])
            require(row['case']==case and key not in current,'Duplicate or misidentified row '+str(path));current[key]=(path,row)
        require(len(record['checks'])==record['rows']==len(current),'Changed audit row inventory '+case)
        seen=set()
        for checked in record['checks']:
            key=(checked['ordering'],checked['budget_units'])
            require(key in current and key not in seen,'Missing or duplicate audited row '+case);seen.add(key)
            path,row=current[key];digest=check_hash(path,checked['row_sha256'])
            require(row['extension_design_sha256']==design_hash,'Row/design binding changed')
            require(row['status'] in ('complete','infeasible'),'Unexpected row status')
            require(checked['status']==('passed' if row['status']=='complete' else 'infeasible'),'Audit/result status mismatch')
            if row['status']=='complete':
                candidates=[contained(root,row['state_file']),contained(root,str(path.parent.relative_to(root)/row['state_file']))]
                states=[p for p in dict.fromkeys(candidates) if p.is_file()]
                require(len(states)==1,'Missing or ambiguous saved state '+str(path))
                require(checked.get('state_sha256')==row['state_sha256'],'Audit/state binding changed')
                check_hash(states[0],checked['state_sha256']);completed+=1
            rows[path.relative_to(extension).as_posix()]=digest
    all_paths=set((root/'rows').rglob('*.json')) if (root/'rows').exists() else set((root/'results').glob('*.b*.json'))
    require(set(rows)=={p.relative_to(extension).as_posix() for p in all_paths},'Untracked campaign rows '+campaign_name)
    require(len(rows)==expected_rows,'Wrong audited row count '+campaign_name)
    return rows,provenance,completed

def input_provenances(qhat,cases):
    values={};heldout=qhat/'numerical_audit/ordering_mechanism/heldout/heldout_validation.zip'
    archive_hash=file_hash(heldout)
    with zipfile.ZipFile(heldout) as z:
        for case in cases:
            if '_s' in case:
                dp=Path('numerical_audit/followup_classes/results/comparison')/(case+'.design.json');sp=dp.with_name(case+'.states.npz')
                values[case]={'design_sha256':file_hash(qhat/dp),'states_sha256':file_hash(qhat/sp),'design_source':str(dp),'state_source':str(sp)}
            else:
                dp='results/'+case+'.design.json';sp='results/'+case+'.states.npz'
                def member_hash(name):
                    h=hashlib.sha256()
                    with z.open(name) as f:
                        for data in iter(lambda:f.read(4*1024**2),b''):h.update(data)
                    return h.hexdigest()
                values[case]={'design_sha256':member_hash(dp),'states_sha256':member_hash(sp),'archive_sha256':archive_hash,'design_source':dp,'state_source':sp}
    return values

def preflight(extension,qhat):
    extension=extension.resolve();qhat=qhat.resolve();summary=read_json(extension/'analysis/summary.json')
    require(summary['rows']==4312 and summary['status']=='complete','Incomplete extension')
    check_hash(extension/'analysis/summarize.py',summary['source_sha256'])
    rows={};all_provenance=[];completed=0
    for campaign,audit,count in CAMPAIGNS:
        observed,provenance,n=campaign_bindings(extension,campaign,audit,CASES,count)
        require(not set(rows).intersection(observed),'Duplicate cross-campaign row');rows.update(observed);all_provenance.append(provenance);completed+=n
    require(summary['input_sha256']==rows,'Analysis input inventory/hashes differ from audited rows')
    require(summary['completed_evolutions']==completed and summary['infeasible']==len(rows)-completed,'Analysis status counts differ')
    expected=input_provenances(qhat,CASES)
    for provenance in all_provenance:require(provenance==expected,'QHAT inputs differ from frozen designs')
    for case in CASES:
        record=read_json(extension/'independent/references'/(case+'.json'))
        require(record['case']==case and record['status']=='passed','Incomplete reference audit')
        check_hash(extension/'independent/check_references.py',record['source_sha256']);check_hash(extension/'independent/core.py',record['core_sha256'])
        require(record['input_design_sha256']==expected[case]['design_sha256'] and record['input_states_sha256']==expected[case]['states_sha256'],'Reference audit inputs changed')
        require((extension/'independent/references'/(case+'.npz')).is_file(),'Missing cached reference')
    fixed=read_json(extension/'fixed/execution_provenance.json')
    check_bindings(qhat,{name:item['sha256'] for name,item in fixed['qhat_sources'].items()})
    execution=read_json(extension/'optimized/EXECUTION_MANIFEST.json')
    check_bindings(qhat,execution['qhat_source_sha256'])
    return summary,{'status':'passed','rows':len(rows),'completed_evolutions':completed,'scope':'Current sources, frozen designs, audited rows/states, analysis inputs, and original QHAT input bytes match previously recorded hashes; no numerical rerun. Historical reference-cache NPZ files lack a recorded output digest and are only checked for presence; published reproduction requires a fresh output directory.'}

def check_output_location(extension,qhat,output):
    output=output.resolve()
    require(not output.is_relative_to(extension.resolve()) and not output.is_relative_to(qhat.resolve()),'Archive output must be outside both input trees')

def main():
    p=argparse.ArgumentParser();p.add_argument('--extension',type=Path,required=True);p.add_argument('--qhat',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.extension=a.extension.resolve();a.qhat=a.qhat.resolve();a.output=a.output.resolve();check_output_location(a.extension,a.qhat,a.output)
    summary,binding_check=preflight(a.extension,a.qhat)
    files={}
    ignored={'__pycache__','cache','numba-cache','mpl-cache','matplotlib','oracle_numba','logs_cache','preview','manuscript-preview','manuscript-build','replay-cache'}
    for path in a.extension.rglob('*'):
        require(not path.is_symlink(),'Symlink archive input: '+str(path))
        if not path.is_file():continue
        rel=path.relative_to(a.extension)
        if any(part in ignored or part.startswith('numba_') or part.startswith('cache_') for part in rel.parts):continue
        if path.name=='initial_replay_check.json':continue
        if path.suffix in ('.pyc','.nbc','.nbi') or '.partial' in path.name or path.name=='v2_exact':continue
        files['extension/'+rel.as_posix()]=path
    tracked=subprocess.check_output(['git','ls-files','-z'],cwd=a.qhat).decode().split('\0')
    for name in tracked:
        if (name.endswith('.py') and (name.split('/')[0] in ('analysis','common','hamiltonian_generator') or name in ('numerical_audit/__init__.py','numerical_audit/run_frozen_audit.py','numerical_audit/reference_tools.py','numerical_audit/followup_classes/__init__.py','numerical_audit/followup_classes/exact_blocks.py'))) or name in ('pyproject.toml','LICENSE.md'):
            files['qhat/'+name]=a.qhat/name
    paths=['numerical_audit/ordering_mechanism/environment-requirements.txt','numerical_audit/ordering_mechanism/heldout/heldout_validation.zip']
    for case in CASES[:8]:
        paths.extend(['numerical_audit/followup_classes/results/comparison/'+case+s for s in ('.design.json','.states.npz')])
        paths.append('numerical_audit/followup_classes/results/bch_supplement/'+case+'.bch_design.json')
        paths.append('numerical_audit/ordering_mechanism/diagonal_refinement/results/'+case+'.design.json')
        for arm in ('class_pauli_JD','class_refined_J0'):paths.append('numerical_audit/ordering_mechanism/objective_ablation/results/'+arm+'/'+case+'.design.json')
    for name in paths:files['qhat/'+name]=a.qhat/name
    require(all(p.is_file() and not p.is_symlink() for p in files.values()),'Missing or symlinked archive input')
    manifest={'scope':'expanded local ordering benchmark; earlier immutable artifacts remain separate','preflight':binding_check,'qhat_base_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.qhat,text=True).strip(),'files':{name:{'bytes':path.stat().st_size,'sha256':file_hash(path)} for name,path in sorted(files.items())}}
    a.output.parent.mkdir(parents=True,exist_ok=True);tmp=a.output.with_suffix('.partial.zip')
    require(tmp!=a.output,'Output cannot use the reserved .partial.zip suffix')
    with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for name,path in sorted(files.items()):z.write(path,name)
        z.writestr('artifact_manifest.json',json.dumps(manifest,indent=2)+'\n')
    verify(tmp) # Detect member changes during packaging before publishing the archive.
    tmp.replace(a.output)
    record={'archive':a.output.name,'sha256':file_hash(a.output),'bytes':a.output.stat().st_size,'files':len(files),'uncompressed_bytes':sum(e['bytes'] for e in manifest['files'].values()),'manifest_sha256':hashlib.sha256((json.dumps(manifest,indent=2)+'\n').encode()).hexdigest(),'benchmark_rows':summary['rows'],'completed_evolutions':summary['completed_evolutions'],'preflight':binding_check}
    a.output.with_suffix('.record.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
if __name__=='__main__':main()
