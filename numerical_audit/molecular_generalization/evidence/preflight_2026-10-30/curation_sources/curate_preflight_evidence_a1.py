"""Bounded exact-file numerical evidence curation; no scientific computation."""
import hashlib,json,subprocess,sys,tarfile
from pathlib import Path
planpath=Path(sys.argv[1]).resolve();plan=json.loads(planpath.read_text());dest=Path(plan['destination']);assert not dest.exists();dest.mkdir(parents=True)
repo=Path('/tmp/qhat-molecular-generalization-2026-10-30/worktree');records=[];total=0;archives=[]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def put(relative,raw,record):
    global total
    target=dest/relative;assert not Path(relative).is_absolute() and '..'not in Path(relative).parts
    total+=len(raw);assert len(records)<plan['maximum_files'] and total<=plan['maximum_bytes']
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('xb')as stream:stream.write(raw)
    assert sha(target.read_bytes())==sha(raw)
    records.append(dict(path=relative,sha256=sha(raw),bytes=len(raw),**record))
for item in plan['items']:
    source=Path(item['source']);assert source.is_file() and not source.is_symlink()
    raw=source.read_bytes();assert len(raw)==item['bytes'] and sha(raw)==item['sha256']
    put(item['target'],raw,dict(original_locator=item['source'],group=item['group'],role=item['role']))
    if source.name.endswith('.tar.gz'):
        inventory=None;seen=set()
        with tarfile.open(source,mode='r|gz')as archive:
            for member in archive:
                assert member.isfile() and member.name.startswith('data/') and '..'not in Path(member.name).parts
                parts=Path(member.name).parts;assert not any(p in ('.git','.aws','.codex','coordination')for p in parts)
                assert Path(member.name).name not in ('TASK.json','CURRENT_STATE.md','RESUME.md')
                relative=member.name[len('data/'):];stream=archive.extractfile(member)
                if relative=='EVIDENCE_INVENTORY.json':
                    assert inventory is None and member.size<=8*1024**2
                    inventory=json.loads(stream.read());continue
                assert inventory is not None and relative not in seen
                expected=inventory['files'][relative];assert member.size==expected['bytes']
                digest=hashlib.sha256();count=0
                for chunk in iter(lambda:stream.read(1024**2),b''):digest.update(chunk);count+=len(chunk)
                assert count==member.size and digest.hexdigest()==expected['sha256'];seen.add(relative)
        assert inventory is not None and seen==set(inventory['files'])
        archives.append(dict(path=item['target'],files=len(seen),every_member_authenticated=True,no_extraction=True))
for entry in plan['source_exports']:
    commit,path=entry['commit'],entry['path'];assert len(commit)==40 and all(c in '0123456789abcdef'for c in commit)
    assert path.startswith('numerical_audit/molecular_generalization/') and '..'not in Path(path).parts
    raw=subprocess.check_output(['git','-C',str(repo),'show',commit+':'+path])
    put('source_exports/'+commit+'/'+path,raw,dict(group='committed_source_export',producer_commit=commit,original_git_path=path,role='Exact committed bytes; not a retrospective claim of contemporaneous test-source freezing'))
manifest=dict(schema='molecular-preflight-evidence-v1',files=records,files_count=len(records),uncompressed_bytes=total,archive_authentication=archives,calibration_reference=plan['calibration_reference'],source_lineage=plan['source_lineage'],known_gaps=plan['known_gaps'],curation_plan_sha256=sha(planpath.read_bytes()),scope=plan['publication_scope'])
with(dest/'MANIFEST.json').open('x')as stream:json.dump(manifest,stream,indent=2,allow_nan=False)
print(json.dumps(dict(complete=True,files=len(records),bytes=total,archive_member_checks=archives)),flush=True)
