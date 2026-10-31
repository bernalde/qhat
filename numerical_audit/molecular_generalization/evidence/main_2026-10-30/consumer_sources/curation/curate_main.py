"""UNEXECUTED stdlib-only exact-byte main evidence curation utility draft.

Root adopts, pins this source and inventory, and wraps ONE invocation in the
original guard20CPU/120wall/512MiB/no-swap. No J, native score, propagation,
project imports, framework tests, installs, client, Git mutation or archive
extraction occurs. Streaming member hashing does not materialize archive files.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import tarfile

R = Path('/home/bernalde/repos/Fermionic_Trotterization/dist/molecular_generalization_2026-10-30')
RUN = R/'runs/main-a1'
NEW = Path('/tmp/qhat-molecular-generalization-2026-10-30/worktree')
TARGET = NEW/'numerical_audit/molecular_generalization/evidence/main_2026-10-30'
FREEZE = 'f3ed635ac7831edcb3feb51d96e6116c2c95f7ad9e2b169b006ca862006776e3'
MAX_BYTES = 512*1024**2


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for raw in iter(lambda:stream.read(1024**2),b''):
            digest.update(raw)
    return digest.hexdigest()


def read(path):
    path=Path(path);size=path.stat().st_size
    require(0<=size<=8*1024**2,'Bounded metadata input')
    with path.open('rb') as stream:
        raw=stream.read(size+1)
    require(len(raw)==size,'Stable metadata bytes')
    return json.loads(raw)


def write_new(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def relative(value):
    require(isinstance(value,str) and value and '\\' not in value and
            not Path(value).is_absolute() and str(Path(value))==value and
            all(part not in ('','.','..') for part in value.split('/')),'Canonical relative destination/member')
    return value


def source_ok(path):
    path=Path(path)
    require(path.is_absolute() and path.resolve(strict=True)==path and not path.is_symlink(),
            'Only canonical real source files')
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode),'Regular source file')
    require(path.is_relative_to(RUN) or path.parent==R/'preflight/replay_main_censored_adopted'
            or path==R/'planning/main_curation_draft/SERVICE_RESOURCE_ACCOUNTING.DRAFT.json',
            'Explicit closed numerical evidence source boundary')
    require(not any(part in ('reviews','.aws','.codex','.agents') or 'handoff' in part.lower()
            or part in ('Task.json','task.json','SOURCE_BINDINGS.json') for part in path.parts),
            'No account/review/task/helper material')


def copy_checked(source,destination,expected_bytes,expected_sha):
    source=Path(source);destination=Path(destination)
    descriptor=os.open(source,os.O_RDONLY|os.O_NOFOLLOW)
    destination.parent.mkdir(parents=True,exist_ok=True)
    digest=hashlib.sha256();count=0
    with os.fdopen(descriptor,'rb') as src,destination.open('xb') as dst:
        before=os.fstat(src.fileno())
        require(stat.S_ISREG(before.st_mode) and before.st_size==expected_bytes,'Expected source length/type')
        while True:
            raw=src.read(min(1024**2,expected_bytes-count+1))
            if not raw:break
            count+=len(raw);require(count<=expected_bytes,'Growing source refused')
            digest.update(raw);dst.write(raw)
        require(count==expected_bytes and digest.hexdigest()==expected_sha,'Exact source bytes preserved')
    os.chmod(destination,0o644)
    require(destination.stat().st_size==expected_bytes and sha(destination)==expected_sha,'Destination byte readback')


def archive_member_proof(archive,index_path,archive_sha,index_sha):
    require(sha(archive)==archive_sha and sha(index_path)==index_sha,'Actual archive/inventory bindings')
    index=read(index_path)
    expected={'data/'+relative(name):value for name,value in index['files'].items()}
    expected['data/EVIDENCE_INVENTORY.json']=dict(bytes=index_path.stat().st_size,sha256=index_sha)
    seen=set();members=[];total=0
    with tarfile.open(archive,'r|gz') as stream:
        for member in stream:
            require(member.isfile() and member.name in expected and member.name not in seen,
                    'Only exact unique declared regular members')
            relative(member.name);value=expected[member.name]
            require(member.size==value['bytes'] and 0<=member.size<=MAX_BYTES,'Exact bounded member size')
            total+=member.size;require(total<=MAX_BYTES+32*1024**2,'Total streaming archive admission')
            digest=hashlib.sha256();count=0
            with stream.extractfile(member) as source:
                for raw in iter(lambda:source.read(1024**2),b''):
                    count+=len(raw);require(count<=member.size,'Growing/truncated archive member')
                    digest.update(raw)
            require(count==member.size and digest.hexdigest()==value['sha256'],'Every member byte matches original sealed inventory')
            seen.add(member.name);members.append(dict(path=member.name,bytes=count,sha256=digest.hexdigest()))
    require(seen==set(expected),'No missing/extra archive member')
    return dict(schema='molecular-main-curation-member-proof-v1',member_bytes_verified=True,
        archive_sha256=archive_sha,inventory_sha256=index_sha,declared_DATA_files=len(index['files']),
        regular_archive_members=len(members),inventory_member_additional=True,total_member_bytes=total,
        members=members,scope='Streaming byte/member proof only; no extraction, J, numerical propagation or scientific gate repair')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',type=Path,required=True)
    parser.add_argument('--inventory-sha256',required=True)
    parser.add_argument('--expected-operator-sha256',required=True)
    parser.add_argument('--target',type=Path,required=True)
    args=parser.parse_args()
    require(sha(__file__)==args.expected_operator_sha256 and sha(args.inventory)==args.inventory_sha256,
            'Actual adopted utility/inventory identity')
    inventory=read(args.inventory);target=args.target.resolve()
    require(target==TARGET and not target.exists() and target.parent.resolve()==target.parent,
            'Fresh exact repository evidence target only')
    require(inventory['target_root']==str(target) and inventory['source_freeze_sha256']==FREEZE,
            'Exact source→target declaration')
    require(sha(RUN/'data/INPUT_PROTOCOL_FROZEN.json')==FREEZE and
            sha(RUN/'data/EVIDENCE_INVENTORY.json')==inventory['sealed_inventory_sha256'],
            'Original sealed freeze/index bytes are unchanged')
    entries=inventory['payload_entries'];names=set();total=0
    require(0<len(entries)<=4096,'Finite explicit whitelist')
    for row in entries:
        source=Path(row['source']);source_ok(source);name=relative(row['target'])
        require(name not in names,'One target per explicit source binding');names.add(name)
        require(Path(name).name not in ('Project.toml','Manifest.toml'),'No loose historical project configuration')
        require(not name.endswith('.toml') or row['role']=='scientific_model_toml',
                'Only scientific model inputTOMLs keep active extension')
        require(source.stat().st_size==row['bytes'] and sha(source)==row['sha256'],'Every source hash/size preflight')
        total+=row['bytes'];require(total<=MAX_BYTES,'Bounded cumulative copy size')
    require(len(entries)==inventory['payload_file_count'] and total==inventory['payload_bytes'],'Full declared whitelist totals')
    require(not any(a!=b and b.startswith(a+'/') for a in names for b in names),'No file/directory destination collision')
    require(not names & {'ARCHIVE_MEMBER_PROOF.json','PAYLOAD_MANIFEST.json','consumer_sources/curation/curate_main.py'},
            'Generated filenames are separate from fixed source payload')
    replay=read(RUN/'relocated-replay-a1/REPLAY_FINAL.json');analysis=read(RUN/'data/ANALYSIS.json')
    require(replay['delivery_complete'] is True and replay['replay_complete'] is True and
            replay['scientific_gate_complete'] is False and replay['raw_claims_verified'] is False and
            analysis['complete'] is False and all(x['decision']=='incomplete' for x in analysis['primary'].values()),
            'Preserve actual delivery/scientific distinction and both incomplete primary decisions')
    target.mkdir()
    actual=[]
    for row in entries:
        copy_checked(row['source'],target/row['target'],row['bytes'],row['sha256']);actual.append(dict(row))
    # This source is now actually executing in the owning curation service.
    own='consumer_sources/curation/curate_main.py';own_path=Path(__file__).resolve()
    copy_checked(own_path,target/own,own_path.stat().st_size,args.expected_operator_sha256)
    actual.append(dict(source=str(own_path),target=own,bytes=own_path.stat().st_size,
        sha256=args.expected_operator_sha256,role='actual_executing_curation_source',byte_policy='Exact executed source'))
    archive=RUN/'relocated-replay-a1/main-a1.data.tar.gz';index=RUN/'data/EVIDENCE_INVENTORY.json'
    proof=archive_member_proof(archive,index,inventory['source_archive_sha256'],inventory['sealed_inventory_sha256'])
    write_new(target/'ARCHIVE_MEMBER_PROOF.json',proof)
    actual.append(dict(source='Generated by this owning curation service',target='ARCHIVE_MEMBER_PROOF.json',
        bytes=(target/'ARCHIVE_MEMBER_PROOF.json').stat().st_size,sha256=sha(target/'ARCHIVE_MEMBER_PROOF.json'),
        role='actual_streaming_member_byte_proof'))
    require({str(q.relative_to(target)) for q in target.rglob('*') if q.is_file()}=={x['target'] for x in actual},
            'No omitted/extra copied or generated payload file')
    write_new(target/'PAYLOAD_MANIFEST.json',dict(schema='molecular-main-curated-payload-v1',
        payload_byte_copy_verified=True,archive_member_bytes_verified=True,
        scientific_gate_complete=False,source_analysis_complete=False,raw_claims_verified=False,
        source_primary_decisions=dict(window_minus_full='incomplete',annealing_minus_full='incomplete'),
        source_commit=inventory['source_commit'],source_freeze_sha256=FREEZE,
        source_inventory=dict(path=str(args.inventory.resolve()),sha256=args.inventory_sha256),
        actual_executing_curation_source=dict(path=str(own_path),sha256=args.expected_operator_sha256),
        copied_existing_files=len(entries),generated_payload_files=2,payload_files=len(actual),
        payload_bytes=sum(x['bytes'] for x in actual),files=actual,
        convenience_subset_is_complete_DATA=False,original_full_DATA_declared_files=proof['declared_DATA_files'],
        historical_embedded_paths_rewritten=False,manifest_excludes_itself=True,
        post_close_root_additions_outside_fixed_payload=['curationPLAN/log/receipt/FINAL','postcurationresourcefields',
            'REPORT','portableREPLAYrecipe','topmanifest'],
        qualification='Exact byte curation and streaming archive/member proof only; source scientific gates and native censor unchanged'))


if __name__=='__main__':
    main()
