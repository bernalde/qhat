#!/usr/bin/env python3
"""Mandatory provenance-checked entry point for interrupted schedule recovery."""
import argparse,json,resource,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import common as c
import run

def prepare_checked(case,adopt=False):
    raw,_=c.load_case(case)
    provenance={'case':case,'input_provenance':raw['extension_input_provenance'],'scientific_source_sha256':run.source_hashes(),'recovery_wrapper_sha256':c.sha(__file__),'algorithm':'two starts, up to ten accepted best insertions; JD/J0 and Pauli budget definitions fixed by PROTOCOL.md'}
    sidecar=HERE/'checkpoints'/f'{case}.provenance.json';existing=list((HERE/'checkpoints').glob(f'{case}.*.start*.json'))
    if sidecar.exists():
        saved=json.loads(sidecar.read_text());c.require(saved['binding']==provenance,'Checkpoint provenance differs from current inputs or source')
    else:
        if existing:
            c.require(adopt,'Unbound checkpoints require an explicit one-time validated adoption')
            design=HERE/'designs'/f'{case}.design.json';c.require(design.exists(),'Cannot adopt unfinished checkpoints without a completed design')
            d=json.loads(design.read_text());c.require(d['source_sha256']==run.source_hashes() and d['extension_input_provenance']==raw['extension_input_provenance'],'Cannot adopt changed source/input')
            for p in existing:
                s=json.loads(p.read_text());c.require(s['complete'],'Cannot adopt an unfinished preguard checkpoint')
            note='Created earlier in this uninterrupted recovery session before the provenance wrapper was added; all searches completed under the unchanged source hashes recorded in the completed design. Development methods additionally reproduced the three archived orders exactly.'
        else:note='Provenance binding created before any checkpoint or search.'
        c.write_atomic(sidecar,{'binding':provenance,'binding_sha256':c.object_sha(provenance),'creation_note':note})
    run.prepare(case)
    # Store completion hashes separately from the immutable input/source binding.
    checkpoints={p.name:c.sha(p) for p in sorted((HERE/'checkpoints').glob(f'{case}.*.start*.json'))}
    c.write_atomic(HERE/'checkpoints'/f'{case}.completed_hashes.json',{'binding_sha256':c.object_sha(provenance),'checkpoints':checkpoints,'design_sha256':c.sha(HERE/'designs'/f'{case}.design.json')})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cases',nargs='*',default=c.CASES);ap.add_argument('--adopt-current-checkpoints',action='store_true');a=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    for case in a.cases:prepare_checked(case,a.adopt_current_checkpoints)
if __name__=='__main__':main()
