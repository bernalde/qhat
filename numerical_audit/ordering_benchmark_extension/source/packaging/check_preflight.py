#!/usr/bin/env python3
"""Tiny byte-binding and archive-location probes; no molecular propagation."""
import argparse,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from build_archive import campaign_bindings,check_output_location,file_hash
    records=[]
    with tempfile.TemporaryDirectory(prefix='ordering-package-preflight-') as tmp:
        base=Path(tmp);ext=base/'extension';root=ext/'fixed';audit=ext/'independent/fixed';case='toy';source=ext/'producer.py';core=ext/'independent/core.py';checker=ext/'independent/check_campaign.py'
        def write(path,data):
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data))
        for f in (source,core,checker):f.parent.mkdir(parents=True,exist_ok=True);f.write_text('# frozen source\n')
        dp=root/'designs/toy.json';write(dp,{'extension_input_provenance':{'design_sha256':'original'}})
        fp=root/'freeze.json';write(fp,{'cases':[case],'design_sha256':{case:file_hash(dp)},'source_sha256':{'producer.py':file_hash(source)}})
        sp=root/'results/toy.method.b1.npz';sp.parent.mkdir(parents=True,exist_ok=True);sp.write_bytes(b'tiny saved state')
        rp=root/'results/toy.method.b1.json';row={'case':case,'ordering':'method','budget_units':1,'extension_design_sha256':file_hash(dp),'status':'complete','state_file':sp.name,'state_sha256':file_hash(sp)};write(rp,row)
        record={'case':case,'status':'passed','design_sha256':file_hash(dp),'source_sha256':file_hash(checker),'core_sha256':file_hash(core),'rows':1,'checks':[{'ordering':'method','budget_units':1,'status':'passed','row_sha256':file_hash(rp),'state_sha256':file_hash(sp)}]};ap=audit/'toy.json';write(ap,record)
        summary={'status':'passed','case_count':1,'rows':1,'cases':[case],'freeze_sha256':file_hash(fp)};su=audit/'summary.json';write(su,summary)
        originals={f:f.read_bytes() for f in (source,core,checker,dp,fp,sp,rp,ap,su)}
        def restore():
            for f,raw in originals.items():f.write_bytes(raw)
        def audit_once():return campaign_bindings(ext,'fixed','fixed',[case],1)
        def test(name,mutate,accept=False):
            restore();mutate()
            try:audit_once();passed=True;message='accepted'
            except (ValueError,KeyError,FileNotFoundError) as e:passed=False;message=str(e)
            if passed!=accept:raise RuntimeError('Unexpected preflight result '+name+': '+message)
            records.append({'probe':name,'accepted':passed,'message':message})
        test('valid',lambda:None,True)
        test('changed_state',lambda:sp.write_bytes(b'changed saved state'))
        test('missing_state',lambda:sp.unlink())
        test('changed_row',lambda:rp.write_text(rp.read_text()+' '))
        test('changed_design',lambda:dp.write_text(dp.read_text()+' '))
        test('changed_producer_source',lambda:source.write_text('# changed\n'))
        test('changed_verifier_source',lambda:checker.write_text('# changed\n'))
        test('changed_freeze',lambda:fp.write_text(fp.read_text()+' '))
        test('stale_summary_count',lambda:write(su,{**summary,'rows':2}))
        test('missing_audited_row',lambda:rp.unlink())
        extra=root/'results/toy.extra.b1.json'
        test('extra_row',lambda:write(extra,{**row,'ordering':'extra'}));extra.unlink()
        test('duplicate_audit_check',lambda:write(ap,{**record,'checks':record['checks']*2}))
        restore()
        for label,path,accept in [('outside',base/'artifact.zip',True),('inside_extension',ext/'artifact.zip',False),('inside_qhat',base/'qhat/artifact.zip',False)]:
            try:check_output_location(ext,base/'qhat',path);passed=True
            except ValueError:passed=False
            if passed!=accept:raise RuntimeError('Bad location guard '+label)
            records.append({'probe':'output_'+label,'accepted':passed})
    result={'status':'passed','python':sys.version,'assertions_disabled':not __debug__,'build_archive_sha256':file_hash(Path(__file__).with_name('build_archive.py')),'probes':records,'scope':'Tiny fixtures for recorded hash bindings and input-tree self-inclusion; no numerical recomputation.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','probes':len(records),'assertions_disabled':not __debug__}))
if __name__=='__main__':main()
