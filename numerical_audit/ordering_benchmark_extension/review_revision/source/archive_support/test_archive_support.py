#!/usr/bin/env python3
"""Small adversarial checks, always invoking Python -O, without propagation."""
import argparse,hashlib,json,os,resource,stat,subprocess,sys,tempfile,warnings,zipfile
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))

def digest(b):return hashlib.sha256(b).hexdigest()
def require(v,m):
    if not v:raise RuntimeError(m)
def main():
    p=argparse.ArgumentParser();p.add_argument('--support',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(not a.output.exists(),'Choose fresh output')
    support=a.support.resolve();checks=[]
    with tempfile.TemporaryDirectory(prefix='archive-support-test-') as td:
        root=Path(td)
        def archive(name,actual=None,entries=None,extra=(),symlink=None,duplicate=False):
            actual=actual or {'a.txt':b'alpha','b.txt':b'beta'}
            entries=entries or {k:{'bytes':len(v),'sha256':digest(v)} for k,v in actual.items()}
            mb=json.dumps({'files':entries},sort_keys=True).encode();zpath=root/name
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',UserWarning)
                with zipfile.ZipFile(zpath,'w') as z:
                    for key,val in actual.items():
                        if key==symlink:
                            info=zipfile.ZipInfo(key);info.create_system=3;info.external_attr=(stat.S_IFLNK|0o777)<<16;z.writestr(info,val)
                        else:z.writestr(key,val)
                    for key,val in extra:z.writestr(key,val)
                    if duplicate:z.writestr('a.txt',b'alpha')
                    z.writestr('artifact_manifest.json',mb)
            return zpath,digest(zpath.read_bytes()),digest(mb)
        valid,sha,msha=archive('valid.zip')
        def check(name,args,passed,env=None):
            done=subprocess.run([sys.executable,'-O',*map(str,args)],cwd=root,env=env,capture_output=True,text=True,timeout=30)
            require((done.returncode==0)==passed,name+': unexpected exit '+str(done.returncode)+' '+done.stderr)
            checks.append({'name':name,'expected_pass':passed,'exit_code':done.returncode,'stdout':done.stdout.strip(),'stderr_tail':done.stderr.strip().splitlines()[-1:]})
        verifier=support/'verify_archive.py'
        check('valid-external-identity',[verifier,valid,'--expected-sha256',sha,'--expected-manifest-sha256',msha,'--expected-members','2'],True)
        check('explicit-integrity-only',[verifier,valid,'--integrity-only'],True)
        check('reject-missing-identity-mode',[verifier,valid],False)
        check('reject-wrong-archive-digest',[verifier,valid,'--expected-sha256','0'*64],False)
        check('reject-wrong-manifest-digest',[verifier,valid,'--expected-sha256',sha,'--expected-manifest-sha256','0'*64],False)
        check('reject-wrong-member-count',[verifier,valid,'--expected-sha256',sha,'--expected-members','3'],False)
        stale={'a.txt':{'bytes':5,'sha256':digest(b'alpha')},'b.txt':{'bytes':4,'sha256':digest(b'beta')}}
        changed,_,_=archive('changed.zip',actual={'a.txt':b'ALPHA','b.txt':b'beta'},entries=stale)
        check('reject-changed-member-under-O',[verifier,changed,'--integrity-only'],False)
        changed_both,_,_=archive('changed-both.zip',actual={'a.txt':b'ALPHA','b.txt':b'beta'})
        check('reject-joint-member-manifest-change',[verifier,changed_both,'--expected-sha256',sha],False)
        missing,_,_=archive('missing.zip',actual={'a.txt':b'alpha'},entries=stale)
        check('reject-missing-member',[verifier,missing,'--integrity-only'],False)
        unlisted,_,_=archive('unlisted.zip',extra=[('c.txt',b'gamma')])
        check('reject-unlisted-member',[verifier,unlisted,'--integrity-only'],False)
        duplicate,_,_=archive('duplicate.zip',duplicate=True)
        check('reject-duplicate-member',[verifier,duplicate,'--integrity-only'],False)
        unsafe,_,_=archive('unsafe.zip',actual={'../a.txt':b'alpha'})
        check('reject-unsafe-path',[verifier,unsafe,'--integrity-only'],False)
        link,_,_=archive('symlink.zip',symlink='a.txt')
        check('reject-symlink',[verifier,link,'--integrity-only'],False)
        arc=root/'extracted';(arc/'qhat').mkdir(parents=True);(arc/'extension').mkdir()
        guard=support/'replay_guarded.py';env=dict(os.environ,PYTHONPATH='',PYTHONDONTWRITEBYTECODE='1')
        check('valid-single-namespace',[guard,'--archive-root',arc,'--check-imports-only','--output',root/'valid_guard.json'],True,env)
        live=root/'live';(live/'qhat').mkdir(parents=True)
        check('reject-merged-namespace',[guard,'--archive-root',arc,'--check-imports-only','--output',root/'bad_guard.json'],False,dict(env,PYTHONPATH=str(live)))
        (arc/'extension/packaging').mkdir()
        stub=arc/'extension/packaging/replay_one.py'
        stub.write_text("from pathlib import Path\nimport sys,json\ndef main():\n p=Path(sys.argv[sys.argv.index('--output')+1]);p.write_text(json.dumps({'status':'passed','scope':'mock runner only'}))\n")
        check('guarded-delegation-mock',[guard,'--archive-root',arc,'--campaign','fixed','--case','mock','--ordering','mock','--budget','2','--output',root/'mock_replay.json'],True,env)
        stub.write_text(stub.read_text()+" import types\n x=types.ModuleType('qhat.late_escape');x.__file__='/tmp/outside-qhat.py';sys.modules[x.__name__]=x\n")
        check('reject-late-import-before-final-output',[guard,'--archive-root',arc,'--campaign','fixed','--case','mock','--ordering','mock','--budget','2','--output',root/'late_escape.json'],False,env)
        require(not (root/'late_escape.json').exists(),'Failed guard published a final success record')
        module_code="import importlib.util,sys,types;from pathlib import Path;s=importlib.util.spec_from_file_location('guard',sys.argv[1]);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);r=Path(sys.argv[2]);sys.path.insert(0,str(r.parent));import qhat;x=types.ModuleType('qhat.escaped');x.__file__=sys.argv[3];sys.modules[x.__name__]=x;m.check_qhat_imports(r)"
        check('reject-module-outside-archive',['-c',module_code,guard,arc/'qhat',live/'outside.py'],False,env)
    result={'status':'passed','checks':checks,'count':len(checks),'python':sys.version,'scope':'adversarial -O verifier and namespace tests; no molecular propagation'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','count':len(checks)}))
if __name__=='__main__':main()
