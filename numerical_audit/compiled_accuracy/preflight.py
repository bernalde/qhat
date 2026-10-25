"""Largest altered-model complete compiled-grid preflight before main errors."""
import argparse,hashlib,importlib.metadata,json,subprocess,tarfile,time
from pathlib import Path
import numpy as np
from compiler import require
from study import prepare,run,write,sources,REPO
from verify import verify
from reporting import report
from delivery import pack,check

def preflight(output,ordering_evidence,ordering_commit):
    require(not output.exists(),'Fresh preflight output required');output.mkdir(parents=True)
    original=sources();commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()
    for name,digest in original.items():
        blob=subprocess.check_output(['git','-C',str(REPO),'show',commit+':'+name])
        require(hashlib.sha256(blob).hexdigest()==digest,'Commit complete preflight pipeline before execution')
    cpu=time.process_time();wall=time.monotonic();phases={}
    for name,call in [
        ('prepare',lambda:prepare(output/'results',ordering_evidence,ordering_commit,preflight=True)),
        ('compiled_propagation',lambda:run(output/'results')),
        ('independent',lambda:verify(output/'results')),
        ('report',lambda:report(output/'results',output/'report')),
        ('package',lambda:pack(output/'results',output/'report',output/'delivery'))]:
        c=time.process_time();w=time.monotonic();call();phases[name]=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.monotonic()-w)
    c=time.process_time();w=time.monotonic();archive=output/'delivery.tar.gz'
    with tarfile.open(archive,'x:gz') as stream:stream.add(output/'delivery',arcname='delivery')
    extracted=output/'extracted';extracted.mkdir()
    with tarfile.open(archive) as stream:stream.extractall(extracted,filter='data')
    check(extracted/'delivery');phases['archive_extract_check']=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.monotonic()-w)
    require(original==sources(),'Source changed during preflight')
    verified=json.loads((output/'results/VERIFIED.json').read_text());require(verified['candidates']==64 and verified['curves']==8 and verified['target_decisions']==16,'Incomplete representative compiled grid')
    write(output/'PREFLIGHT.json',dict(status='passed',source_commit=commit,sources=original,
        environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')},longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),case='preflight_altered_H2O_scale1p10',coefficient_scale='7/5',modes=8,pauli_words=104,T=2.11,
        curves=8,candidates=64,feasible=verified['feasible'],censored=verified['censored'],target_decisions=16,phases=phases,
        algorithm_CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        retained_delivery_bytes=sum(p.stat().st_size for p in (output/'delivery').rglob('*') if p.is_file()),archive_bytes=archive.stat().st_size,
        total_preflight_bytes=sum(p.stat().st_size for p in output.rglob('*') if p.is_file()),
        scope='All8fixed30Utrainingarms/all8stepgridpoints, exactcompleteorpartialcircuits, fullindependentreplay/DAG/targetdecisions,report,archive/extract/check; altered target only'))
    print('COMPILED PREFLIGHT PASSED',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ordering-evidence',type=Path,required=True);parser.add_argument('--ordering-commit',required=True)
    a=parser.parse_args();preflight(a.output.resolve(),a.ordering_evidence.resolve(),a.ordering_commit)
