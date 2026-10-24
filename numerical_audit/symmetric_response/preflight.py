"""Largest altered-target end-to-end preflight before prospective errors."""
import argparse,hashlib,importlib.metadata,json,subprocess,tarfile,time
from pathlib import Path
from study import prepare,run,write
from verify import verify
from reporting import generate
from delivery import pack,check,source_hashes

def preflight(output):
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True);identity=source_hashes();root=Path(__file__).resolve().parents[2]
    commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    for name,digest in identity.items():
        blob=subprocess.check_output(['git','-C',str(root),'show',commit+':numerical_audit/'+name])
        if hashlib.sha256(blob).hexdigest()!=digest:raise ValueError('Commit complete preflight source before execution')
    cpu=time.process_time();wall=time.monotonic();phases={}
    for name,call in [('prepare',lambda:prepare(output/'results',True)),('propagation',lambda:run(output/'results')),
        ('independent',lambda:verify(output/'results')),('report',lambda:generate(output/'results',output/'report')),
        ('package',lambda:pack(output/'results',output/'report',output/'delivery'))]:
        c=time.process_time();w=time.monotonic();call();phases[name]=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.monotonic()-w)
    c=time.process_time();w=time.monotonic();archive=output/'delivery.tar.gz'
    with tarfile.open(archive,'x:gz') as t:t.add(output/'delivery',arcname='delivery')
    extracted=output/'extracted';extracted.mkdir()
    with tarfile.open(archive) as t:t.extractall(extracted,filter='data')
    check(extracted/'delivery');phases['archive_extract_check']=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.monotonic()-w)
    if identity!=source_hashes():raise ValueError('Preflight source changed')
    write(output/'PREFLIGHT.json',dict(status='passed',source_commit=commit,sources=identity,
        environment={p:importlib.metadata.version(p) for p in ('numpy','scipy')},case='H2O_scale1p10',coefficient_scale='7/5',
        modes=8,pauli_words=104,T=2.11,circuits=16,phases=phases,
        algorithm_CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        retained_bytes=sum(p.stat().st_size for p in (output/'delivery').rglob('*') if p.is_file()),archive_bytes=archive.stat().st_size,
        total_preflight_bytes=sum(p.stat().st_size for p in output.rglob('*') if p.is_file()),
        scope='All8arms/bothcaps,newobjective/fullindependentverification/report/archive/extract/check; altered target only'))
    print('PREFLIGHT PASSED',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    preflight(p.parse_args().output.resolve())
