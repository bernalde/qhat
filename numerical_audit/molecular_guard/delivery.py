"""Package producer bytes and independently verify the extracted numerical data.

Only the live, committed checker is executed. The archive's Python is retained
as evidence and is never imported. One full extracted audit also checks the
original delivery because every payload byte is compared before that audit.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import subprocess
import zipfile
from independent import verify

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MAX_FILES,MAX_BYTES,MAX_FILE_BYTES = 512,64*1024**2,16*1024**2


def check(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def package(producer, output):
    producer,output = producer.resolve(),output.resolve()
    check(not output.exists() and not output.is_relative_to(producer),'Use a fresh separate delivery directory')
    files = {}
    for path in sorted(producer.rglob('*')):
        check(not path.is_symlink(),'No symlinks in certificate delivery')
        if not path.is_file():
            continue
        name = path.relative_to(producer).as_posix()
        check(not name.endswith('.partial'),'Retain/fix interrupted producer outside this delivery')
        size = path.stat().st_size
        check(size <= MAX_FILE_BYTES and len(files) < MAX_FILES,'Per-file/count delivery guard')
        files[name] = dict(bytes=size,sha256=digest(path))
    check(files and sum(row['bytes'] for row in files.values()) <= MAX_BYTES,'Total delivery guard')
    output.mkdir(parents=True)
    manifest = output/'DELIVERY_MANIFEST.json'
    manifest.write_text(json.dumps(files,indent=2)+'\n')
    archive = output/'certificate.zip'
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,path in [('DELIVERY_MANIFEST.json',manifest)]+[('results/'+name,producer/name) for name in files]:
            info = zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (stat.S_IFREG | 0o600) << 16
            z.writestr(info,path.read_bytes(),compresslevel=6)
    extracted = output/'extracted'
    with zipfile.ZipFile(archive) as z:
        seen,total = set(),0
        for info in z.infolist():
            name = PurePosixPath(info.filename)
            check(not name.is_absolute() and '..' not in name.parts and str(name) == info.filename,
                  'Unsafe archive path')
            check(info.filename not in seen and stat.S_IFMT(info.external_attr >> 16) == stat.S_IFREG,
                  'Duplicate or nonregular archive member')
            check(info.file_size <= MAX_FILE_BYTES,'Extracted member bound')
            seen.add(info.filename)
            total += info.file_size
        check(len(seen) == len(files)+1 and total <= MAX_BYTES+MAX_FILE_BYTES,'Extraction count/size guard')
        check(seen == {'DELIVERY_MANIFEST.json'}|{'results/'+name for name in files},'Archive inventory differs')
        z.extractall(extracted)
    check((extracted/'DELIVERY_MANIFEST.json').read_bytes() == manifest.read_bytes(),'Extracted manifest changed')
    actual = {p.relative_to(extracted/'results').as_posix() for p in (extracted/'results').rglob('*') if p.is_file()}
    check(actual == set(files),'Extracted payload inventory differs')
    for name,record in files.items():
        check(digest(producer/name) == digest(extracted/'results'/name) == record['sha256'] and
              (extracted/'results'/name).stat().st_size == record['bytes'],'Delivered payload bytes changed')
    return files,archive,extracted/'results'


def deliver(producer, output):
    head = subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    checks = {}
    for name in ('delivery.py','independent.py'):
        path = HERE/name
        blob = subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+str(path.relative_to(ROOT))])
        check(blob == path.read_bytes(),'Commit trusted delivery/checker before audit')
        checks[str(path.relative_to(ROOT))] = digest(path)
    files,archive,extracted = package(producer,output)
    (output/'TRUSTED_VERIFIER_SOURCES.json').write_text(json.dumps(dict(head=head,sha256=checks),indent=2)+'\n')
    verify(extracted,output/'INDEPENDENT_VERIFIED.json')
    audit = json.loads((output/'INDEPENDENT_VERIFIED.json').read_text())
    record = dict(status='passed',payload_files=len(files),payload_bytes=sum(x['bytes'] for x in files.values()),
        archive_bytes=archive.stat().st_size,archive_sha256=digest(archive),
        producer_freeze_sha256=audit['producer_freeze_sha256'],trusted_checker_head=head,
        circuits=audit['circuits'],reference_evolutions=audit['reference_evolutions'],
        replay='One complete independent audit on extracted data;all original/extracted bytes are identical',
        archived_python_executed=False,scope=audit['scope'])
    (output/'DELIVERY_VERIFIED.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    deliver(args.producer,args.output)
