#!/usr/bin/env python3
"""Streaming SHA-256/size inventory check; all gates remain active under -O.

This checks recorded bytes, not mathematical or numerical correctness. An
expected archive hash must come from a separately retained record to detect
joint alteration of a member and its internal manifest.
"""
import argparse,hashlib,json,stat,zipfile
from pathlib import Path,PurePosixPath

def require(ok,message):
    if not ok:raise ValueError(message)
def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
def safe_name(name):
    p=PurePosixPath(name);require(not p.is_absolute() and '..' not in p.parts and '\\' not in name and str(p)==name,'Unsafe member '+name)
def verify(path,expected=None):
    digest=file_hash(path)
    if expected is not None:require(digest==expected,'Archive SHA-256 mismatch')
    with zipfile.ZipFile(path) as z:
        names=z.namelist();require(len(names)==len(set(names)),'Duplicate member')
        for info in z.infolist():
            safe_name(info.filename);require(not stat.S_ISLNK(info.external_attr>>16),'Symlink member')
        manifest=json.loads(z.read('artifact_manifest.json'));entries=manifest['files'];require(set(names)==set(entries)|{'artifact_manifest.json'},'Unlisted or missing members')
        total=0
        for name,e in entries.items():
            safe_name(name);require(z.getinfo(name).file_size==e['bytes'],'Declared size mismatch '+name)
            h=hashlib.sha256();size=0
            with z.open(name) as f:
                for b in iter(lambda:f.read(4*1024**2),b''):h.update(b);size+=len(b)
            require(size==e['bytes'] and h.hexdigest()==e['sha256'],'Member digest/size mismatch '+name);total+=size
    return {'status':'passed','archive_sha256':digest,'members':len(entries),'uncompressed_bytes':total,'scope':'byte integrity; no simulation or proof checking'}
def main():
    p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--expected-sha256');a=p.parse_args();print(json.dumps(verify(a.archive,a.expected_sha256)))
if __name__=='__main__':main()
