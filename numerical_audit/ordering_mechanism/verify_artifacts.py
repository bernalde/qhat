#!/usr/bin/env python3
"""Verify saved research bytes, including the embedded held-out archive.

Integrity checks only: this does not execute simulations or prove results.
No validation depends on Python assert or optimization mode.
"""
import argparse,hashlib,json
from pathlib import Path,PurePosixPath
from zipfile import ZipFile

def check_name(name):
 p=PurePosixPath(name)
 if p.is_absolute() or '..' in p.parts:raise ValueError(f'Unsafe relative path: {name}')

def check_bytes(data,entry,name):
 if len(data)!=entry['bytes']:raise ValueError(f'Size mismatch: {name}')
 if hashlib.sha256(data).hexdigest()!=entry['sha256']:raise ValueError(f'Hash mismatch: {name}')

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--manifest',type=Path);a=p.parse_args();root=a.root.resolve();m=json.loads((a.manifest or root/'artifact_manifest.json').read_text())
 for name,entry in m['files'].items():
  check_name(name);file=(root/name).resolve()
  if not file.is_relative_to(root):raise ValueError(f'Escaping artifact: {name}')
  check_bytes(file.read_bytes(),entry,name)
 zipped_count=0
 for name in m.get('embedded_archives',[]):
  with ZipFile(root/name) as z:
   names=z.namelist()
   if len(names)!=len(set(names)):raise ValueError('Duplicate ZIP member')
   for n in names:check_name(n)
   inner=json.loads(z.read('artifact_manifest.json'))
   for n,entry in inner['files'].items():check_bytes(z.read(n),entry,name+':'+n)
   if set(names)!=set(inner['files'])|{'artifact_manifest.json'}:raise ValueError('Unlisted ZIP member')
   zipped_count+=len(inner['files'])
 print(json.dumps({'status':'passed','artifact_files':len(m['files']),'embedded_verified_files':zipped_count,'scope':'byte integrity only'}))
if __name__=='__main__':main()
