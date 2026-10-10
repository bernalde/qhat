#!/usr/bin/env python3
"""Archive this review's additions without rewriting the prior experiment."""
import argparse,hashlib,importlib.util,json,shutil,zipfile
from pathlib import Path
BASE_SHA='d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--revision',type=Path,required=True);p.add_argument('--qhat',type=Path,required=True);a=p.parse_args();dest=a.qhat/'numerical_audit/ordering_benchmark_extension/review_revision';dest.mkdir(parents=True,exist_ok=True)
 v=json.loads((a.revision/'strang_controls/verification.json').read_text())
 if v.get('status') not in ('passed','pass'):raise ValueError('Independent replay incomplete')
 if not (a.revision/'strang_controls/record_audit.json').is_file():raise ValueError('Record chains unaudited')
 files={}
 for sub in ('strang_controls','archive_support','analysis'):
  for src in sorted((a.revision/sub).rglob('*')):
   rel=src.relative_to(a.revision)
   if not src.is_file() or any(x in ('__pycache__','cache','mpl-cache','preview') for x in rel.parts) or src.suffix in ('.pyc','.nbc','.nbi') or '.partial' in src.name:continue
   files[rel.as_posix()]=src
 # The archive carries full inputs and states. The checkout keeps the readable
 # sources/summaries plus the same JSON inventory; large arrays stay in the ZIP.
 source={}
 for name,src in files.items():
  if src.suffix=='.npz' or '/verification/' in name or '/rows/' in name or '/designs/' in name:continue
  target=dest/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,target);source[str(target.relative_to(dest))]={'bytes':target.stat().st_size,'sha256':sha(target)}
 archive=dest/'review_supplement.zip';manifest={'scope':'Supplement to the unchanged original extension archive, which is required for archived QHAT source and legacy inputs. New control states, portable reproduction, missing historical arrays and recomputed review analysis.','requires_original_archive_sha256':BASE_SHA,'files':{name:{'bytes':src.stat().st_size,'sha256':sha(src)} for name,src in files.items()}}
 manifest_bytes=(json.dumps(manifest,indent=2)+'\n').encode()
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name,src in files.items():
   info=zipfile.ZipInfo(name,date_time=(2026,10,10,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
   with src.open('rb') as inp,z.open(info,'w') as out:shutil.copyfileobj(inp,out,1024**2)
  info=zipfile.ZipInfo('artifact_manifest.json',date_time=(2026,10,10,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,manifest_bytes)
 record=dict(archive_file=archive.name,archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),members=len(files),requires_original_archive_sha256=BASE_SHA,source_sha256=sha(Path(__file__)))
 (dest/'archive_record.json').write_text(json.dumps(record,indent=2)+'\n');(dest/'artifact_manifest.json').write_bytes(manifest_bytes);(dest/'source_manifest.json').write_text(json.dumps(source,indent=2)+'\n')
 verifier=a.qhat/'numerical_audit/ordering_benchmark_extension/verify_archive.py';spec=importlib.util.spec_from_file_location('archive_verifier',verifier);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);r=m.verify(archive,record['archive_sha256'],expected_manifest_sha256=record['manifest_sha256'],expected_members=record['members']);(dest/'archive_check.json').write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(record,indent=2))
if __name__=='__main__':main()
