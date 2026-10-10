#!/usr/bin/env python3
"""Package a completed control; dependency archive is identified explicitly."""
import resource
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
import json,zipfile
from pathlib import Path
import run as r
ROOT=Path(__file__).resolve().parent
r.require((ROOT/'READY.json').exists(),'Cannot package unfinished verification')
r.require(json.loads((ROOT/'verification.json').read_text())['replays']==504,'Incomplete independent replays')
source_names={p.name for p in ROOT.iterdir() if p.suffix=='.py'}|{'PROTOCOL.md','REPRODUCE.md','CHECKER_CORRECTION.md','environment-requirements.txt'}
files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and not any(x in {'cache','__pycache__'} for x in p.relative_to(ROOT).parts) and p.suffix!='.zip' and p.name not in {'artifact_record.json','artifact_manifest.json'} and not p.name.endswith('.partial')]
entries={str(p.relative_to(ROOT)):{'sha256':r.sha(p),'bytes':p.stat().st_size} for p in files}
manifest={'files':entries,'rows':504,'new_evolutions':420,'reused_baselines':84,'independent_replays':504,'dependency_archive_sha256':'d7f2ae431585e3ed93bb612c06d4558fdc162b196f0111f78b1ae7cc983f8f43','dependency':'Original extension_validation.zip, containing the immutable input and QHAT source snapshot; required for fresh production runs, not for independent verification of the saved control.'}
r.write(ROOT/'artifact_manifest.json',manifest)
records={}
for name,chosen in [('strang_control_validation.zip',files),('strang_control_source.zip',[ROOT/n for n in sorted(source_names)])]:
 out=ROOT/name;mini={'files':{str(p.relative_to(ROOT)):entries[str(p.relative_to(ROOT))] for p in chosen},'dependency_archive_sha256':manifest['dependency_archive_sha256'],'scope':'Complete validation artifact' if 'validation' in name else 'Source-only portable runner and independent checker; run in a fresh directory'}
 with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in chosen:z.write(p,str(p.relative_to(ROOT)))
  z.writestr('artifact_manifest.json',json.dumps(mini,indent=2)+'\n')
 records[name]={'sha256':r.sha(out),'bytes':out.stat().st_size,'files':len(chosen)}
r.write(ROOT/'artifact_record.json',records)
print(json.dumps(records,indent=2))
