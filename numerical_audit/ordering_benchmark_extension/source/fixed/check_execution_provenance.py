#!/usr/bin/env python3
"""Record import resolution, unchanged QHAT sources, and worker process identity.

This is a separate observer; it does not claim to introspect another Python
process's sys.modules. Import resolution is probed with the runner's identical
import sequence and executable, and worker commands/PIDs are recorded.
"""
import importlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
from common import *
MODULES=('analysis.benchmark_b2_active_spaces_matrix_free','analysis.benchmark_b2_signed_coefficient_baseline','analysis.benchmark_L_sweep_trotter','analysis.state_metrics','numerical_audit.run_frozen_audit','numerical_audit.reference_tools','numerical_audit.followup_classes.exact_blocks','analysis.fermionic_coefficients')
record={'scope':'Same-executable, same-source import-resolution probe plus procfs worker identity; not remote sys.modules introspection','executable':sys.executable,'qhat_root':str(QHAT),'freeze_sha256':sha(HERE/'freeze.json'),'qhat_sources':{},'workers':[]}
for name in MODULES:
 module=importlib.import_module('qhat.'+name);path=Path(module.__file__).resolve();require(path.is_relative_to(QHAT),'Imported unintended QHAT checkout '+str(path));rel=str(path.relative_to(QHAT));record['qhat_sources'][rel]={'resolved_file':str(path),'sha256':sha(path)}
record['qhat_commit']=subprocess.check_output(['git','-C',str(QHAT),'rev-parse','HEAD'],text=True).strip();record['qhat_status_porcelain']=subprocess.check_output(['git','-C',str(QHAT),'status','--porcelain'],text=True)
for path in sorted(HERE.glob('job_*.json')):
 raw=path.read_bytes();job=json.loads(raw);snapshot=HERE/'provenance_job_snapshots'/path.name;snapshot.parent.mkdir(exist_ok=True)
 if snapshot.exists():require(snapshot.read_bytes()==raw,'Worker start snapshot already differs')
 else:snapshot.write_bytes(raw)
 p=Path('/proc')/str(job['pid']);entry={'job_file_snapshot':str(snapshot.relative_to(HERE)),'job_sha256':sha(snapshot),'pid':job['pid'],'recorded_command':job['command']}
 if p.exists():
  actual=(p/'cmdline').read_bytes().decode().rstrip('\0').split('\0');require(actual==job['command'],'Worker command changed');entry.update(actual_command=actual,actual_executable=os.readlink(p/'exe'),actual_cwd=os.readlink(p/'cwd'),process_start_ticks=(p/'stat').read_text().split()[21])
 else:entry['procfs_status']='worker no longer running'
 record['workers'].append(entry)
write_atomic(HERE/'execution_provenance.json',record);print('PASS',len(record['qhat_sources']),'QHAT modules;',len(record['workers']),'worker identities')
