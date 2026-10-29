"""Freeze every completed attempt from RUN_ROOT into its contained DATA_ROOT.

--campaign RUN_ROOT contains RESOURCE_POOLS.json, jobs/, logs/, and campaign/
with resource_guard.py, run_guarded_job.py and optional dispatcher_versions/.
--root is DATA_ROOT (normally RUN_ROOT/data). Guard --output remains RUN_ROOT:
disk receipts measure cumulative science, metadata, logs and source storage.
Only the standard library is used; archived code is never executed.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil

MAX_BYTES = 32 * 1024**2
LABEL = re.compile(r'[A-Za-z0-9][A-Za-z0-9-]{0,127}\Z', re.ASCII)
HASH = re.compile(r'[0-9a-f]{64}\Z', re.ASCII)
DISK_SCOPE = ('cumulative RUN_ROOT including DATA_ROOT, jobs, logs, guards and '
              'dispatcher versions; not per-job output bytes')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    require(path.is_file() and path.stat().st_size <= MAX_BYTES, 'Missing/oversized evidence: '+str(path))
    digest = hashlib.sha256()
    total = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            total += len(chunk)
            require(total <= MAX_BYTES, 'Evidence grew beyond byte cap')
            digest.update(chunk)
    return digest.hexdigest()


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, 'Duplicate JSON key: '+key)
        result[key] = value
    return result


def read(path):
    require(path.is_file() and path.stat().st_size <= MAX_BYTES, 'Bounded JSON evidence required')
    with path.open('rb') as stream:
        raw = stream.read(MAX_BYTES+1)
    require(len(raw) <= MAX_BYTES, 'JSON evidence grew beyond byte cap')
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(raw, object_pairs_hook=unique_object, parse_constant=invalid)


def contained(base, spelling):
    require(isinstance(spelling, str) and '\\' not in spelling, 'Relative POSIX evidence path required')
    relative = PurePosixPath(spelling)
    require(not relative.is_absolute() and str(relative) == spelling and relative.parts and
            all(part not in ('.', '..') for part in relative.parts), 'Canonical relative evidence path required')
    path = (base/spelling).resolve(strict=True)
    require(path.is_relative_to(base) and path.is_file(), 'Evidence source escapes RUN_ROOT: '+spelling)
    return path


def destination(data_root, spelling):
    path = data_root/spelling
    require(path.resolve().is_relative_to(data_root), 'Archive destination escapes DATA_ROOT')
    if path.exists() or path.is_symlink():
        require(path.is_file(), 'Existing archive destination is not a file')
    require(not path.with_name(path.name+'.partial').exists(), 'Preserve partial archive output')
    return path


def write_new(path, value):
    temporary = path.with_name(path.name+'.partial')
    require(not path.exists() and not path.is_symlink() and not temporary.exists(), 'Preserve original ledger cutoff')
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    os.link(temporary, path)
    temporary.unlink()


def copy_evidence(source, target, expected):
    if target.exists():
        require(sha(target) == expected, 'Stale archived evidence differs: '+str(target))
        return
    temporary = target.with_name(target.name+'.partial')
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as stream, temporary.open('xb') as output:
        shutil.copyfileobj(stream, output, length=1024**2)
    require(sha(temporary) == expected, 'Evidence changed during archive copy')
    os.link(temporary, target)
    temporary.unlink()


def artifact_labels(directory, suffix):
    require(directory.is_dir(), 'Missing job evidence directory')
    labels = set()
    for path in directory.glob('*'+suffix):
        name = path.name[:-len(suffix)]
        require(LABEL.fullmatch(name) is not None, 'Bounded ASCII alphanumeric/dash job label required')
        labels.add(name)
    return labels


def archive(root, campaign):
    data_root, run_root = Path(root).resolve(strict=True), Path(campaign).resolve(strict=True)
    require(data_root.is_dir() and run_root.is_dir() and data_root != run_root and
            data_root.is_relative_to(run_root), 'RUN_ROOT must contain DATA_ROOT for cumulative disk accounting')
    cutoff = data_root/'campaign/LEDGER_CUTOFF.json'
    require(not cutoff.exists() and not cutoff.is_symlink() and
            not cutoff.with_name(cutoff.name+'.partial').exists(), 'Preserve original ledger cutoff')
    pools_path = contained(run_root, 'RESOURCE_POOLS.json')
    pools = read(pools_path)
    require(pools['schema'] == 'cost-matched-resource-pools-v1' and
            pools['guard']['path'] == 'campaign/resource_guard.py' and
            HASH.fullmatch(pools['guard']['sha256']) is not None,
            'Canonical RESOURCE_POOLS guard path/schema required')
    pools_sha = sha(pools_path)
    guard = contained(run_root, pools['guard']['path'])
    require(sha(guard) == pools['guard']['sha256'], 'Frozen guard bytes differ')
    contained(run_root, 'campaign/run_guarded_job.py')
    plans = artifact_labels(run_root/'jobs', '.PLAN.json')
    require(plans == artifact_labels(run_root/'jobs', '.FINAL.json') ==
            artifact_labels(run_root/'logs', '.resources.json') == artifact_labels(run_root/'logs', '.log'),
            'Unfinished/orphan PLAN, FINAL, receipt or log blocks scientific cutoff')
    copies = {}
    def retain(source_name, target_name, expected=None):
        source = contained(run_root, source_name)
        actual = sha(source)
        require(expected is None or actual == expected, 'Frozen evidence digest differs: '+source_name)
        target = destination(data_root, target_name)
        require(target_name not in copies or copies[target_name][2] == actual, 'Conflicting archive source identity')
        require(not target.exists() or sha(target) == actual, 'Stale archived evidence differs: '+target_name)
        copies[target_name] = source, target, actual
        return actual
    retain('RESOURCE_POOLS.json', 'campaign/RESOURCE_POOLS.json', pools_sha)
    retain('campaign/resource_guard.py', 'campaign/resource_guard.py', pools['guard']['sha256'])
    retain('campaign/run_guarded_job.py', 'campaign/run_guarded_job.py')
    entries = []
    for label in sorted(plans):
        plan_name, final_name = 'jobs/'+label+'.PLAN.json', 'jobs/'+label+'.FINAL.json'
        plan = read(contained(run_root, plan_name))
        final = read(contained(run_root, final_name))
        require(plan['label'] == label and PurePosixPath(plan_name).name == plan['label']+'.PLAN.json',
                'PLAN filename/label differs')
        require(final['label'] == label and final['phase'] == plan['phase'] and
                plan['phase'] in pools['phases'] and plan['pools_sha256'] == pools_sha,
                'PLAN/FINAL resource pool identity differs')
        dispatcher = plan.get('dispatcher_archive_path', 'campaign/run_guarded_job.py')
        dispatcher_sha = plan['dispatcher_sha256']
        require(isinstance(dispatcher_sha, str) and HASH.fullmatch(dispatcher_sha) is not None,
                'Dispatcher SHA-256 spelling')
        require(dispatcher == 'campaign/run_guarded_job.py' or
                dispatcher == 'campaign/dispatcher_versions/'+dispatcher_sha+'.py', 'Canonical dispatcher archive path')
        entry = dict(label=label, phase=plan['phase'])
        for kind, name in (('plan', plan_name), ('final', final_name),
                           ('receipt', 'logs/'+label+'.resources.json'), ('log', 'logs/'+label+'.log')):
            entry[kind+'_path'] = 'campaign/'+name
            entry[kind+'_sha256'] = retain(name, 'campaign/'+name)
        require(final['receipt_sha256'] == entry['receipt_sha256'], 'FINAL receipt digest differs')
        entry['dispatcher_path'] = dispatcher
        entry['dispatcher_sha256'] = retain(dispatcher, dispatcher, dispatcher_sha)
        entries.append(entry)
    # Check the entire source graph before changing any archive destination.
    for source, target, expected in copies.values():
        copy_evidence(source, target, expected)
    result = dict(schema='cost-matched-ledger-v1', cutoff_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        resource_pools=dict(path='campaign/RESOURCE_POOLS.json', sha256=pools_sha), jobs=entries,
        disk_accounting_scope=DISK_SCOPE, data_root_relative_to_run_root=str(data_root.relative_to(run_root)),
        scope='All scientific jobs through cutoff, including failures; later analysis/report/review/packaging recorded separately')
    write_new(cutoff, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='DATA_ROOT contained within RUN_ROOT')
    parser.add_argument('--campaign', type=Path, required=True, help='Cumulative RUN_ROOT containing jobs/, logs/, campaign/')
    args = parser.parse_args()
    archive(args.root, args.campaign)
