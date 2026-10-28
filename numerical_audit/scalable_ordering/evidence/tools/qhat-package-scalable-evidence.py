"""Package completed scientific evidence without private coordination records.

Run only after the independent replay and report have completed, under the
campaign's reporting resource allowance. Archives retain source bytes inertly.
"""
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile

N = Path('/tmp/qhat-scalable-ordering-2026-10-27')
Q = N / 'worktree'
DEST = Q / 'numerical_audit/scalable_ordering/evidence'
MAX_BYTES = 400 * 1024**2
SECRET = re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{25,}|github_pat_[A-Za-z0-9_]{30,}|sk-[A-Za-z0-9]{35,})\b')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def inventory(roots):
    result = {}
    for root in roots:
        require(root.exists() and root.is_relative_to(N), 'Missing/outside retained input')
        paths = sorted(root.rglob('*')) if root.is_dir() else [root]
        for path in paths:
            require(not path.is_symlink(), 'Symlink in retained evidence')
            if path.is_dir():
                continue
            require(path.is_file(), 'Nonregular retained evidence')
            name = path.relative_to(N).as_posix()
            require(name not in result and '..' not in Path(name).parts, 'Duplicate/unsafe retained path')
            size = path.stat().st_size
            require(size <= 32 * 1024**2, 'Evidence member too large')
            if path.suffix in ('.json', '.py', '.md', '.log', '.txt'):
                require(not SECRET.search(path.read_bytes()), 'Possible credential in ' + name)
            result[name] = dict(bytes=size, sha256=sha(path))
    require(len(result) < 10000 and sum(v['bytes'] for v in result.values()) <= MAX_BYTES,
            'Retained evidence inventory exceeds bounds')
    return result


def archive(name, roots):
    members = inventory(roots)
    path = DEST / name
    with path.open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=1, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w|', format=tarfile.PAX_FORMAT) as target:
                for relative, entry in sorted(members.items()):
                    source = N / relative
                    info = tarfile.TarInfo(relative)
                    info.size = entry['bytes']
                    info.mode = 0o644
                    info.mtime = 0
                    info.uid = info.gid = 0
                    with source.open('rb') as stream:
                        target.addfile(info, stream)
                    require(sha(source) == entry['sha256'], 'Source changed during archive')
    checked = {}
    with tarfile.open(path, 'r:gz') as source:
        for member in source:
            require(member.isfile() and member.name in members and member.name not in checked,
                    'Unexpected/duplicate archive member')
            require(member.size == members[member.name]['bytes'], 'Archive member size differs')
            value = hashlib.sha256()
            with source.extractfile(member) as stream:
                for block in iter(lambda: stream.read(1024**2), b''):
                    value.update(block)
            require(value.hexdigest() == members[member.name]['sha256'], 'Archive member bytes differ')
            checked[member.name] = True
    require(set(checked) == set(members), 'Archive member inventory differs')
    write(DEST / (name + '.inventory.json'), dict(archive=name, bytes=path.stat().st_size,
        sha256=sha(path), members=members, member_bytes_verified=True,
        storage='Regular-file tar members; deterministic zero metadata; source files preserved byte-for-byte'))


def main():
    study = N / 'study-a2'
    run = read(study / 'RUN_COMPLETE.json')
    independent = read(N / 'verification-a1/INDEPENDENT_VERIFIED.json')
    tables = read(N / 'report-a1/report_tables.json')
    require(len(run['rows']) == 96 and len(run['references']) == 8 and
            independent['status'] == 'passed' and len(independent['pairs']) == 240 and
            tables['completion']['independent_status'] == 'passed', 'Incomplete scientific payload')
    require(not DEST.exists(), 'Fresh evidence destination required; preserve earlier packages')
    DEST.mkdir(parents=True)
    archive('main-study.tar.gz', [study])
    task = read(N / 'TASK.json')
    retained = [N / name for name in ('inputs-a2', 'preflight-a1', 'preflight-a2', 'study-a1',
        'pauli-setup-probe-a1', 'stream-probe-a1')]
    retained += sorted(N.glob('*PROFILE*.json'))
    retained += [N / name for name in ('PAULI_PAIR_COUNT.a1.json', 'PREFLIGHT_CONTINUATION.a2.a1.json',
        'MAIN_EXECUTION_FROZEN.json', 'PHASE_ACCOUNTING_BEFORE_MAIN.json', 'FULL_ANALYSIS_TESTS.a1.json',
        'INDEPENDENT_ALLOCATION.a1.json')]
    retained += sorted(N.glob('TEST_RESULT*.json')) + sorted(N.glob('TESTED_SOURCE_INPUTS*.json'))
    # The current packaging log is still being written. Include terminal jobs
    # from the ledger only; do not archive an active log as finished evidence.
    for row in task['guarded_attempts']:
        receipt_path = Path(row['receipt'])
        retained += [receipt_path, Path(read(receipt_path)['log'])]
    archive('retained-attempts.tar.gz', retained)
    shutil.copytree(N / 'report-a1', DEST / 'report')
    shutil.copyfile(N / 'verification-a1/INDEPENDENT_VERIFIED.json', DEST / 'INDEPENDENT_VERIFIED.json')
    shutil.copyfile(study / 'SUMMARY.json', DEST / 'SUMMARY.json')
    costs = []
    for row in task['guarded_attempts']:
        receipt = read(Path(row['receipt']))
        costs.append(dict(label=Path(row['receipt']).name.removesuffix('.resources.json'),
            receipt_sha256=row['sha256'], CPU_seconds=row['CPU_seconds'], wall_seconds=row['wall_seconds'],
            exit_code=row['exit_code'], termination=row['termination'], budget=receipt['budget'],
            source_receipt_archive='retained-attempts.tar.gz',
            archive_member='logs/' + Path(row['receipt']).name,
            process_measurements=receipt.get('process_measurements', {})))
    write(DEST / 'CAMPAIGN_COSTS.json', dict(scope='Completed guarded attempts through report generation; '
        'includes failed attempts. Coordinator/metadata overhead and remote review/CI are not included. '
        'Later packaging, publication and review costs are outside this snapshot.',
        CPU_seconds=sum(r['CPU_seconds'] for r in costs), wall_seconds=sum(r['wall_seconds'] for r in costs),
        whole_CPU_limit=task['resources']['whole_CPU_seconds'],
        whole_wall_limit=task['resources']['whole_wall_seconds'],
        machine=task['resources']['machine'], attempts=costs))
    tools = DEST / 'tools'
    tools.mkdir()
    for name in ('qhat-scalable-render-report.py', 'qhat-package-scalable-evidence.py',
                 'qhat-amend-main-source.py', 'qhat-resume-preflight.py', 'qhat-profile-current.py',
                 'qhat-count-pauli-pairs.py', 'qhat-pauli-setup-probe.py', 'qhat-probe-stream.py',
                 'qhat-profile-cubic-seed.py', 'qhat-scalable-run-tests.py',
                 'qhat-scalable-training-before-optimization.py', 'claude-guarded.py'):
        source = Path('/tmp') / name
        require(source.is_file(), 'Missing retained producer tool: ' + name)
        shutil.copyfile(source, tools / name)
    payload = {p.relative_to(DEST).as_posix(): dict(bytes=p.stat().st_size, sha256=sha(p))
               for p in sorted(DEST.rglob('*')) if p.is_file()}
    write(DEST / 'PACKAGE_MANIFEST.json', dict(schema='qhat-scalable-ordering-evidence-v1',
        scientific_source_commit=read(study / 'ALL_ORDERS_FROZEN.json')['source_commit'],
        original24order_freeze_sha256=sha(N / 'study-a1/ALL_ORDERS_FROZEN.json'),
        effective24order_freeze_sha256=sha(study / 'ALL_ORDERS_FROZEN.json'),
        original24orders_reused=True, independent_candidates=96, independent_pairs=240,
        payload=payload, excludes=['Private TASK/coordination records', 'model-client account metadata'],
        note='Additional prose and review reports may follow; this manifest authenticates this scientific payload.'))
    print(json.dumps(dict(status='packaged', destination=str(DEST), files=len(payload))))


if __name__ == '__main__':
    main()
