"""Archive immutable scientific runs and completed receipts without executing them.

The current packaging attempt is outside its own completed-receipt cutoff.
Task/handoff records and reviewer packets/reports are not inputs.
"""
import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import tarfile


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            value.update(block)
    return value.hexdigest()


def package(task, output):
    task, output = task.resolve(strict=True), output.resolve(strict=True)
    runs = ('main-a1', 'preflight-a1', 'preflight-a2', 'preflight-a3')
    archive_directory = output/'archives'
    assert not archive_directory.exists(), 'Fresh archive destination required'
    collections = {}
    counts = {}
    for name in runs:
        source = task/'runs'/name
        pairs = []
        for folder in ('data', 'campaign', 'validation'):
            base = source/folder
            if base.exists():
                pairs.extend((path, Path(folder)/path.relative_to(base))
                             for path in sorted(base.rglob('*')) if path.is_file())
        for pattern in ('RESOURCE_POOLS.json', 'CAMPAIGN_FREEZE_A.json', 'PREFLIGHT_ADMISSION*.json',
                        'PREFLIGHT_REVIEW_LIMITATIONS*.json'):
            pairs.extend((path, Path(path.name)) for path in sorted(source.glob(pattern)))
        complete = []
        for final in sorted((source/'jobs').glob('*.FINAL.json')):
            label = final.name.removesuffix('.FINAL.json')
            for spelling in ('jobs/'+label+'.PLAN.json', 'jobs/'+label+'.FINAL.json',
                             'logs/'+label+'.resources.json', 'logs/'+label+'.log'):
                path = source/spelling
                assert path.is_file(), 'Incomplete retained job: '+spelling
                pairs.append((path, Path(spelling)))
            complete.append(label)
        counts[name] = dict(completed_jobs=complete,
            incomplete_job_labels=sorted(path.name.removesuffix('.PLAN.json')
                for path in (source/'jobs').glob('*.PLAN.json')
                if not path.with_name(path.name.replace('.PLAN.json', '.FINAL.json')).exists()))
        for source_file, spelling in pairs:
            assert not source_file.is_symlink() and source_file.resolve().is_relative_to(source)
            assert '__pycache__' not in spelling.parts and source_file.stat().st_size <= 32*1024**2
        collections[name] = pairs
    collections['operators'] = []
    for name in ('run_panel.py', 'test_controller.py', 'run_guarded_job.py',
                 'PREFLIGHT_CUMULATIVE_AMENDMENT.a2.json', 'PREFLIGHT_CUMULATIVE_AMENDMENT.a3.json'):
        source_file = task/'coordination'/name
        collections['operators'].append((source_file, Path(name)))
    copies = [(path, root+'/'+str(spelling)) for root, pairs in collections.items() for path, spelling in pairs]
    assert len({target for _, target in copies}) == len(copies), 'Duplicate output path'
    assert len(copies) <= 10000 and sum(path.stat().st_size for path, _ in copies) <= 128*1024**2
    assert not (output/'COPY_MANIFEST.a1.json').exists()
    archive_directory.mkdir()
    archives = []
    for root, pairs in collections.items():
        target = archive_directory/(root+'.tar.gz')
        members = []
        with target.open('xb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w|', format=tarfile.PAX_FORMAT) as archive:
                for source_file, spelling in sorted(pairs, key=lambda pair: str(pair[1])):
                    before = digest(source_file)
                    size = source_file.stat().st_size
                    info = tarfile.TarInfo(str(spelling))
                    info.size, info.mode, info.mtime = size, 0o644, 0
                    with source_file.open('rb') as stream:
                        archive.addfile(info, stream)
                    assert digest(source_file) == before, 'Source changed during archive'
                    members.append(dict(path=str(spelling), source=str(source_file), bytes=size, sha256=before))
        expected = {row['path']: row for row in members}
        observed = set()
        with tarfile.open(target, mode='r|gz') as archive:
            for member in archive:
                assert member.isfile() and member.name in expected and member.name not in observed
                row = expected[member.name]
                assert member.size == row['bytes']
                value = hashlib.sha256()
                with archive.extractfile(member) as stream:
                    for block in iter(lambda: stream.read(1024**2), b''):
                        value.update(block)
                assert value.hexdigest() == row['sha256'], 'Archived member differs'
                observed.add(member.name)
        assert observed == set(expected), 'Incomplete archive'
        archives.append(dict(path=str(target.relative_to(output)), root=root,
            bytes=target.stat().st_size, sha256=digest(target), members=members,
            total_uncompressed_bytes=sum(row['bytes'] for row in members)))
    value = dict(schema='cost-matched-evidence-archive-v1',
        UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scope='Complete DATA_ROOTs, operator/validation snapshots and every completed job at packaging cutoff; no current-attempt FINAL exists yet',
        scientific_cutoff='main-a1/data/campaign/LEDGER_CUTOFF.json; later receipts stay outside that immutable ledger',
        counts=counts, archives=archives,
        total_member_bytes=sum(row['total_uncompressed_bytes'] for row in archives),
        total_archive_bytes=sum(row['bytes'] for row in archives),
        package_source_sha256=digest(Path(__file__)))
    with (output/'COPY_MANIFEST.a1.json').open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(files=len(copies), total_member_bytes=value['total_member_bytes'],
                         archive_bytes=value['total_archive_bytes'],
                         manifest_sha256=digest(output/'COPY_MANIFEST.a1.json'))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    package(args.task_root, args.output)
