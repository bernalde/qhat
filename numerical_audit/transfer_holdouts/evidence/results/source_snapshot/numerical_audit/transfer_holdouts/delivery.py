"""Complete transfer payload, deterministic archive, extraction and Git checks."""
import gzip
import hashlib
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import reporting

audit = reporting.audit
HERE = Path(__file__).resolve().parent
audit.require(Path(reporting.__file__).resolve() == HERE / 'reporting.py', 'Wrong delivery reporter binding')


def file_hashes(root):
    return {name: audit.sha(Path(root) / name) for name in sorted(audit.inventory(root))}


def pack(results, verification, resource_ledger, destination):
    results, verification, destination = Path(results), Path(verification), Path(destination)
    verified = reporting.require_verified(results, verification)
    audit.require(not destination.exists(), 'Fresh delivery required;retain earlier attempts')
    audit.inventory(verification)
    destination.mkdir(parents=True)
    shutil.copytree(results, destination / 'results')
    shutil.copytree(verification, destination / 'verification')
    reporting.report(destination / 'results', destination / 'verification', resource_ledger, destination / 'report')
    audit.write(destination / 'SOURCE_PROVENANCE.json', dict(
        producer_source_commit=verified['producer_source_commit'],
        producer_sources=verified['producer_sources'], current_verifier_sources=verified['current_verifier_sources'],
        packaging_sources={name: audit.sha(HERE / name) for name in ('reporting.py', 'delivery.py')},
        practical_dependency=audit.read(results / 'PRACTICAL_INTERFACE.json'),
        required_publication_gate='Coordinator verifies actual numerical prerequisite merge and exact source/model bytes before publication',
        review_status='PairedAstra/Argonconvergenceandowningpublication/separatemerge/CIremainrequired',
        scope='Completefrozenpanel,allreferences/constructors/traces/rows/censors/pairs,andindependentreceipts;private source-review reports excluded'))
    audit.write(destination / 'MANIFEST.json', dict(files=file_hashes(destination),
        circuits=verified['unique_circuits'], logical_circuits=verified['logical_circuits'],
        comparisons=verified['comparisons'], all_payload_files_included=True))


def git_blobs(root, commit):
    root = Path(root).resolve()
    audit.require(root.is_relative_to(audit.ROOT.resolve()), 'Git delivery must be inside its repository')
    prefix = str(root.relative_to(audit.ROOT.resolve())) + '/'
    raw = subprocess.check_output(['git', '-C', str(audit.ROOT), 'ls-tree', '-rz', commit, '--', prefix])
    stored = {}
    for entry in raw.split(b'\0'):
        if entry:
            metadata, name = entry.split(b'\t', 1)
            mode, kind, oid = metadata.split()
            audit.require(kind == b'blob' and mode in (b'100644', b'100755'), 'Unexpected Git payload type')
            stored[name.decode()[len(prefix):]] = oid.decode()
    audit.require(set(stored) == audit.inventory(root), 'Complete Git delivery inventory differs,including binary archives')
    for name, expected in stored.items():
        path = root / name
        digest = hashlib.sha1(b'blob ' + str(path.stat().st_size).encode() + b'\0')
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        audit.require(digest.hexdigest() == expected, 'Delivered Git blob changed: ' + name)


def check(root, commit=None):
    root = Path(root)
    manifest = audit.read(root / 'MANIFEST.json')
    audit.require(audit.inventory(root) == set(manifest['files']) | {'MANIFEST.json'}, 'Strict complete delivery inventory differs')
    for name, digest in manifest['files'].items():
        audit.require(audit.sha(root / audit.safe(name)) == digest, 'Delivered scientific/report/verification bytes changed')
    verified = reporting.require_verified(root / 'results', root / 'verification')
    audit.require((manifest['circuits'], manifest['logical_circuits'], manifest['comparisons']) ==
                  (verified['unique_circuits'], verified['logical_circuits'], verified['comparisons']), 'Delivery panel incomplete')
    with tempfile.TemporaryDirectory(prefix='transfer-report-check-', dir=root.parent) as temporary:
        rebuilt = Path(temporary) / 'report'
        reporting.report(root / 'results', root / 'verification', root / 'report/EXECUTION_RESOURCES_SNAPSHOT.json', rebuilt)
        audit.require(audit.inventory(rebuilt) == audit.inventory(root / 'report'), 'Deterministic report file inventory differs')
        for name in audit.inventory(rebuilt):
            audit.require(audit.sha(rebuilt / name) == audit.sha(root / 'report' / name), 'Rebuilt report differs')
    provenance = audit.read(root / 'SOURCE_PROVENANCE.json')
    audit.require(provenance['producer_sources'] == verified['producer_sources'] and
                  provenance['current_verifier_sources'] == verified['current_verifier_sources'], 'Source provenance differs')
    for name, digest in provenance['packaging_sources'].items():
        audit.require(audit.sha(HERE / name) == digest, 'Trusted packaging source changed')
    if commit is not None:
        git_blobs(root, commit)
    return dict(status='passed', files=len(manifest['files']) + 1, logical_circuits=verified['logical_circuits'])


def archive(root, path):
    root, path = Path(root), Path(path)
    check(root)
    audit.require(not path.exists(), 'Fresh deterministic archive required')
    with path.open('xb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=1, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w') as stream:
                for entry in [root] + sorted(root.rglob('*')):
                    member = stream.gettarinfo(str(entry), arcname=str(Path('delivery') / entry.relative_to(root)))
                    audit.require(member.isdir() or member.isfile(), 'Archive links/special files refused')
                    member.uid = member.gid = 0
                    member.uname = member.gname = ''
                    member.mtime = 0
                    member.mode = 0o755 if member.isdir() else 0o644
                    member.pax_headers = {}
                    if member.isfile():
                        with entry.open('rb') as contents:
                            stream.addfile(member, contents)
                    else:
                        stream.addfile(member)
    return audit.sha(path)


def extract_check(path, destination):
    path, destination = Path(path), Path(destination)
    audit.require(not destination.exists(), 'Fresh extraction required;retain existing evidence')
    destination.mkdir(parents=True)
    with tarfile.open(path, 'r:gz') as stream:
        names, total, members = set(), 0, []
        for member in stream:
            audit.safe(member.name)
            total += member.size
            audit.require(len(members) < 17024 and member.name not in names and total <= 512 * 1024**2 and
                          (member.isdir() or member.isfile()) and
                          (member.name == 'delivery' or member.name.startswith('delivery/')), 'Unsafe or oversized archive member')
            names.add(member.name)
            members.append(member)
        stream.extractall(destination, members=members, filter='data')
    return check(destination / 'delivery')
