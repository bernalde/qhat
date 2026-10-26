"""Complete byte inventories, independent-result bindings and Git blob checks."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile
import common as k
import study
from reporting import report

MAX_FILES = 16000
MAX_BYTES = 512 * 1024 ** 2


def inventory(root):
    paths = list(root.rglob('*'))
    k.require(not any(p.is_symlink() for p in paths), 'Delivery symlink refused')
    files = [p for p in paths if p.is_file()]
    k.require(len(files) <= MAX_FILES and sum(p.stat().st_size for p in files) <= MAX_BYTES, 'Delivery storage ceiling')
    return {str(p.relative_to(root)) for p in files}


def pack(results, destination):
    k.require(not destination.exists(), 'Preserve deliveries; fresh destination required')
    inventory(results)
    destination.mkdir(parents=True)
    shutil.copytree(results, destination / 'results')
    report(destination / 'results', destination / 'report')
    k.write(destination / 'MANIFEST.json', dict(files={n: k.sha(destination / n) for n in sorted(inventory(destination))},
            packaging_sources=k.sources(), scope='Complete exact primitive circuits,partialprefixes,rows,states,alltargetchoices,allcomparisons,source/upstreamprovenance'))


def git_blobs(delivery, commit):
    prefix = str(delivery.relative_to(k.ROOT)) + '/'
    listing = subprocess.check_output(['git', '-C', str(k.ROOT), 'ls-tree', '-rz', commit, '--', prefix])
    blobs = {}
    for entry in listing.split(b'\0'):
        if entry:
            metadata, name = entry.split(b'\t', 1)
            blobs[name.decode()[len(prefix):]] = metadata.split()[2].decode()
    k.require(set(blobs) == inventory(delivery), 'Delivered Git inventory differs, including binary archives')
    for name, expected in blobs.items():
        data = (delivery / name).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        k.require(actual == expected, 'Delivered Git blob differs: ' + name)


def check(delivery, commit=None):
    manifest = k.read(delivery / 'MANIFEST.json')
    k.require(inventory(delivery) == set(manifest['files']) | {'MANIFEST.json'}, 'Exact delivery inventory differs')
    for name, digest in manifest['files'].items():
        k.require(k.sha(delivery / k.safe(name)) == digest, 'Delivery bytes changed: ' + name)
    output = delivery / 'results'
    protocol, frozen = study.bind(output, independent=True)
    verified = k.read(output / 'VERIFIED.json')
    k.require(verified['status'] == 'passed' and verified['decision_agreement'], 'Unverified study')
    k.require(verified['producer_sources'] == protocol['sources'], 'Producing source identity differs')
    k.require(verified['schedule_freeze_sha256'] == k.sha(output / 'SCHEDULES_FROZEN.json'), 'Verification freeze differs')
    k.require((verified['curves'], verified['candidates'], verified['complete'], verified['target_rows'], verified['comparisons']) ==
              (protocol['expected_curves'], protocol['expected_candidates'], frozen['complete'],
               protocol['expected_target_rows'], protocol['expected_comparisons']), 'Verified panel incomplete')
    k.require(k.read(output / 'INDEPENDENT_SELECTIONS.json') == k.read(output / 'SELECTIONS.json'), 'Independent target choices differ')
    k.require(k.sha(output / 'INDEPENDENT_ERRORS.json') == verified['independent_errors_sha256'] and
              k.sha(output / 'INDEPENDENT_SELECTIONS.json') == verified['independent_selections_sha256'], 'Independent records changed')
    physical = {folder + '/' + n for folder in ('rows', 'states') for n in inventory(output / folder)}
    physical |= {'SELECTIONS.json', 'COMPARISONS.json', 'RUN_COMPLETE.json', 'RUN_COSTS.json'}
    k.require(set(verified['physical_artifacts']) == physical, 'Verified physical inventory differs')
    for name, digest in verified['physical_artifacts'].items():
        k.require(k.sha(output / k.safe(name)) == digest, 'Verified physical file changed')
    expected = set(frozen['artifacts']) | physical | {'SCHEDULES_FROZEN.json', 'VERIFIED.json',
                'INDEPENDENT_ERRORS.json', 'INDEPENDENT_SELECTIONS.json'}
    k.require(inventory(output) == expected, 'Unexpected or omitted result artifact')
    with tempfile.TemporaryDirectory(prefix='practical-report-check-', dir=delivery.parent) as temporary:
        rebuilt = Path(temporary) / 'report'
        report(output, rebuilt)
        k.require(inventory(rebuilt) == inventory(delivery / 'report'), 'Report inventory differs')
        for name in inventory(rebuilt):
            k.require((rebuilt / name).read_bytes() == (delivery / 'report' / name).read_bytes(), 'Deterministic report differs')
    if commit is not None:
        git_blobs(delivery, commit)
    print('PRACTICAL_DELIVERY_PASSED', len(manifest['files']), 'files;', verified['candidates'], 'candidates', flush=True)
    return len(manifest['files'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('pack', 'check'))
    parser.add_argument('--delivery', type=Path, required=True)
    parser.add_argument('--results', type=Path)
    parser.add_argument('--git-commit')
    args = parser.parse_args()
    if args.phase == 'pack':
        k.require(args.results is not None, 'Results required')
        pack(args.results.resolve(), args.delivery.resolve())
    else:
        check(args.delivery.resolve(), args.git_commit)
