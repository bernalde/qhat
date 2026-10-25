"""Altered-model representative end-to-end check, including full extraction."""
import argparse
from pathlib import Path
import tarfile
import time
import common as k
import study
import verify
import delivery


def run(root):
    k.require(not root.exists(), 'Fresh preflight root required')
    root.mkdir(parents=True)
    output = root / 'results'
    began, wall = time.process_time(), time.monotonic()
    study.prepare(output, preflight=True)
    study.run(output)
    verified = verify.verify(output)
    packed = root / 'delivery'
    delivery.pack(output, packed)
    delivery.check(packed)
    archive = root / 'delivery.tar.gz'
    with tarfile.open(archive, 'w:gz', compresslevel=1) as stream:
        stream.add(packed, arcname='delivery')
    extraction = root / 'extracted'
    extraction.mkdir()
    with tarfile.open(archive, 'r:gz') as stream:
        for member in stream.getmembers():
            k.safe(member.name)
            k.require(member.isdir() or member.isfile(), 'Unexpected archive member type')
        stream.extractall(extraction, filter='data')
    delivery.check(extraction / 'delivery')
    # The full independent numerical replay already covers every candidate;
    # extraction checks all identical bytes against that verified inventory.
    total_bytes = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    record = dict(status='passed', CPU_seconds=time.process_time() - began,
                  wall_seconds=time.monotonic() - wall, output_bytes=total_bytes,
                  complete_candidates=verified['complete'], candidates=verified['candidates'],
                  verified_target_rows=verified['target_rows'],
                  full_archive_extraction_authenticated=True,
                  independent_replay_before_extract=True,
                  main_estimation_scope='Measured high-mode altered-model full path;main3repetitionsand3models',
                  extrapolated_main_candidate_multiplier=1296 / 176,
                  inherited_reference_costs_not_new_measurements=True)
    k.write(root / 'PREFLIGHT.json', record)
    print('PRACTICAL_PREFLIGHT_PASSED', record, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve())
