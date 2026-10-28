"""Carry identical frozen orders across the row-packaging-only source fix."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')
sys.path.insert(0, str(ROOT / 'numerical_audit/scalable_ordering'))
import common as k


def main():
    cpu, wall = time.process_time(), time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    old, new = args.prior.resolve(strict=True), args.output.resolve()
    k.require(not new.exists() and not new.is_relative_to(old) and not old.is_relative_to(new), 'Fresh separate output required')
    prior_path = old / 'ALL_ORDERS_FROZEN.json'
    prior = k.read(prior_path)
    k.require(len(prior['constructors']) == 24 and len(prior['models']) == 4, 'Complete prior order freeze')
    k.require(not list((old / 'rows').glob('*.json')) and not (old / 'RUN_COMPLETE.json').exists(), 'Only a zero-result-row interrupted attempt may use this amendment')
    head, sources = k.committed_sources(), k.source_map()
    changed = {name for name in set(prior['sources']) | set(sources) if prior['sources'].get(name) != sources.get(name)}
    k.require(changed == {'numerical_audit/scalable_ordering/study.py'}, 'Only the row-packaging source may change')
    k.require(prior['environment'] == k.environment(), 'Environment must remain unchanged')
    changed_paths = set(subprocess.check_output(['git', '-C', str(ROOT), 'diff', '--name-only', prior['source_commit'], head], text=True).splitlines())
    k.require(changed_paths == {'numerical_audit/scalable_ordering/study.py', 'numerical_audit/scalable_ordering/test_pipeline.py'}, 'Unexpected addressing delta')
    for name, digest in prior['sources'].items():
        k.require(k.sha(old / 'source_snapshot' / name) == digest, 'Prior source snapshot changed')
    for name, digest in prior['constructors'].items():
        k.require(k.sha(old / name) == digest, 'Prior frozen order changed')
    for case, digest in prior['models'].items():
        k.require(k.sha(old / 'chemistry/inputs' / case / 'model.json') == digest, 'Prior model changed')
    prior_reference_files = sorted(str(p.relative_to(old)) for p in (old / 'references').glob('*.npz'))
    prior_traces = {str(p.relative_to(old)): k.sha(p) for p in (old / 'traces').glob('*.gz')}
    new.mkdir(parents=True)
    shutil.copytree(old / 'chemistry', new / 'chemistry')
    shutil.copytree(old / 'constructors', new / 'constructors')
    shutil.copyfile(prior_path, new / 'PRIOR_ALL_ORDERS_FROZEN.json')
    for name in sources:
        target = new / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    for name, digest in prior['constructors'].items():
        k.require(k.sha(new / name) == digest, 'Copied order differs')
    frozen = dict(prior)
    frozen.update(source_commit=head, sources=sources,
        original_order_source_commit=prior['source_commit'],
        prior_order_freeze_sha256=k.sha(prior_path),
        note='All 24 original orders remain byte-identical and were frozen before any validation. '
             'This source amendment removes a duplicate row metadata keyword only. '
             'The prior attempt created one reference and one trace but evaluated no ordering error or result row. '
             'No constructor was rerun and no order or scientific parameter was retuned.')
    k.write(new / 'ALL_ORDERS_FROZEN.json', frozen)
    k.write(new / 'SOURCE_AMENDMENT.json', dict(status='completed',
        purpose='Packaging-only correction: cap_multiple is supplied once by the compilation result.',
        prior_study=str(old), prior_source_commit=prior['source_commit'], new_source_commit=head,
        prior_order_freeze_sha256=k.sha(prior_path), new_order_freeze_sha256=k.sha(new / 'ALL_ORDERS_FROZEN.json'),
        changed_scientific_sources=sorted(changed), original_constructors=prior['constructors'],
        same_all_24_order_hashes=True, retraining=False, retuning=False,
        ordering_errors_evaluated_in_prior_attempt=False, prior_result_rows=0,
        prior_references_retained_at_old_path=prior_reference_files,
        prior_trace_hashes_retained_at_old_path=prior_traces,
        prior_attempt_costs_retained_separately=True,
        amendment_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        CPU_seconds=time.process_time() - cpu, wall_seconds=time.monotonic() - wall))
    print('SOURCE_AMENDMENT_PASS same24orders no_retraining no_error_feedback', flush=True)


if __name__ == '__main__':
    main()
