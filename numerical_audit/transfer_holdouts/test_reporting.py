"""Truthful resource ledgers and refusal to report unverified scientific data."""
import copy
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pytest
import reporting


def ledger():
    pools = {name: dict(CPU_seconds=limits[0], wall_seconds=limits[1], disk_MiB=limits[2],
                       attempts=[], consumed_CPU_seconds=0., consumed_wall_seconds=0.)
             for name, limits in reporting.POOLS.items()}
    receipt = dict(exit_code=0, elapsed=3., final_accounting=dict(CPU_seconds=1., output_bytes=100),
                   settings=dict(CPUQuota='25%', Nice='19', CPUWeight='5', IOWeight='10',
                                 MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0'),
                   environment={name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                     'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
    raw = json.dumps(receipt, indent=2) + '\n'
    pools['chemistry.tests'].update(attempts=[dict(receipt='unit-fixture.resources.json', CPU_seconds=1., wall_seconds=3., exit_code=0)],
                                    consumed_CPU_seconds=1., consumed_wall_seconds=3.)
    return dict(pools=pools, receipts={'unit-fixture.resources.json': dict(text=raw, sha256=hashlib.sha256(raw.encode()).hexdigest())},
                snapshot_scope='Synthetic unit fixture;no real study or review attestation')


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def test_resource_totals_bind_every_original_receipt_without_resetting_pools(tmp_path):
    path, value = tmp_path / 'ledger.json', ledger()
    save(path, value)
    _, totals = reporting.resources(path)
    assert totals['chemistry.tests']['CPU_seconds'] == 1.
    value['pools']['chemistry.tests']['CPU_seconds'] = 240
    save(path, value)
    with pytest.raises(ValueError, match='Original transfer phase ceiling changed'):
        reporting.resources(path)


def test_undercharged_attempt_and_unlisted_receipt_are_rejected(tmp_path):
    path, value = tmp_path / 'ledger.json', ledger()
    original = copy.deepcopy(value)
    value['pools']['chemistry.tests']['attempts'][0]['CPU_seconds'] = .5
    save(path, value)
    with pytest.raises(ValueError, match='Resource attempt undercharged'):
        reporting.resources(path)
    original['receipts']['extra.resources.json'] = original['receipts']['unit-fixture.resources.json']
    save(path, original)
    with pytest.raises(ValueError, match='Uncharged or omitted raw resource receipt'):
        reporting.resources(path)


def test_proof_memory_approval_cannot_be_applied_to_numeric_jobs(tmp_path):
    path, value = tmp_path / 'ledger.json', ledger()
    row = value['receipts']['unit-fixture.resources.json']
    receipt = json.loads(row['text'])
    receipt['settings']['MemoryMax'] = '6G'
    row['text'] = json.dumps(receipt, indent=2) + '\n'
    row['sha256'] = hashlib.sha256(row['text'].encode()).hexdigest()
    save(path, value)
    with pytest.raises(ValueError, match='Numerical worker resource profile differs'):
        reporting.resources(path)


def test_unverified_producer_cannot_create_a_report(tmp_path):
    verification = tmp_path / 'verification'
    verification.mkdir()
    save(verification / 'VERIFIED.json', dict(status='producer_complete'))
    destination = tmp_path / 'report'
    with pytest.raises(ValueError, match='Complete independent audit required'):
        reporting.report(tmp_path / 'not-read', verification, tmp_path / 'not-read-ledger', destination)
    assert not destination.exists()
