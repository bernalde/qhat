"""Root-only, one-job dispatcher with append-only cumulative phase accounting."""
import argparse
import datetime
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--cpu', type=float, required=True)
    parser.add_argument('--wall', type=int, required=True)
    parser.add_argument('--cwd', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    root = args.root.resolve()
    pools_path = root / 'RESOURCE_POOLS.json'
    pools = json.loads(pools_path.read_text())
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.label.replace('-', '').isalnum() or not math.isfinite(args.cpu) or min(args.cpu, args.wall) <= 0:
        raise ValueError('Nonempty command, unique label and positive caps required')
    caps = pools['phases'][args.phase]
    guard = Path(pools['guard']['path'])
    if digest(guard) != pools['guard']['sha256']:
        raise ValueError('Frozen guard changed')
    jobs = root / 'jobs'
    logs = root / 'logs'
    jobs.mkdir(exist_ok=True)
    logs.mkdir(exist_ok=True)
    with (root / 'DISPATCH.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        used_cpu = used_wall = 0.0
        for path in sorted(jobs.glob('*.PLAN.json')):
            plan = json.loads(path.read_text())
            final = path.with_name(path.name.replace('.PLAN.json', '.FINAL.json'))
            if not final.exists():
                raise RuntimeError('Unreconciled attempt blocks another launch: ' + str(path))
            result = json.loads(final.read_text())
            if plan['phase'] == args.phase:
                used_cpu += result['charged_CPU_seconds']
                used_wall += result['absolute_guarded_wall_seconds']
        if used_cpu + args.cpu > caps['CPU_seconds'] or used_wall + args.wall > caps['wall_seconds']:
            raise RuntimeError('Attempt does not fit remaining frozen phase allowance')
        plan_path = jobs / (args.label + '.PLAN.json')
        save(plan_path, dict(phase=args.phase, label=args.label,
            attempt_CPU_seconds=args.cpu, attempt_wall_seconds=args.wall,
            prior_phase_CPU_seconds=used_cpu, prior_phase_wall_seconds=used_wall,
            cwd=str(args.cwd.resolve()), command=command,
            pools_sha256=digest(pools_path), dispatcher_sha256=digest(__file__),
            started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        run = subprocess.run([sys.executable, str(guard), '--label', args.label,
            '--logs', str(logs), '--cpu-seconds', str(args.cpu),
            '--wall-seconds', str(args.wall), '--disk-mib', str(pools['output_MiB']),
            '--output', str(root), '--', *command], cwd=args.cwd)
        receipt_path = logs / (args.label + '.resources.json')
        if not receipt_path.is_file():
            raise RuntimeError('Missing receipt; keep reserved attempt unreconciled')
        receipt = json.loads(receipt_path.read_text())
        charged = receipt.get('final_accounting', {}).get('CPU_seconds')
        cpu = args.cpu if charged is None else float(charged)
        log = logs / (args.label + '.log')
        factors = {'h': 3600, 'min': 60, 's': 1, 'ms': .001,
                   'us': .000001, 'µs': .000001, 'ns': .000000001}
        for summary in re.findall(r'^CPU time consumed: (.+)$', log.read_text(), re.MULTILINE):
            tokens = re.findall(r'([0-9.]+)(h|min|ms|us|µs|ns|s)', summary)
            if tokens:
                cpu = max(cpu, sum(float(value) * factors[unit] for value, unit in tokens) + .001)
        # Failed services commonly remain available after systemd-run exits.
        snapshot = subprocess.run(['systemctl', '--user', 'show', receipt['unit'],
            '-p', 'CPUUsageNSec', '-p', 'MemoryPeak', '-p', 'Result',
            '-p', 'ActiveState', '-p', 'SubState'], text=True, capture_output=True)
        fields = dict(row.split('=', 1) for row in snapshot.stdout.splitlines() if '=' in row)
        if fields.get('CPUUsageNSec', '').isdigit():
            cpu = max(cpu, int(fields['CPUUsageNSec']) / 1e9)
        save(jobs / (args.label + '.FINAL.json'), dict(phase=args.phase, label=args.label,
            exit_code=run.returncode, receipt_sha256=digest(receipt_path),
            charged_CPU_seconds=cpu, absolute_guarded_wall_seconds=receipt['elapsed'],
            final_service_snapshot=fields, final_service_query_exit=snapshot.returncode,
            accounting_scope='max guard process/service measurement and retained final service CPU; missing CPU charges full reserved allowance',
            completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        return run.returncode


if __name__ == '__main__':
    raise SystemExit(main())
