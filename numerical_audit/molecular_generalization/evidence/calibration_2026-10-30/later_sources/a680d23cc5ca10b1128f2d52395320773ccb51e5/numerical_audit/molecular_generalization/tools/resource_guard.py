"""Run one owned systemd service at low priority with monitored memory headroom."""
import argparse
import datetime
import json
import math
import os
import re
from pathlib import Path
import subprocess
import time


def available():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    raise RuntimeError('MemAvailable unavailable')


def output_bytes(path):
    """Charge a named output file or all files below an output directory."""
    if path.is_file():
        return path.stat().st_size
    if not path.exists():
        return 0
    total = 0
    for entry in path.rglob('*'):
        try:
            if entry.is_file():
                total += entry.stat().st_size
        except FileNotFoundError:
            # A producer can replace a temporary file between list and stat.
            continue
    return total


def process_measurements(contents):
    measured = {}
    for field, pattern in (
            ('user_cpu_seconds', r'User time \(seconds\): ([0-9.]+)'),
            ('system_cpu_seconds', r'System time \(seconds\): ([0-9.]+)'),
            ('max_rss_kib', r'Maximum resident set size \(kbytes\): ([0-9]+)')):
        match = re.search(pattern, contents)
        if match:
            measured[field] = float(match.group(1))
    return measured


def service_cpu_seconds(contents):
    fields = dict(line.split('=', 1) for line in contents.splitlines() if '=' in line)
    cpu = fields.get('CPUUsageNSec', '')
    return int(cpu) / 1e9 if cpu.isdigit() else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('--logs', type=Path, required=True)
    parser.add_argument('--cpu-seconds', type=float, required=True)
    parser.add_argument('--wall-seconds', type=int, required=True)
    parser.add_argument('--disk-mib', type=int, required=True)
    parser.add_argument('--memory-profile', choices=('numerical', 'calibration-reference'),
                        default='numerical', help='Reference-only extension requires a declared allowance')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.label.replace('-', '').isalnum():
        raise ValueError('Invalid command or service label')
    args.logs.mkdir(parents=True, exist_ok=True)
    unit = 'fermionic-future-' + args.label + '-' + str(os.getpid())
    report = args.logs / (args.label + '.resources.json')
    if report.exists():
        raise FileExistsError(report)
    if not math.isfinite(args.cpu_seconds) or min(args.cpu_seconds, args.wall_seconds, args.disk_mib) <= 0:
        raise ValueError('Finite positive resource allowances required')
    while available() <= 5.5 * 1024**3:
        print('WAITING FOR MEMORY HEADROOM', available(), flush=True)
        time.sleep(15)
    high, maximum = (('640M', '768M') if args.memory_profile == 'calibration-reference'
                     else ('448M', '512M'))
    settings = dict(CPUQuota='25%', CPUWeight='5', IOWeight='10', Nice='19',
                    MemoryHigh=high, MemoryMax=maximum, MemorySwapMax='0',
                    TasksMax='32', MemoryAccounting='yes', CPUAccounting='yes',
                    RuntimeMaxSec=str(args.wall_seconds), TimeoutStopSec='5s')
    invocation = ['systemd-run', '--user', '--wait', '--pipe', '--unit', unit,
                  '--working-directory', str(Path.cwd())]
    for key, value in settings.items():
        invocation += ['--property', key + '=' + value]
    env = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMBA_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1', PYTHONUNBUFFERED='1',
               NUMBA_CACHE_DIR=str(args.logs / 'numba-cache'),
               MPLCONFIGDIR=str(args.logs / 'matplotlib-cache'))
    env['LC_ALL'] = 'C'
    invocation += ['/usr/bin/env'] + [k + '=' + v for k, v in env.items()] + ['/usr/bin/time', '--verbose'] + command
    logpath = args.logs / (args.label + '.log')
    if logpath.exists():
        raise FileExistsError(logpath)
    started = time.monotonic()
    samples = []
    frozen = False
    frozen_at = None
    frozen_seconds = 0.0
    termination = None
    guard_error = None
    cleanup_errors = []
    child = None
    supervisor_status = None
    status = 1
    try:
        with logpath.open('x') as log:
            child = subprocess.Popen(invocation, stdout=log, stderr=subprocess.STDOUT)
            while child.poll() is None:
                headroom = available()
                sample = {'elapsed': time.monotonic() - started,
                          'available_bytes': headroom, 'frozen': frozen}
                result = subprocess.run(['systemctl', '--user', 'show', unit,
                                         '-p', 'MemoryCurrent', '-p', 'MemoryPeak',
                                         '-p', 'CPUUsageNSec'], capture_output=True, text=True, check=False)
                sample['service'] = result.stdout.strip()
                cpu = service_cpu_seconds(sample['service'])
                disk = output_bytes(args.output)
                sample['output_bytes'] = disk
                samples.append(sample)
                if cpu is not None and cpu > args.cpu_seconds:
                    termination = 'CPU budget exhausted'
                if disk > args.disk_mib * 1024**2:
                    termination = 'Disk budget exhausted'
                # Frozen time remains part of the absolute wall ceiling, as it
                # does in RuntimeMaxSec. Headroom pauses never extend the job.
                if time.monotonic() - started > args.wall_seconds:
                    termination = 'Wall budget exhausted'
                if termination:
                    subprocess.run(['systemctl', '--user', 'stop', unit], check=False)
                    break
                action = ('freeze' if headroom < 5 * 1024**3 and not frozen else
                          'thaw' if headroom > 5.5 * 1024**3 and frozen else None)
                if action:
                    transition = subprocess.run(['systemctl', '--user', action, unit], check=False)
                    if transition.returncode != 0:
                        # The transient unit can disappear after poll(). A
                        # finished child needs no transition; a live child must
                        # not continue without the requested headroom control.
                        if child.poll() is None:
                            termination = 'Headroom state transition failed: ' + action
                        break
                    now = time.monotonic()
                    if action == 'freeze':
                        frozen = True
                        frozen_at = now
                        print('FROZE OWN SERVICE FOR HEADROOM', unit, flush=True)
                    else:
                        frozen_seconds += now - frozen_at
                        frozen = False
                        frozen_at = None
                        print('RESUMED OWN SERVICE', unit, flush=True)
                time.sleep(3)
    except Exception as exception:
        guard_error = type(exception).__name__ + ': ' + str(exception)
        termination = termination or 'Resource guard failure'
    finally:
        # Keep the receipt path reachable after launch, polling, or transition
        # failures. Cleanup is best effort and never replaces the original cause.
        if child is not None:
            try:
                if child.poll() is None:
                    subprocess.run(['systemctl', '--user', 'stop', unit], check=False)
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        subprocess.run(['systemctl', '--user', 'kill', '--kill-whom=all', unit], check=False)
                        child.kill()
                        child.wait(timeout=5)
                supervisor_status = child.wait(timeout=5)
            except Exception as exception:
                cleanup_errors.append(type(exception).__name__ + ': ' + str(exception))
                termination = termination or 'Resource guard cleanup failure'
        finished = time.monotonic()
        elapsed = finished - started
        if frozen_at is not None:
            frozen_seconds += max(0.0, finished - frozen_at)
        measured = {}
        try:
            measured = process_measurements(logpath.read_text())
        except Exception as exception:
            cleanup_errors.append('log accounting: ' + str(exception))
        process_cpu_available = all(k in measured for k in ('user_cpu_seconds', 'system_cpu_seconds'))
        cpu_observations = [service_cpu_seconds(sample['service']) for sample in samples]
        cpu_observations = [value for value in cpu_observations if value is not None]
        if process_cpu_available:
            cpu_observations.append(measured['user_cpu_seconds'] + measured['system_cpu_seconds'])
        # None is unmeasured; a real observed zero remains numeric zero.
        final_cpu = max(cpu_observations) if cpu_observations else None
        final_disk = None
        try:
            final_disk = output_bytes(args.output)
        except Exception as exception:
            cleanup_errors.append('output accounting: ' + str(exception))
            termination = termination or 'Output accounting failure'
        if termination is None:
            if final_cpu is not None and final_cpu > args.cpu_seconds:
                termination = 'CPU budget exhausted'
            elif final_disk is not None and final_disk > args.disk_mib * 1024**2:
                termination = 'Disk budget exhausted'
            elif elapsed > args.wall_seconds:
                termination = 'Wall budget exhausted'
        status = (supervisor_status if supervisor_status not in (None, 0) else
                  0 if supervisor_status == 0 and termination is None else 1)
        record = {'unit': unit, 'settings': settings, 'memory_profile': args.memory_profile, 'environment': env,
                  'command': command, 'exit_code': status, 'supervisor_exit_code': supervisor_status,
                  'elapsed': elapsed, 'samples': samples,
                  'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'budget': dict(cpu_seconds=args.cpu_seconds, wall_seconds=args.wall_seconds, disk_mib=args.disk_mib),
                  'termination': termination, 'log': str(logpath), 'process_measurements': measured,
                  'final_accounting': dict(CPU_seconds=final_cpu, output_bytes=final_disk, elapsed_seconds=elapsed,
                      process_CPU_timing_available=process_cpu_available, CPU_accounting_available=final_cpu is not None),
                  'frozen_seconds': frozen_seconds, 'guard_error': guard_error, 'cleanup_errors': cleanup_errors}
        with report.open('x') as receipt:
            receipt.write(json.dumps(record, indent=2) + '\n')
        print(json.dumps({k: record[k] for k in ['unit', 'exit_code', 'elapsed', 'termination', 'log']}), flush=True)
    raise SystemExit(status)


if __name__ == '__main__':
    main()
