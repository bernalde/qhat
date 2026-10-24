"""Run one owned systemd service at low priority with monitored memory headroom."""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import time


def available():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    raise RuntimeError('MemAvailable unavailable')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('--logs', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.label.replace('-', '').isalnum():
        raise ValueError('Invalid command or service label')
    args.logs.mkdir(parents=True, exist_ok=True)
    unit = 'fermionic-prospective-' + args.label + '-' + str(os.getpid())
    report = args.logs / (args.label + '.resources.json')
    if report.exists():
        raise FileExistsError(report)
    while available() < 5.5 * 1024**3:
        print('WAITING FOR MEMORY HEADROOM', available(), flush=True)
        time.sleep(15)
    settings = dict(CPUQuota='25%', CPUWeight='5', IOWeight='10', Nice='19',
                    MemoryHigh='448M', MemoryMax='512M', MemorySwapMax='0',
                    TasksMax='32', MemoryAccounting='yes', CPUAccounting='yes',
                    RuntimeMaxSec='3600')
    invocation = ['systemd-run', '--user', '--wait', '--pipe', '--unit', unit,
                  '--working-directory', str(Path.cwd())]
    for key, value in settings.items():
        invocation += ['--property', key + '=' + value]
    env = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMBA_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1', PYTHONUNBUFFERED='1',
               NUMBA_CACHE_DIR=str(args.logs / 'numba-cache'),
               MPLCONFIGDIR=str(args.logs / 'matplotlib-cache'))
    invocation += ['/usr/bin/env'] + [k + '=' + v for k, v in env.items()] + command
    logpath = args.logs / (args.label + '.log')
    started = time.monotonic()
    samples = []
    frozen = False
    with logpath.open('x') as log:
        child = subprocess.Popen(invocation, stdout=log, stderr=subprocess.STDOUT)
        try:
            while child.poll() is None:
                headroom = available()
                sample = {'elapsed': time.monotonic() - started,
                          'available_bytes': headroom, 'frozen': frozen}
                result = subprocess.run(['systemctl', '--user', 'show', unit,
                                         '-p', 'MemoryCurrent', '-p', 'MemoryPeak',
                                         '-p', 'CPUUsageNSec'], capture_output=True, text=True)
                sample['service'] = result.stdout.strip()
                samples.append(sample)
                if headroom < 5 * 1024**3 and not frozen:
                    subprocess.run(['systemctl', '--user', 'freeze', unit], check=True)
                    frozen = True
                    print('FROZE OWN SERVICE FOR HEADROOM', unit, flush=True)
                elif headroom >= 5.5 * 1024**3 and frozen:
                    subprocess.run(['systemctl', '--user', 'thaw', unit], check=True)
                    frozen = False
                    print('RESUMED OWN SERVICE', unit, flush=True)
                time.sleep(3)
        finally:
            if child.poll() is None:
                subprocess.run(['systemctl', '--user', 'stop', unit], check=False)
                child.wait(timeout=30)
        status = child.wait()
    record = {'unit': unit, 'settings': settings, 'environment': env,
              'command': command, 'exit_code': status,
              'elapsed': time.monotonic() - started, 'samples': samples,
              'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'log': str(logpath)}
    report.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: record[k] for k in ['unit', 'exit_code', 'elapsed', 'log']}), flush=True)
    raise SystemExit(status)


if __name__ == '__main__':
    main()
