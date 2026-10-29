"""Small producer-side joins between worker arguments and owning OS evidence.

Original absolute argument paths are compared lexically, so an archived run can
be relocated. The analysis consumer implements these checks independently.
"""
import math
from pathlib import PurePosixPath

SCRIPT = 'numerical_audit/cost_matched_continuation/construct.py'
VALUES = frozenset(('--root', '--case', '--repetition', '--arm', '--budget',
                    '--raw', '--cpu-cap', '--wall-cap', '--julia-executable'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def absolute(value):
    require(isinstance(value, str) and len(value) <= 4096 and '\\' not in value and '\x00' not in value,
            'Bounded POSIX path')
    path = PurePosixPath(value)
    require(path.is_absolute() and str(path) == value and '..' not in path.parts,
            'Canonical original absolute path')
    return path


def relative(value):
    require(isinstance(value, str), 'Relative worker path required')
    path = PurePosixPath(value)
    require(not path.is_absolute() and str(path) == value and
            all(part not in ('.', '..') for part in path.parts) and bool(path.parts), 'Canonical archived worker path')
    return path


def worker_arguments(plan, receipt):
    command = plan['command']
    require(isinstance(command, list) and 3 <= len(command) <= 32 and
            all(isinstance(value, str) and len(value) <= 4096 and '\x00' not in value for value in command),
            'Bounded fixed worker command')
    require(receipt['command'] == command, 'PLAN and guard command differ')
    cwd = absolute(plan['cwd'])
    script = PurePosixPath(command[1])
    script = script if script.is_absolute() else cwd/script
    require(script == cwd/SCRIPT and '..' not in script.parts, 'PLAN runs a different worker source')
    result = {}
    index = 2
    while index < len(command):
        key = command[index]
        require(key not in result, 'Duplicate worker argument')
        if key == '--worker':
            result[key] = True
            index += 1
        else:
            require(key in VALUES and index+1 < len(command), 'Unknown or incomplete worker argument')
            result[key] = command[index+1]
            index += 2
    require(result.get('--worker') is True and
            {'--root', '--case', '--repetition', '--arm', '--raw', '--cpu-cap', '--wall-cap'} <= result.keys(),
            'Required worker identity arguments')
    return result


def raw_cpu(final, receipt):
    if 'terminal_CPU_measurement_available' in final:
        known = final['terminal_CPU_measurement_available']
        require(type(known) is bool and 'raw_terminal_CPU_seconds' in final, 'Terminal CPU declaration')
        raw = final['raw_terminal_CPU_seconds']
        require((raw is None) == (not known), 'Unknown terminal CPU must remain null')
        return raw
    # Legacy guard CPU includes live polls. It is a lower-bound diagnostic
    # until process timing or a successful post-exit service query is known.
    def valid(value):
        try:
            return type(value) in (int, float) and math.isfinite(value) and value >= 0
        except OverflowError:
            return False
    accounting = receipt.get('final_accounting', {})
    process = receipt.get('process_measurements', {})
    measured = []
    if accounting.get('process_CPU_timing_available') is True:
        user, system = process.get('user_cpu_seconds'), process.get('system_cpu_seconds')
        if valid(user) and valid(system) and valid(user+system):
            measured.append(user+system)
    service = final.get('final_service_snapshot', {}).get('CPUUsageNSec')
    if final.get('final_service_query_exit') == 0 and isinstance(service, str) and service.isascii() and service.isdigit():
        try:
            value = int(service)/1e9
        except (ValueError, OverflowError):
            value = None
        if valid(value):
            measured.append(value)
    if not measured:
        return None
    guard = accounting.get('CPU_seconds')
    require(guard is None or valid(guard), 'Finite nonnegative legacy guard CPU')
    return max(measured+([] if guard is None else [guard]))


def bind_worker(plan, receipt, owner, parent, worker, config):
    args = worker_arguments(plan, receipt)
    require(plan['phase'] == 'constructors', 'Owning job is not a constructor')
    require((args['--case'], int(args['--repetition']), args['--arm'], int(args.get('--budget', '0'))) ==
            (owner['case'], owner['repetition'], owner['arm'], owner['total_CPU_budget']),
            'PLAN constructor identity differs from selected worker')
    original_root = absolute(args['--root'])
    require(absolute(args['--raw']) == original_root/relative(owner['worker_output_path']),
            'PLAN writes a different worker output')
    parent_arm = owner['arm'] == 'parent'
    if not parent_arm:
        require(isinstance(parent, dict) and
                (parent.get('case'), parent.get('repetition'), parent.get('arm')) ==
                (owner['case'], owner['repetition'], 'parent'), 'Frozen parent identity differs from worker')
    expected_cpu = config['parent']['CPU_seconds'] if parent_arm else (
        owner['total_CPU_budget']-parent['costs']['total']['charged_CPU_seconds'])
    expected_wall = config['parent']['wall_seconds'] if parent_arm else math.floor(
        config['accounting']['total_wall_multiplier']*owner['total_CPU_budget']-
        parent['costs']['total']['wall_seconds'])
    cpu, wall = float(args['--cpu-cap']), int(args['--wall-cap'])
    require(math.isfinite(cpu) and 0 < cpu <= expected_cpu and expected_cpu-cpu <= 1e-12 and wall == expected_wall,
            'PLAN constructor allowance differs from frozen parent remainder')
    require(cpu == plan['attempt_CPU_seconds'] == receipt['budget']['cpu_seconds'] and
            wall == plan['attempt_wall_seconds'] == receipt['budget']['wall_seconds'], 'PLAN/guard allowance differs')
    if worker is not None:
        require(worker['guarded_CPU_cap'] == cpu and worker['guarded_wall_cap'] == wall,
                'Worker and owning PLAN allowance differ')
    return dict(original_root=str(original_root), original_worker_path=args['--raw'],
                CPU_allowance=cpu, wall_allowance=wall, verified=True)
