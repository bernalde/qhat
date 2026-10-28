#!/usr/bin/env python3
"""Fixed-canonical serialization/replay parity probe on a retained preflight.

Run under the root's existing resource guard, only after the full preflight is
complete and the byte-equivalent stream patch is committed. No old artifacts
are written, no constructor is called, and no independent series is repeated.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMBA_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.dont_write_bytecode = True
PROJECT_ROOT = Path('/tmp/qhat-scalable-ordering-2026-10-27/worktree')
ARMS = ('canonical', 'magnitude', 'current', 'first', 'pauli', 'cubic')
FIXED_ARM = 'canonical'
FIXED_ORIGINAL_CASE = 'CH4_r1p09'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, label, *, positive=False):
    require(type(value) in (int, float) and math.isfinite(value) and
            (value > 0 if positive else value >= 0), 'Invalid finite value: ' + label)
    return value


def nonnegative_integer(value, label, *, maximum=None):
    require(type(value) is int and value >= 0 and (maximum is None or value <= maximum),
            'Invalid count: ' + label)
    return value


def artifact(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'Expected regular artifact: ' + str(path))
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=value.hexdigest())


def ratio(numerator, denominator):
    if denominator == 0:
        return dict(value=None, reason='both_zero' if numerator == 0 else 'zero_old_measurement')
    value = numerator / denominator
    return dict(value=value if math.isfinite(value) else None,
                reason=None if math.isfinite(value) else 'nonfinite_quotient')


class Probe:
    def __init__(self, output, entry_cpu, entry_wall):
        self.output, self.cpu, self.wall = output, entry_cpu, entry_wall
        self.active_cpu = self.active_wall = None
        self.data = dict(schema='qhat-fixed-canonical-stream-probe-v1', status='initializing',
            started_UTC=datetime.now(timezone.utc).isoformat(), fixed_arm=FIXED_ARM,
            fixed_original_case=FIXED_ORIGINAL_CASE, phases={}, consumed_artifacts={},
            probe_source=artifact(Path(__file__).resolve()), main_errors_exposed=False,
            new_training=False, independent_series_repeated=False,
            scope='One fixed retained altered-CH4 canonical circuit. Serialization CPU ratio is '
                'descriptive for this one pair; no method selection, general speed claim, new '
                'scientific benchmark, or formal review is implied.',
            outer_cost_scope='Whole outer CPU is this process plus waited child CPU measured before '
                'the final metadata write. It includes stdlib imports and current-module imports. '
                'Final journal writing, file closure, remaining object destruction/interpreter '
                'shutdown and the parent guard overhead are excluded; the external ledger is '
                'authoritative for the complete charge.')
        self.save()

    def save(self):
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.data['whole_outer_CPU_seconds_before_metadata_write'] = (
            time.process_time() + children.ru_utime + children.ru_stime)
        self.data['self_CPU_seconds_since_main_entry'] = time.process_time() - self.cpu
        self.data['wall_seconds_since_main_entry'] = time.monotonic() - self.wall
        path = self.output / 'PROBE.json'
        temporary = path.with_suffix('.json.partial')
        with temporary.open('w', encoding='utf-8') as stream:
            json.dump(self.data, stream, indent=2, allow_nan=False)
            stream.write('\n')
        temporary.replace(path)

    def begin(self, phase):
        require(self.active_cpu is None, 'Previous probe phase is still active')
        self.data.update(status='running', active_phase=phase)
        self.active_cpu, self.active_wall = time.process_time(), time.monotonic()
        self.active_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.save()

    def end(self, **details):
        child = resource.getrusage(resource.RUSAGE_CHILDREN)
        phase = self.data['active_phase']
        self.data['phases'][phase] = dict(CPU_seconds=time.process_time() - self.active_cpu,
            child_CPU_seconds=(child.ru_utime + child.ru_stime - self.active_children.ru_utime -
                               self.active_children.ru_stime),
            wall_seconds=time.monotonic() - self.active_wall, **details)
        self.active_cpu = self.active_wall = None
        self.data['active_phase'] = None
        self.save()

    def bind(self, key, path):
        value = artifact(path)
        previous = self.data['consumed_artifacts'].get(key)
        require(previous is None or previous == value, 'Previously bound artifact changed: ' + key)
        self.data['consumed_artifacts'][key] = value
        return value


def validate_resources(value, compiler):
    require(isinstance(value, dict), 'Missing old compiled resources')
    for key in ('gates', 'single_qubit_gates', 'two_qubit_gates',
                'unit_duration_ASAP_depth', 'CX_dependency_depth'):
        nonnegative_integer(value[key], 'old resources/' + key, maximum=compiler.MAX_STORED_GATES)
    require(value['gates'] == value['single_qubit_gates'] + value['two_qubit_gates'] and
            value['CX_dependency_depth'] <= value['unit_duration_ASAP_depth'] <= value['gates'],
            'Old gate/depth counts differ')


def run(args, entry_cpu, entry_wall):
    preflight = args.preflight.resolve(strict=True)
    project = args.project_root.resolve(strict=True)
    output = args.output.resolve()
    require(preflight.is_dir() and not output.exists(), 'Existing preflight and fresh output required')
    for root in (project, preflight):
        require(not output.is_relative_to(root) and not root.is_relative_to(output),
                'Probe output must be outside project and retained preflight')
    output.mkdir(parents=True, exist_ok=False)
    probe = Probe(output, entry_cpu, entry_wall)
    try:
        probe.begin('imports_and_authentication')
        source_dir = project / 'numerical_audit/scalable_ordering'
        sys.path.insert(0, str(source_dir))
        import common as k
        import compiler
        import compilation
        import reference as reference_module
        import numpy as np
        for module in (k, compiler, compilation, reference_module):
            require(Path(module.__file__).resolve().parent == source_dir, 'Wrong current-module import')
        frozen_path = preflight / 'PREFLIGHT_ORDERS_FROZEN.json'
        frozen = k.read(frozen_path)
        frozen_binding = probe.bind('old/PREFLIGHT_ORDERS_FROZEN.json', frozen_path)
        current_sources, current_environment = k.source_map(), k.environment()
        current_head = k.committed_sources()
        require(current_environment == frozen['environment'], 'Old/current environment differs')
        old_head = frozen['source_commit']
        require(isinstance(old_head, str) and len(old_head) == 40 and
                all(c in '0123456789abcdef' for c in old_head), 'Invalid frozen old head')
        require(current_head != old_head, 'Commit the stream optimization before probing')
        require(subprocess.run(['git', '-C', str(project), 'merge-base', '--is-ancestor', old_head, current_head],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0,
                'Old source head is not an ancestor of the current committed source')
        for relative, expected in frozen['sources'].items():
            path = Path(relative)
            require(not path.is_absolute() and '..' not in path.parts, 'Unsafe frozen source path')
            raw = subprocess.check_output(['git', '-C', str(project), 'show', old_head + ':' + relative])
            require(len(raw) <= 32 * 1024**2 and hashlib.sha256(raw).hexdigest() == expected,
                    'Frozen old committed-source SHA-256 differs: ' + relative)
        probe.data['source_identity'] = dict(old_commit=old_head, current_commit=current_head,
            frozen_old_source_sha256=frozen['sources'], current_source_sha256=current_sources,
            environment=current_environment, old_committed_blobs_authenticated=True,
            changed_sources={name: dict(old=frozen['sources'].get(name), current=current_sources.get(name))
                for name in sorted(set(frozen['sources']) | set(current_sources))
                if frozen['sources'].get(name) != current_sources.get(name)})
        require(set(frozen['orders']) == set(ARMS) and frozen['original_case'] == FIXED_ORIGINAL_CASE and
                frozen['scale'] == '7/5' and frozen['T'] == 1.5 and frozen['cap_multiple'] == 30,
                'Probe requires the fixed retained largest altered-CH4 preflight')
        complete_path = preflight / 'PREFLIGHT_COMPLETE.json'
        complete = k.read(complete_path)
        probe.bind('old/PREFLIGHT_COMPLETE.json', complete_path)
        require(complete['status'] == 'passed' and complete['original_case'] == FIXED_ORIGINAL_CASE and
                complete['main_cases_evaluated'] is False and complete['main_errors_exposed'] is False and
                len(complete['arms']) == 6 and {r['arm'] for r in complete['arms']} == set(ARMS) and
                all(r['status'] == 'complete' for r in complete['arms']), 'Full retained preflight has not passed')

        model_path = preflight / 'model.json'
        model = k.read(model_path)
        model_binding = probe.bind('old/model.json', model_path)
        require(model_binding['sha256'] == frozen['model_sha256'] == complete['model_sha256'] and
                model['case'] == 'preflight_' + FIXED_ORIGINAL_CASE + '_scaled7over5' and
                model['n'] == 12 and model['electrons'] == 6 and model['preflight']['scale'] == '7/5',
                'Old altered-model hash/identity differs')
        k.blocks(model)
        order_path = preflight / 'order-canonical.json'
        order = k.read(order_path)
        order_binding = probe.bind('old/order-canonical.json', order_path)
        require(order_binding['sha256'] == frozen['orders'][FIXED_ARM] and order['arm'] == FIXED_ARM and
                order['model_sha256'] == frozen['model_sha256'] and order['partition'] == 'conserving' and
                order['order'] == list(range(len(model['blocks']))), 'Frozen canonical order differs')
        compilation.validate_record(model, order)
        old_path = preflight / 'result-canonical.json'
        old = k.read(old_path)
        probe.bind('old/result-canonical.json', old_path)
        require(next(row for row in complete['arms'] if row['arm'] == FIXED_ARM) == old,
                'Old result differs from the completed preflight row')
        steps = compilation.affordable_steps(order['unit_lengths'], 30)
        old_identity = dict(case=model['case'], arm=FIXED_ARM, mode='preflight', T=1.5,
            steps=steps, cap_multiple=30, formula=2, preparation='HF_prefix', compilation_status='complete')
        require(old['arm'] == FIXED_ARM and old['status'] == 'complete' and old['identity'] == old_identity and
                old['compiler'] == compiler.VERSION and old['steps'] == steps and
                old['ordering_CPU_seconds'] == order['CPU_seconds'] and
                old['ordering_status'] == order['status'] and old['ordering_evaluations'] == order['evaluations'],
                'Old canonical result identity/ordering metadata differs')
        validate_resources(old['compiled_resources'], compiler)
        require(old['partial_compiled_resources'] is None and
                old['ideal_rotations'] == old['expected_ideal_rotations'] ==
                    compilation.rotation_count(order['unit_lengths'], steps) <= old['rotation_cap'] and
                old['rotation_cap'] == 30 * len(model['coefficients']), 'Old rotation/gate bounds differ')
        nonnegative_integer(old['generated'], 'old generated gates', maximum=compiler.MAX_GENERATED_GATES)
        old_serialization = finite(old['serialization']['CPU_seconds'], 'old serialization CPU')
        finite(old['serialization']['wall_seconds'], 'old serialization wall')
        for key in ('error', 'squared_infidelity', 'number_leakage', 'joint_sector_leakage'):
            finite(old[key], 'old ' + key)
        finite(old['norm'], 'old norm', positive=True)
        old_trace = preflight / 'canonical.ndjson.gz'
        trace_binding = probe.bind('old/canonical.ndjson.gz', old_trace)
        require(trace_binding['sha256'] == old['trace_sha256'] and trace_binding['bytes'] == old['trace_bytes'],
                'Old trace SHA-256/byte count differs')
        reference_path = preflight / 'reference.npz'
        reference_binding = probe.bind('old/reference.npz', reference_path)

        # A continuation journals the previously authenticated reference, so
        # retain that stronger prior binding when the artifact is available.
        if complete.get('continuation'):
            journal_path = Path(complete['continuation_record']).resolve(strict=True)
            continuation = k.read(journal_path)
            probe.bind('old/continuation_record.json', journal_path)
            require(continuation['status'] == 'passed' and
                    continuation['original_frozen_sha256'] == frozen_binding['sha256'] and
                    continuation['same_frozen_orders'] is True and continuation['new_training'] is False and
                    continuation['main_errors_exposed'] is False and continuation['complete_arms'] == 6 and
                    continuation['reference_check']['status'] == 'passed' and
                    continuation['reused_artifacts']['reference.npz']['sha256'] == reference_binding['sha256'],
                    'Retained continuation/reference authentication differs')
        probe.end(gates=old['compiled_resources']['gates'], rotations=old['ideal_rotations'],
                  fixed_case=FIXED_ORIGINAL_CASE, fixed_arm=FIXED_ARM)

        probe.begin('load_old_gate_list')
        reader = compiler.read_trace(old_trace)
        gate = None
        gates = []
        try:
            header = next(reader)
            require(header['n'] == model['n'] and header['compiler'] == compiler.VERSION and
                    header['identity'] == old_identity, 'Loaded old trace header differs')
            for gate in reader:
                require(len(gates) < compiler.MAX_STORED_GATES, 'Stored-gate ceiling while loading old trace')
                gates.append(gate)
        finally:
            reader.close()
        gate_count = len(gates)
        require(gate_count == old['compiled_resources']['gates'], 'Loaded old gate inventory differs')
        probe.end(gates=gate_count, note='Load/parse time is separate from serialization timing.')

        probe.begin('write_new_trace')
        new_trace = output / 'canonical-stream-probe.ndjson.gz'
        measurements = {}
        new_digest = compiler.write_trace(new_trace, model['n'], gates, old_identity, measurements=measurements)
        require(new_digest == old['trace_sha256'], 'New compressed trace is not byte-identical to old trace')
        new_binding = probe.bind('new/canonical-stream-probe.ndjson.gz', new_trace)
        require(new_binding['sha256'] == trace_binding['sha256'] and new_binding['bytes'] == trace_binding['bytes'],
                'New trace hash/size authentication differs')
        new_serialization = finite(measurements['serialization']['CPU_seconds'], 'new serialization CPU')
        probe.data['serialization_comparison'] = dict(old=old['serialization'],
            new=measurements['serialization'], new_hashing=measurements['hashing'],
            new_over_old_CPU_ratio=ratio(new_serialization, old_serialization),
            compressed_bytes_identical=True,
            scope='One fixed circuit, old retained measurement versus this fresh measurement. '
                'Gate-list loading, independent replay and reference work are excluded from this ratio.')
        probe.end(gates=gate_count, serialization=measurements['serialization'], hashing=measurements['hashing'])

        probe.begin('release_gate_storage')
        del gates
        del gate
        del reader
        collected = gc.collect()
        probe.end(gc_objects_collected=collected,
                  note='All Gate-list and generator references deleted before allocating Reference.')

        probe.begin('reference_allocation')
        context = reference_module.Reference(model)
        probe.end(spin_modes=model['n'], dimension=context.dimension)

        probe.begin('audit_new_trace')
        actual, counts = context.audit_trace(new_trace, old_identity, compiler.VERSION)
        actual_rotations = counts.pop('actual_Rz_rotations')
        require(counts == old['compiled_resources'], 'New primitive DAG resources differ from retained old row')
        probe.end(gates=counts['gates'], actual_Rz_rotations=actual_rotations,
                  DAG_resources=counts, old_resources_equal=True)

        probe.begin('ideal_state_parity')
        ideal = context.ideal_formula(order, old['steps'], 1.5)
        overlap = np.vdot(ideal, actual)
        require(abs(overlap) > 0, 'New actual/ideal overlap is zero')
        max_difference = float(np.max(abs(actual - ideal * overlap / abs(overlap))))
        require(math.isfinite(max_difference) and max_difference < 2e-9, 'New actual/ideal state parity failed')
        probe.end(max_state_difference=max_difference, tolerance=2e-9)

        probe.begin('saved_reference_metrics_parity')
        exact = k.load_vector(reference_path, 'exact', model['n'])
        values = reference_module.metrics(exact, actual, model)
        differences = {key: abs(values[key] - old[key]) for key in ('error', 'squared_infidelity')}
        require(all(math.isfinite(value) and value <= 1e-12 for value in differences.values()),
                'New/old error or squared-infidelity parity failed')
        require(abs(values['norm'] - 1) < 1e-8, 'New state norm drift exceeds fixed tolerance')
        require(values['joint_sector_leakage'] < 1e-16 and values['number_leakage'] < 1e-16,
                'New conserving state sector leakage exceeds fixed tolerance')
        probe.data['metric_parity'] = dict(old={key: old[key] for key in values}, new=values,
            absolute_differences=differences, error_parity_tolerance=1e-12,
            norm_drift_tolerance=1e-8, sector_leakage_tolerance=1e-16,
            reference_source='Saved bounded NPZ; independent series already passed in completed preflight.')
        probe.end(error_and_squared_infidelity_parity=True, norm_and_sector_checks=True,
                  independent_series_repeated=False)

        probe.begin('final_identity_check')
        require(k.source_map() == current_sources and k.environment() == current_environment,
                'Current execution source/environment changed during probe')
        require(subprocess.check_output(['git', '-C', str(project), 'rev-parse', 'HEAD'], text=True).strip() == current_head,
                'Current committed head changed during probe')
        for key, path in (
            ('old/PREFLIGHT_ORDERS_FROZEN.json', frozen_path), ('old/PREFLIGHT_COMPLETE.json', complete_path),
            ('old/model.json', model_path), ('old/order-canonical.json', order_path),
            ('old/result-canonical.json', old_path), ('old/canonical.ndjson.gz', old_trace),
            ('old/reference.npz', reference_path)):
            probe.bind(key, path)
        probe.end(old_artifacts_unchanged=True)
        probe.data.update(status='passed', active_phase=None, gates=gate_count,
                          actual_Rz_rotations=actual_rotations, main_errors_exposed=False)
        probe.save()
        print(json.dumps(dict(status='passed', gates=gate_count, actual_Rz_rotations=actual_rotations,
            old_serialization_CPU_seconds=old_serialization, new_serialization_CPU_seconds=new_serialization,
            new_over_old_serialization_CPU_ratio=probe.data['serialization_comparison']['new_over_old_CPU_ratio'],
            whole_outer_CPU_seconds_before_metadata_write=probe.data['whole_outer_CPU_seconds_before_metadata_write']),
            sort_keys=True), flush=True)
        return probe.data
    except BaseException as error:
        probe.data.update(status='failed', exception=dict(type=type(error).__name__, message=str(error)))
        if probe.active_cpu is not None:
            child = resource.getrusage(resource.RUSAGE_CHILDREN)
            probe.data['interrupted_phase_charge'] = dict(phase=probe.data['active_phase'],
                CPU_seconds=time.process_time() - probe.active_cpu,
                child_CPU_seconds=(child.ru_utime + child.ru_stime - probe.active_children.ru_utime -
                                   probe.active_children.ru_stime),
                wall_seconds=time.monotonic() - probe.active_wall)
        probe.save()
        raise


def main():
    cpu, wall = time.process_time(), time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Fresh directory outside Q and old artifacts')
    parser.add_argument('--project-root', type=Path, default=PROJECT_ROOT)
    run(parser.parse_args(), cpu, wall)


if __name__ == '__main__':
    main()
