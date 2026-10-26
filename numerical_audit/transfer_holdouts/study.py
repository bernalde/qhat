"""Freeze the complete transfer construction panel before measuring its errors.

Preparation requires a coordinator-attested reviewed practical interface. The
attestation binds source bytes; it does not automate reviewer judgment. The
separate final publication phase must replace a development pin by its verified
actual numerical merge.
"""
from collections import defaultdict
import hashlib
import io
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile
import numpy as np
from scipy.sparse.linalg import expm_multiply
import transfer_common as k
import compilation
import cubic_response
import panel
import state_preparation
import training

for _module in (k, compilation, cubic_response, panel, state_preparation, training):
    k.require(Path(_module.__file__).resolve() == Path(__file__).resolve().parent /
              (_module.__name__ + '.py'), 'Wrong transfer source binding: ' + _module.__name__)

MAX_FILES = 16000
MAX_BYTES = 512 * 1024**2
PHASE_LIMITS = {
    'tests': dict(CPU_seconds=120, wall_seconds=900, disk_MiB=64),
    'full_preflight': dict(CPU_seconds=600, wall_seconds=3600, disk_MiB=256),
    'measurement': dict(CPU_seconds=2400, wall_seconds=14400, disk_MiB=512),
    'independent': dict(CPU_seconds=2400, wall_seconds=14400, disk_MiB=1024),
    'reporting': dict(CPU_seconds=300, wall_seconds=1800, disk_MiB=1024),
}


def read(path):
    path = Path(path)
    k.require(path.stat().st_size <= 8 * 1024**2, 'JSON byte bound before parsing')
    return k.read(path)


def inventory(root, max_files=MAX_FILES, max_bytes=MAX_BYTES):
    paths, total, entries = {}, 0, 0
    for path in Path(root).rglob('*'):
        entries += 1
        k.require(entries <= max_files + 1024, 'Evidence directory-entry ceiling')
        k.require(not path.is_symlink(), 'Evidence symlink refused')
        if path.is_file():
            size = path.stat().st_size
            k.require(size <= 32 * 1024**2, 'Individual evidence file bound before hashing')
            total += size
            k.require(len(paths) < max_files and total <= max_bytes, 'Evidence file/byte ceiling')
            paths[str(path.relative_to(root))] = k.sha(path)
    return {name: paths[name] for name in sorted(paths)}


def practical_interface(path):
    record = read(path)
    k.require(record['schema'] == 1 and record['repository'] == 'bernalde/qhat' and record['issue'] == 16 and
              record['state'] == 'reviewed_frozen_interface' and
              record['automated_handoff'] == 'READY FOR HUMAN REVIEW',
              'Reviewed practical interface required before transfer execution')
    k.require(record['requested_reviewers'] == ['gpt-6-astra/high', 'gemini-4-argon-high/high'] and
              record['actionable_findings_remaining'] == 0 and record['controller_policy_version'] == 2,
              'Required paired-review convergence is not attested')
    k.require(record['cost_convention'] == 'caller_through_return_including_serialization_hash_cleanup',
              'Practical accounting convention differs')
    k.require(isinstance(record['head'], str) and re.fullmatch('[0-9a-f]{40}', record['head']) is not None,
              'Exact reviewed commit SHA required')
    sources = record['source_sha256']
    k.require(all(isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest) is not None
                  for digest in sources.values()), 'Exact source SHA256 required')
    required = {'numerical_audit/practical_accuracy/' + name for name in
                ('common.py', 'construction.py', 'study.py', 'verify.py', 'PROTOCOL.md')}
    required.add('numerical_audit/compiled_accuracy/compiler.py')
    k.require(required <= set(sources), 'Practical interface source inventory incomplete')
    for name, digest in sources.items():
        k.safe(name)
        raw = subprocess.check_output(['git', '-C', str(k.ROOT), 'show', record['head'] + ':' + name])
        k.require(hashlib.sha256(raw).hexdigest() == digest, 'Practical reviewed Git source differs')
    k.require(k.sha(k.AUDIT / 'compiled_accuracy/compiler.py') ==
              sources['numerical_audit/compiled_accuracy/compiler.py'], 'Compiler differs from reviewed interface')
    return record


class References:
    """Read already independently checked references, with bounded archive parsing."""
    def __init__(self, root, preflight):
        self.root = Path(root)
        actual = inventory(self.root, max_files=1024, max_bytes=128 * 1024**2)
        self.protocol = read(self.root / 'PROTOCOL_FROZEN.json')
        frozen = read(self.root / 'REFERENCE_STATES_FROZEN.json')
        verified = read(self.root / 'REFERENCE_VERIFIED.json')
        k.require(set(actual) == set(frozen['artifacts']) | {'REFERENCE_STATES_FROZEN.json',
                  'REFERENCE_VERIFIED.json', 'INDEPENDENT_REFERENCE_LABELS.json'}, 'Unexpected reference evidence inventory')
        k.require(verified['status'] == 'passed' and
                  verified['reference_freeze_sha256'] == k.sha(self.root / 'REFERENCE_STATES_FROZEN.json') and
                  verified['producer_sources'] == self.protocol['sources'] and
                  verified['independent_labels_sha256'] == k.sha(self.root / 'INDEPENDENT_REFERENCE_LABELS.json'),
                  'Independent reference audit binding differs')
        k.require(self.protocol['preflight'] == preflight and self.protocol['states'] == list(k.STATES) and
                  self.protocol['validation_times'] == list(k.TIMES) and
                  self.protocol['training_time'] == k.TRAIN_TIME, 'Frozen reference panel differs')
        for name, digest in frozen['artifacts'].items():
            k.require(k.sha(self.root / k.safe(name)) == digest, 'Frozen reference bytes changed')
        registry = read(k.HERE / 'REGISTRY.json')
        specs = [registry['chemistry_preflight']] if preflight else registry['cases']
        self.cases = [row['case'] for row in specs]
        self.inputs, _ = state_preparation.authenticate_inputs(self.root / 'inputs_evidence', preflight)
        k.require({row['case'] for row in verified['models']} == set(self.cases) and
                  verified['states'] == 5 * len(self.cases) == frozen['states'] and
                  frozen['models'] == len(self.cases), 'Complete independently checked reference cohort required')
        self.records, self.vectors = {}, {}
        for case in self.cases:
            record = read(self.root / 'references' / (case + '.json'))
            k.require(record['case'] == case and set(record['states']) == set(k.STATES) and
                      record['model_sha256'] == self.inputs['successes'][case] and
                      record['vectors_sha256'] == k.sha(self.root / k.safe(record['vectors_file'])),
                      'Reference state/model identity differs')
            for state in k.STATES:
                cost = record['states'][state]['CPU_seconds']
                k.require(math.isfinite(cost) and cost >= 0, 'Invalid reference preparation charge')
            self.records[case] = record

    def model_path(self, case):
        k.require(case in self.cases, 'Undeclared reference case')
        return self.root / 'inputs_evidence/inputs' / case / 'model.json'

    def vector(self, case, state):
        k.require(case in self.cases and state in k.STATES, 'Undeclared reference state')
        if case not in self.vectors:
            model = read(self.model_path(case))
            k.require(type(model['n']) is int and 1 <= model['n'] <= 6, 'Register bound before allocation')
            path = self.root / k.safe(self.records[case]['vectors_file'])
            with path.open('rb') as stream:
                packed = stream.read(32769)
            k.require(len(packed) <= 32768 and hashlib.sha256(packed).hexdigest() ==
                      self.records[case]['vectors_sha256'], 'Reference archive bytes changed before loading')
            with zipfile.ZipFile(io.BytesIO(packed)) as archive:
                entries = archive.infolist()
                k.require(len(entries) == 1 and entries[0].filename == 'vectors.npy' and
                          entries[0].file_size <= 8192, 'Reference archive shape/byte container guard')
                raw = archive.read(entries[0])
            stream = io.BytesIO(raw)
            version = np.lib.format.read_magic(stream)
            k.require(version in ((1, 0), (2, 0)), 'Supported bounded NPY header required')
            reader = (np.lib.format.read_array_header_1_0 if version == (1, 0) else
                      np.lib.format.read_array_header_2_0)
            shape, _, dtype = reader(stream)
            k.require(shape == (1 << model['n'], 5) and dtype == np.dtype('complex128') and
                      len(raw) - stream.tell() == (1 << model['n']) * 5 * dtype.itemsize,
                      'Reference array shape/dtype before allocation')
            values = np.load(io.BytesIO(raw), allow_pickle=False)
            k.require(np.all(np.isfinite(values)), 'Nonfinite reference archive')
            self.vectors[case] = values
        return self.vectors[case][:, k.STATES.index(state)].copy()

    def cost(self, case, state):
        return self.records[case]['states'][state]['CPU_seconds']


def constructor(reference, specification, *, warmup=False):
    case, state, method = (specification[name] for name in ('case', 'state', 'method'))
    model_path, digest = reference.model_path(case), reference.inputs['successes'][case]
    psi = reference.vector(case, state)
    charge = 0. if warmup else reference.cost(case, state)
    cpu, wall = time.process_time(), time.monotonic()
    if method == 'current':
        k.require(k.sha(model_path) == digest, 'Current-placement model changed')
        model = read(model_path)
        placed = training.current_construction(model)
        order = placed['order']
        result = dict(record=training.layout(model, order, 'conserving'), current_placement=placed, training_time=1.,
                      reference_preparation_CPU_charge=charge, phases={}, scope='Legacy fixed-T1 HF placement')
        del model, order
    elif method == 'cubic_JD':
        result = training.cubic(model_path, digest, psi, reference_CPU_charge=charge)
    else:
        k.require(method in ('first_JD', 'pauli_JD'), 'Unknown training method')
        result = training.first(model_path, digest, psi, 'pauli' if method == 'pauli_JD' else 'conserving',
                                reference_CPU_charge=charge)
    elapsed = time.process_time() - cpu
    result.update(caller_CPU_seconds=elapsed, caller_wall_seconds=time.monotonic() - wall,
                  standalone_reference_and_ordering_CPU_charge=charge + elapsed,
                  timing='Outer caller through complete constructor return;reference archive IO is separate',
                  constructor_identity=specification, warmup=warmup)
    return result


def prepare(reference_root, output, interface_path, *, preflight=False):
    cpu, wall = time.process_time(), time.monotonic()
    output = Path(output)
    k.require(not output.exists(), 'Fresh transfer output required;retain earlier attempts')
    dependency = practical_interface(interface_path)
    commit = k.committed_sources()
    reference = References(reference_root, preflight)
    plan = panel.build_plan(reference.cases, preflight=preflight)
    output.mkdir(parents=True)
    shutil.copytree(reference.root, output / 'reference_evidence')
    reference = References(output / 'reference_evidence', preflight)
    for name in k.sources():
        target = output / 'source_snapshot' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(k.ROOT / name, target)
    k.write(output / 'PRACTICAL_INTERFACE.json', dependency)
    k.write(output / 'PLAN_FROZEN.json', plan)
    k.write(output / 'PROTOCOL_FROZEN.json', dict(source_commit=commit, sources=k.sources(), environment=k.environment(),
            preflight=preflight, plan_sha256=k.sha(output / 'PLAN_FROZEN.json'),
            reviewed_practical_interface_sha256=k.sha(output / 'PRACTICAL_INTERFACE.json'),
            reference_freeze_sha256=k.sha(reference.root / 'REFERENCE_STATES_FROZEN.json'),
            reference_verified_sha256=k.sha(reference.root / 'REFERENCE_VERIFIED.json'),
            phase_limits=PHASE_LIMITS, quantum_nonHF_preparation_circuit_cost='unknown',
            warmup_policy=dict(training_state='HF', training_methods=('current',) + panel.METHODS,
                compilation_time=k.TRAIN_TIME, compilation_cap_multiple=10,
                formula_paths={'current': [1], 'first_JD': [1, 2], 'pauli_JD': [2], 'cubic_JD': [2]},
                actual_cost_scope='campaign only;not a repeated standalone construction charge')))
    ordering_actual = warm_actual = warm_compilation_actual = compilation_actual = 0.
    for case in reference.cases:
        warm = {}
        for method in ('current',) + panel.METHODS:
            specification = dict(case=case, state='HF', method=method)
            warm[method] = constructor(reference, specification, warmup=True)
            warm_actual += warm[method]['caller_CPU_seconds']
            for formula in ((1, 2) if method == 'first_JD' else (1,) if method == 'current' else (2,)):
                path = output / 'warmup_traces' / (case + '.' + method + '.f' + str(formula) + '.ndjson.gz')
                path.parent.mkdir(exist_ok=True)
                began = time.process_time()
                compiled = compilation.construct_candidate(reference.model_path(case), reference.inputs['successes'][case],
                    dict(warm[method]['record'], formula=formula), k.TRAIN_TIME, 10, 'HF', path,
                    dict(purpose='compilation_warmup', case=case, method=method))
                cost = time.process_time() - began
                compiled['outer_compilation_CPU_seconds'] = cost
                compiled['trace_file'] = str(path.relative_to(output))
                k.require(compiled['status'] == 'complete', 'Full compilation warmup required')
                warm_compilation_actual += cost
                warm[method].setdefault('compiled', []).append(compiled)
        k.write(output / 'warmups' / (case + '.json'), warm)
    for specification in plan['training']:
        result = constructor(reference, specification)
        ordering_actual += result['caller_CPU_seconds']
        k.write(output / 'constructors' / (specification['constructor'] + '.json'), result)
        print('FROZE_TRANSFER_CONSTRUCTOR', specification['constructor'], flush=True)
    k.write(output / 'ORDERS_FROZEN.json', dict(constructors=plan['expected_constructors'],
            artifacts={str(path.relative_to(output)): k.sha(path) for path in sorted((output / 'constructors').glob('*.json'))},
            scope='All main orders precede every validation-time error'))
    for candidate in plan['circuits']:
        case, state = candidate['case'], candidate['state']
        path = output / 'constructors' / (candidate['constructor'] + '.json')
        trained = read(path)
        trace = output / 'traces' / (candidate['candidate'] + '.ndjson.gz')
        trace.parent.mkdir(exist_ok=True)
        began = time.process_time()
        result = compilation.construct_candidate(reference.model_path(case), reference.inputs['successes'][case],
            dict(trained['record'], formula=candidate['formula']), candidate['T'], candidate['cap_multiple'], state,
            trace, dict(candidate=candidate['candidate'], case=case, constructor=candidate['constructor'],
                        constructor_sha256=k.sha(path), model_sha256=reference.inputs['successes'][case],
                        reference_vectors_sha256=reference.records[case]['vectors_sha256']))
        elapsed = time.process_time() - began
        compilation_actual += elapsed
        charge_states = panel.reference_charge_states(candidate['training_state'], state)
        reference_charge = sum(reference.cost(case, label) for label in charge_states)
        input_charge = read(reference.model_path(case))['generation_CPU_seconds']
        result.update(**candidate, trace_file=str(trace.relative_to(output)),
                      outer_compilation_CPU_seconds=elapsed, ordering_CPU_charge=trained['caller_CPU_seconds'],
                      reference_charge_states=charge_states, standalone_reference_CPU_charge=reference_charge,
                      input_generation_CPU_charge=input_charge,
                      standalone_without_input_CPU_charge=reference_charge + trained['caller_CPU_seconds'] + elapsed,
                      standalone_including_input_CPU_charge=input_charge + reference_charge + trained['caller_CPU_seconds'] + elapsed)
        k.write(output / 'candidates' / (candidate['candidate'] + '.json'), result)
    k.write(output / 'PREPARE_COSTS.json', dict(CPU_seconds=time.process_time() - cpu,
            wall_seconds=time.monotonic() - wall, actual_outer_training_CPU_seconds=ordering_actual,
            actual_training_warmup_CPU_seconds=warm_actual, actual_compilation_warmup_CPU_seconds=warm_compilation_actual,
            actual_compilation_CPU_seconds=compilation_actual,
            original_reference_preparation=read(reference.root / 'PREPARATION_COSTS.json'),
            repeated_standalone_charges_are_campaign_work=False,
            scope='Actual interface/reference/source authentication,copies,warmups,unique constructors/circuits and metadata;imports/process startup in resource receipt'))
    k.write(output / 'SCHEDULES_FROZEN.json', dict(artifacts=inventory(output),
            unique_candidates=plan['expected_unique_circuits'], logical_candidates=plan['expected_logical_circuits'],
            policy='All complete/partial exact primitive circuits precede all new validation errors'))


def bind(output):
    output = Path(output)
    protocol, frozen = read(output / 'PROTOCOL_FROZEN.json'), read(output / 'SCHEDULES_FROZEN.json')
    for name, digest in frozen['artifacts'].items():
        k.require(k.sha(output / k.safe(name)) == digest, 'Frozen transfer construction changed')
    k.require(protocol['plan_sha256'] == k.sha(output / 'PLAN_FROZEN.json') and
              protocol['reviewed_practical_interface_sha256'] == k.sha(output / 'PRACTICAL_INTERFACE.json'),
              'Transfer protocol/plan/dependency differs')
    return protocol, read(output / 'PLAN_FROZEN.json'), frozen


def run(output):
    output = Path(output)
    protocol, plan, frozen = bind(output)
    k.require(protocol['sources'] == k.sources(), 'Producer sources changed after construction freeze')
    k.require(not (output / 'RUN_STARTED.json').exists(), 'Preserve partial evaluations;explicit continuation required')
    reference = References(output / 'reference_evidence', protocol['preflight'])
    k.write(output / 'RUN_STARTED.json', dict(schedule_freeze_sha256=k.sha(output / 'SCHEDULES_FROZEN.json')))
    cpu, wall = time.process_time(), time.monotonic()
    physical, complete, censored = {}, 0, 0
    for case in reference.cases:
        model = read(reference.model_path(case))
        H = state_preparation.hamiltonian(model)
        for state in k.STATES:
            psi = reference.vector(case, state)
            replay = compilation.Replay(model, psi)
            for ti, T in enumerate(plan['validation_times']):
                began = time.process_time()
                exact = expm_multiply((-1j * T) * H, psi)
                reference_cost = time.process_time() - began
                ref_path = output / 'exact_references' / (case + '.' + state + '.t' + str(ti) + '.npz')
                ref_path.parent.mkdir(exist_ok=True)
                np.savez_compressed(ref_path, vector=exact)
                for specification in plan['circuits']:
                    if (specification['case'], specification['state'], specification['time_index']) != (case, state, ti):
                        continue
                    name = specification['candidate']
                    candidate = read(output / 'candidates' / (name + '.json'))
                    row = dict(candidate)
                    if candidate['status'] == 'complete':
                        began = time.process_time()
                        actual = replay.replay(output / candidate['trace_file'], candidate['identity'])
                        replay_cost = time.process_time() - began
                        row.update(**compilation.metrics(exact, actual, model),
                                   reference_CPU_seconds=reference_cost, replay_CPU_seconds=replay_cost)
                        path = output / 'states' / (name + '.npz')
                        path.parent.mkdir(exist_ok=True)
                        np.savez_compressed(path, vector=actual)
                        row.update(state_file=str(path.relative_to(output)), state_sha256=k.sha(path),
                                   at_numerical_floor=row['error'] <= k.FLOOR)
                        complete += 1
                    else:
                        censored += 1
                    row['candidate_sha256'] = k.sha(output / 'candidates' / (name + '.json'))
                    k.write(output / 'physical' / (name + '.json'), row)
                    physical[name] = row
        print('SCREENED_TRANSFER_CASE', case, flush=True)
    k.require(complete + censored == plan['expected_unique_circuits'], 'Incomplete unique circuit panel')
    logical, groups = [], defaultdict(dict)
    for index, specification in enumerate(plan['logical']):
        row = dict(physical[specification['candidate']])
        row.update(**specification, correlation=reference.records[specification['case']]['states'][specification['state']]['metadata'])
        k.write(output / 'rows' / ('logical-' + str(index).zfill(4) + '.json'), row)
        logical.append(row)
        groups[row['case'], row['state'], row['time_index'], row['cap_multiple']][row['arm']] = row
    comparisons = []
    for (case, state, ti, cap), arms in groups.items():
        k.require(set(arms) == set(k.ARMS), 'Logical arm cohort incomplete')
        for pair in plan['pairs']:
            comparisons.append(dict(case=case, state=state, T=plan['validation_times'][ti], cap_multiple=cap,
                                    **pair, **panel.compare(arms[pair['a']], arms[pair['b']])))
    k.require(len(logical) == plan['expected_logical_circuits'] and
              len(comparisons) == plan['expected_comparisons'], 'Logical/pair panel incomplete')
    k.write(output / 'COMPARISONS.json', comparisons)
    k.write(output / 'RUN_COSTS.json', dict(CPU_seconds=time.process_time() - cpu,
            wall_seconds=time.monotonic() - wall, scope='Full-space references and primitive screening;offline evaluation,not training'))
    k.write(output / 'RUN_COMPLETE.json', dict(status='producer_complete_independent_audit_pending',
            unique_complete=complete, unique_censored=censored, logical_circuits=len(logical), comparisons=len(comparisons),
            schedule_freeze_sha256=k.sha(output / 'SCHEDULES_FROZEN.json')))
