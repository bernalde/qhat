"""State-aware training adapters; no validation-time propagation is available.

This is preparation for the reviewed practical interface. Real panel search and
final publication remain gated on issue16's frozen review/merge dependency.
"""
import transfer_common as k
import numpy as np
import importlib
from pathlib import Path
import sys
import time
from cubic_response import Context as CubicContext, prior, validate

_saved = list(sys.path)
sys.path.insert(0, str(k.AUDIT / 'prospective_reference'))
for _name in ('engine', 'one_pass'):
    _module = importlib.import_module(_name)
    k.require(Path(_module.__file__).resolve() == k.AUDIT / 'prospective_reference' / (_name + '.py'),
              'Wrong current-placement source')
placement = k.load('transfer_original_placement', k.AUDIT / 'prospective_reference/campaign.py')
sys.path[:] = _saved


class FirstOrderContext:
    def __init__(self, model, psi, partition, T=k.TRAIN_TIME):
        blocks = validate(model)
        k.require(partition in ('conserving', 'pauli') and np.isfinite(T) and T > 0, 'First-order context domain')
        psi = np.asarray(psi, complex)
        k.require(psi.shape == (1 << model['n'],) and np.all(np.isfinite(psi)) and
                  abs(np.vdot(psi, psi) - 1) < 1e-10 and len(np.flatnonzero(psi != 0)) <= 9,
                  'Declared normalized state/support domain')
        self.blocks = blocks if partition == 'conserving' else [
            {w: c} for w, c in sorted(prior.decode(model['coefficients']).items())]
        # All six frozen targets contain33--77 terms. Preserve the inherited
        # feature guard rather than changing historical scientific source.
        k.require(len(self.blocks) <= 128, 'Inherited first-order feature domain exceeded')
        D = {w: c for b in blocks if k.engine.descriptor(b)[0] == 0 for w, c in b.items()}
        self.features = k.response.Features(self.blocks, D, psi, T, int(model['denominator']), model['n'])
        self.m = len(self.blocks)
        self.T = T
        self.words = len(model['coefficients'])
        self.lengths = [len(b) for b in self.blocks]

    def score(self, order, T, rotation_cap, deadline=float('inf'), clock=time.process_time):
        k.require(T == self.T and type(rotation_cap) is int and rotation_cap > 0, 'Training time/cap changed')
        if clock() >= deadline:
            raise prior.Interrupted(0)
        prior.permutation(order, self.m)
        value = self.features.score(self.features.amplitude(order))
        steps = rotation_cap // self.words
        k.require(steps >= 1 and np.isfinite(value) and value >= 0, 'Invalid first-order score/affordability')
        return dict(order=list(order), JD1=value, score=value / steps**2, steps=steps)


def current_construction(model):
    blocks = validate(model)
    diagonal = sorted((i for i, b in enumerate(blocks) if k.engine.descriptor(b)[0] == 0),
                      key=lambda i: (-abs(next(iter(blocks[i].values()))), next(iter(blocks[i]))))
    off = sorted((i for i, b in enumerate(blocks) if k.engine.descriptor(b)[0]),
                 key=lambda i: (-sum(abs(c) for c in blocks[i].values()), k.engine.descriptor(blocks[i])[0]))
    result = placement.current_placement(blocks, diagonal, off, (1 << model['electrons']) - 1,
                                         int(model['denominator']))
    result.update(diagonal_priority=diagonal, off_priority=off)
    return result


def current_order(model):
    return current_construction(model)['order']


def layout(model, order, partition):
    blocks = validate(model)
    if partition == 'pauli':
        blocks = [{w: c} for w, c in sorted(prior.decode(model['coefficients']).items())]
    prior.permutation(order, len(blocks))
    return dict(order=list(order), partition=partition,
                unit_lengths=[len(blocks[i]) for i in order],
                schedule=[[[[q, p] for q, p in w], str(c)]
                          for i in order for w, c in sorted(blocks[i].items())])


def _first(model_path, digest, psi, partition):
    phases = {}
    cpu = time.process_time()
    k.require(k.sha(model_path) == digest, 'Training model changed')
    model = k.read(model_path)
    validate(model)
    placed = current_construction(model)
    current = placed['order']
    phases['input_and_current_placement_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    context = FirstOrderContext(model, psi, partition)
    phases['context_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    if partition == 'pauli':
        blocks = validate(model)
        index = {w: i for i, w in enumerate(sorted(prior.decode(model['coefficients'])))}
        start = [index[w] for i in current for w in sorted(blocks[i])]
        del blocks
    else:
        start = current
    starts = [start, list(reversed(start)), list(range(context.m))]
    phases['initializer_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    cap = 30 * len(model['coefficients'])
    seeds = [context.score(order, k.TRAIN_TIME, cap) for order in starts]
    phases['seed_scores_CPU_seconds'] = time.process_time() - cpu
    allowance = 8. if partition == 'pauli' else 2.
    evaluations = min(325120, 20 * context.m * (context.m - 1)) if partition == 'pauli' else 3000
    moves = 10 if partition == 'pauli' else 5
    cpu = time.process_time()
    result = prior.insertion(context, seeds, k.TRAIN_TIME, cap, allowance, evaluations, moves)
    phases['search_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    record = layout(model, result['record']['order'], partition)
    del context, model
    phases['layout_and_cleanup_CPU_seconds'] = time.process_time() - cpu
    return dict(record=record, seeds=seeds, result=result, phases=phases,
                training_time=k.TRAIN_TIME, legacy_placement_initializer_time=1.,
                nominal_search_CPU_seconds=allowance, search_evaluation_cap=evaluations,
                max_moves=moves, current_order=current, current_placement=placed)


def first(model_path, digest, psi, partition, reference_CPU_charge=0.):
    k.require(np.isfinite(reference_CPU_charge) and reference_CPU_charge >= 0, 'Invalid reference charge')
    cpu, wall = time.process_time(), time.monotonic()
    result = _first(model_path, digest, psi, partition)
    elapsed = time.process_time() - cpu
    result.update(actual_local_CPU_seconds=elapsed, wall_seconds=time.monotonic() - wall,
                  accounting_overhead_CPU_seconds=elapsed - sum(result['phases'].values()),
                  reference_preparation_CPU_charge=reference_CPU_charge,
                  standalone_reference_and_ordering_CPU_charge=reference_CPU_charge + elapsed,
                  scope='Caller-through-return ordering charge;one-circuitcompilationandSCF/archiveIOreportedseparately')
    return result


def _cubic(model_path, digest, psi):
    phases = {}
    cpu = time.process_time()
    # This complete warm constructor is inside the outer cubic caller timer.
    # Charge reference preparation once at the outer boundary.
    warm = first(model_path, digest, psi, 'conserving', reference_CPU_charge=0.)
    phases['charged_first_order_warm_constructor_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    k.require(k.sha(model_path) == digest, 'Cubic training model changed')
    model = k.read(model_path)
    context = CubicContext(model, psi)
    phases['model_and_cubic_context_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    start = warm['record']['order']
    starts = [start, list(reversed(start)), warm['current_order']]
    cap = 30 * len(model['coefficients'])
    seeds = [context.score(order, k.TRAIN_TIME, cap) for order in starts]
    phases['cubic_seeds_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    result = prior.insertion(context, seeds, k.TRAIN_TIME, cap, 2., 3000, 5)
    phases['cubic_search_CPU_seconds'] = time.process_time() - cpu
    cpu = time.process_time()
    record = layout(model, result['record']['order'], 'conserving')
    del context, model
    phases['layout_and_cleanup_CPU_seconds'] = time.process_time() - cpu
    return dict(record=record, warm_start=warm, seeds=seeds, result=result, phases=phases,
                training_time=k.TRAIN_TIME, nominal_cubic_search_CPU_seconds=2.,
                scope='State-coherent cubic response;fullfirst-orderwarmconstructioncharged')


def cubic(model_path, digest, psi, reference_CPU_charge=0.):
    k.require(np.isfinite(reference_CPU_charge) and reference_CPU_charge >= 0, 'Invalid reference charge')
    cpu, wall = time.process_time(), time.monotonic()
    result = _cubic(model_path, digest, psi)
    elapsed = time.process_time() - cpu
    result.update(actual_local_CPU_seconds=elapsed, wall_seconds=time.monotonic() - wall,
                  accounting_overhead_CPU_seconds=elapsed - sum(result['phases'].values()),
                  reference_preparation_CPU_charge=reference_CPU_charge,
                  standalone_reference_and_ordering_CPU_charge=reference_CPU_charge + elapsed,
                  timing='Caller through complete warm search,cubic context/search and cleanup return')
    return result
