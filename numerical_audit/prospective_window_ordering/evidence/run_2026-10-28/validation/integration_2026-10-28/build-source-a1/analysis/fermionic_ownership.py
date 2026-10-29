"""Opt-in first-occurrence ownership with an explicit coefficient policy.

This constructs the real, number-conserving molecular schedule studied in the
manuscript. It does not replace or mutate the imported historical benchmarks.
Raw ladder words are normalized without an intermediate cutoff. Final JW
coefficients are accumulated directly from the raw Hamiltonian and cut off once.
Every retained Pauli is assigned to its first nonzero contributing canonical
parent. Scalar energy, discarded weight, and floating symmetry residuals are
returned explicitly; no Hermitian averaging or real projection is performed.

A nonzero reported number-commutator residual is a numerical diagnostic, not an
exact Lean certificate. Likewise, first ownership is discontinuous at zero:
retaining tiny parent coefficients can change the schedule even if it changes
its total Hamiltonian by very little.
"""

from collections import defaultdict
from dataclasses import dataclass
import math

from openfermion import FermionOperator

from qhat.analysis.fermionic_coefficients import (
    PRECISION_POLICY_VERSION, jordan_wigner_molecular, normal_ordered_molecular,
)


@dataclass(frozen=True)
class CanonicalMolecularParent:
    index: int
    first_raw_index: int
    operator: FermionOperator
    component_keys: tuple


@dataclass(frozen=True)
class OwnedPauliBlock:
    parent_index: int
    pauli_keys: tuple


@dataclass(frozen=True)
class MolecularOwnership:
    hamiltonian: FermionOperator
    scalar: complex
    parents: tuple
    parent_order: tuple
    final_coefficients: dict
    owned_blocks: tuple
    fallback_pauli_keys: tuple
    pauli_order: tuple
    diagnostics: dict


def _fsum(values):
    return complex(math.fsum(complex(v).real for v in values),
                   math.fsum(complex(v).imag for v in values))


def particle_number_commutator_l1(coefficients):
    """Pauli coefficient l1 norm of [N, H], an upper bound on operator norm.

    N = sum_j (I-Z_j)/2. For each transverse site, [N,X_j] = -iY_j
    and [N,Y_j] = iX_j. Contributions to the same word are consolidated.
    """
    contributions = defaultdict(list)
    for key, coefficient in coefficients.items():
        for position, (mode, label) in enumerate(key):
            if label == 'X':
                replacement, phase = 'Y', -1j
            elif label == 'Y':
                replacement, phase = 'X', 1j
            else:
                continue
            target = key[:position] + ((mode, replacement),) + key[position + 1:]
            contributions[target].append(phase * coefficient)
    return math.fsum(abs(_fsum(values)) for values in contributions.values())


def _commute(left, right):
    right_labels = dict(right)
    return sum(label != right_labels[mode] for mode, label in left if mode in right_labels) % 2 == 0


def _adjoint_key_and_phase(key):
    adjoint = tuple((mode, 1 - action) for mode, action in reversed(key))
    result = normal_ordered_molecular(FermionOperator(adjoint), atol=0)
    if len(result.terms) != 1:
        raise ValueError('Expected one nonzero normal-ordered adjoint component.')
    return next(iter(result.terms.items()))


def build_molecular_ownership(raw_operator, *, atol, ordering='signed_ascending'):
    """Build a molecular first-owner schedule and report its numerical scope.

    ``ordering`` is ``raw``, ``signed_ascending``, or ``magnitude_descending``.
    Signed/magnitude ties use the sorted canonical component keys. Inputs with
    more than four ladder operators per word, particle-changing canonical
    components, imaginary coefficients above atol, or adjoint mismatches above
    atol are rejected. Smaller mismatches are retained and reported, not fixed.
    The full identity coefficient is returned in ``scalar`` and excluded from
    ``final_coefficients``/``pauli_order``.
    """
    if ordering not in {'raw', 'signed_ascending', 'magnitude_descending'}:
        raise ValueError(f'Unsupported canonical-parent ordering {ordering!r}.')
    # These calls also validate finite coefficients, rank, and the final cutoff.
    canonical = normal_ordered_molecular(raw_operator, atol=0)
    full_jw = jordan_wigner_molecular(raw_operator, atol=0)
    retained_jw = jordan_wigner_molecular(raw_operator, atol=atol)
    final = {key: value for key, value in retained_jw.terms.items() if key}
    raw_indices = {key: index for index, key in enumerate(canonical.terms)}
    input_imaginary_max = max((abs(complex(value).imag) for value in canonical.terms.values()), default=0.0)
    if input_imaginary_max > atol:
        raise ValueError('Only real canonical fermionic coefficients are supported (up to reported atol residuals).')

    parents, visited = [], set()
    residuals = []
    missing_partners = 0
    for key, coefficient in canonical.terms.items():
        if sum(action for _, action in key) * 2 != len(key):
            raise ValueError('Only number-conserving canonical fermionic components are supported.')
        if key in visited:
            continue
        partner, phase = _adjoint_key_and_phase(key)
        partner_coefficient = canonical.terms.get(partner, 0j)
        mismatch = abs(partner_coefficient - complex(coefficient).conjugate() * phase)
        if mismatch > atol:
            raise ValueError(f'Hermitian-adjoint coefficient mismatch {mismatch} exceeds atol={atol}: {key}.')
        residuals.append(mismatch if partner == key else 2 * mismatch)
        if partner not in canonical.terms:
            missing_partners += 1
        visited.update((key, partner))
        if not key:  # Fermionic scalar is part of full_jw.scalar, not a parent.
            continue
        components = (key,) if partner == key else (key, partner)
        operator = FermionOperator()
        operator.terms = {component: canonical.terms.get(component, 0j) for component in components}
        parents.append(CanonicalMolecularParent(
            index=len(parents), first_raw_index=min(raw_indices.get(component, raw_indices[key]) for component in components),
            operator=operator, component_keys=components,
        ))

    def priority(index):
        keys = tuple(sorted(parents[index].component_keys))
        coefficient = complex(parents[index].operator.terms[keys[0]]).real
        return (coefficient if ordering == 'signed_ascending' else -abs(coefficient), keys)

    parent_order = list(range(len(parents)))
    if ordering != 'raw':
        parent_order.sort(key=priority)

    seen, blocks = set(), []
    for index in parent_order:
        mapped = jordan_wigner_molecular(parents[index].operator, atol=0)
        keys = tuple(key for key in mapped.terms if key in final and key not in seen)
        seen.update(keys)
        if keys:
            blocks.append(OwnedPauliBlock(index, keys))
    fallback = tuple(key for key in final if key not in seen)
    pauli_order = tuple(key for block in blocks for key in block.pauli_keys) + fallback

    noncommuting_blocks = sum(
        not all(_commute(a, b) for i, a in enumerate(block.pauli_keys) for b in block.pauli_keys[i + 1:])
        for block in blocks
    )
    block_number_residuals = [particle_number_commutator_l1({key: final[key] for key in block.pauli_keys}) for block in blocks]
    canonical_jw = jordan_wigner_molecular(canonical, atol=0)
    all_keys = full_jw.terms.keys() | canonical_jw.terms.keys()
    diagnostics = {
        'precision_policy_version': PRECISION_POLICY_VERSION,
        'ordering': ordering,
        'coefficient_cutoff': float(atol),
        'input_max_imaginary_coefficient': float(input_imaginary_max),
        'parent_hermiticity_residual_l1': math.fsum(residuals),
        'max_parent_hermiticity_residual_l1': max(residuals, default=0.0),
        'missing_adjoint_partners_below_tolerance': missing_partners,
        'pauli_imaginary_coefficient_l1': math.fsum(abs(complex(value).imag) for value in final.values()),
        'discarded_pauli_coefficient_l1': math.fsum(abs(value) for key, value in full_jw.terms.items() if key and key not in final),
        'normal_order_jw_residual_l1': math.fsum(abs(full_jw.terms.get(key, 0) - canonical_jw.terms.get(key, 0)) for key in all_keys),
        'number_of_parents': len(parents),
        'number_of_owned_blocks': len(blocks),
        'fallback_pauli_count': len(fallback),
        'noncommuting_owned_block_count': noncommuting_blocks,
        'hamiltonian_number_commutator_l1': particle_number_commutator_l1(final),
        'max_owned_block_number_commutator_l1': max(block_number_residuals, default=0.0),
    }
    return MolecularOwnership(canonical, full_jw.terms.get((), 0j), tuple(parents), tuple(parent_order),
                              final, tuple(blocks), fallback, pauli_order, diagnostics)
