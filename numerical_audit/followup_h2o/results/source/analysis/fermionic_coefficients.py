"""Explicit coefficient cutoffs for molecular fermion/Jordan--Wigner transforms.

OpenFermion's operator additions can discard coefficients below its own default
``EQ_TOLERANCE`` before a later ``compress(abs_tol=...)`` call. These opt-in
routines transform unit-weight monomials, accumulate their physical coefficients
in dictionaries, and apply only the caller's final absolute cutoff. They leave
OpenFermion's process-wide behavior and the historical benchmark unchanged.

The supported input consists of words of at most four ladder operators (a scalar,
one-body, or two-body molecular Hamiltonian). Unit-word normal ordering has
integer coefficients; its JW expansion has exactly representable dyadic
coefficients of magnitude at least 1/16 when nonzero. Thus OpenFermion's internal
1e-8 arithmetic pruning cannot discard any such unit-word coefficient. Floating
point summation remains floating point; these routines do not certify exact
arithmetic or an error bound for discarded terms.
"""

from collections import defaultdict
from functools import lru_cache
import math

import numpy as np
from openfermion import FermionOperator, InteractionOperator, QubitOperator
from openfermion import jordan_wigner, normal_ordered


PRECISION_POLICY_VERSION = "molecular-unit-weight-v1"
MAX_LADDER_OPERATORS = 4


def _coefficient(value):
    result = complex(value)
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError("Molecular operator coefficients must be finite.")
    return result


def _validate_cutoff(atol):
    if not math.isfinite(atol) or atol < 0:
        raise ValueError("atol must be a finite, nonnegative absolute cutoff.")


@lru_cache(maxsize=16384)
def _unit_transform(key, mapping):
    if len(key) > MAX_LADDER_OPERATORS:
        raise ValueError(
            "Explicit molecular transforms support at most four ladder "
            f"operators per monomial; received {len(key)}."
        )
    unit = FermionOperator(key, 1.0)
    transformed = normal_ordered(unit) if mapping == "normal" else jordan_wigner(unit)
    return tuple(transformed.terms.items())


def _transform(operator, mapping, atol):
    if not isinstance(operator, FermionOperator):
        raise TypeError("Expected an OpenFermion FermionOperator.")
    _validate_cutoff(atol)
    contributions = defaultdict(list)
    for key, value in operator.terms.items():
        coefficient = _coefficient(value)
        # Validate rank even for zero coefficients: the supported domain is
        # explicit and must not depend on accidental numerical cancellation.
        for target, unit_coefficient in _unit_transform(key, mapping):
            contributions[target].append(coefficient * unit_coefficient)

    result = FermionOperator() if mapping == "normal" else QubitOperator()
    for key, values in contributions.items():
        coefficient = complex(
            math.fsum(value.real for value in values),
            math.fsum(value.imag for value in values),
        )
        if abs(coefficient) > atol:
            # Do not use +=: that would reintroduce OpenFermion's cutoff.
            result.terms[key] = coefficient
    return result


def normal_ordered_molecular(operator: FermionOperator, *, atol: float) -> FermionOperator:
    """CAR-normalize words of length <=4, then cut off final coefficients."""
    return _transform(operator, "normal", atol)


def jordan_wigner_molecular(operator: FermionOperator, *, atol: float) -> QubitOperator:
    """Map words of length <=4 to JW Paulis, then cut off final coefficients.

    Inputs need not already be normal ordered. Applying this directly to the
    raw molecular operator avoids an intermediate normal-order cutoff.
    """
    return _transform(operator, "jw", atol)


def fermion_from_interaction(interaction: InteractionOperator) -> FermionOperator:
    """Read molecular integrals without a coefficient cutoff or implicit 1/2.

    Coefficients follow OpenFermion InteractionOperator's convention:
    constant + h[p,q] a_p^dagger a_q + g[p,q,r,s] a_p^dagger a_q^dagger a_r a_s.
    Only exact zero entries are omitted. Choose any truncation explicitly in a
    subsequent transform, after contributions have been combined.
    """
    if not isinstance(interaction, InteractionOperator):
        raise TypeError("Expected an OpenFermion InteractionOperator.")
    one_body = np.asarray(interaction.one_body_tensor)
    two_body = np.asarray(interaction.two_body_tensor)
    if one_body.ndim != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError("one_body_tensor must be a square matrix.")
    if two_body.shape != (one_body.shape[0],) * 4:
        raise ValueError("two_body_tensor must have four matching orbital axes.")
    result = FermionOperator()
    constant = _coefficient(interaction.constant)
    if constant != 0:
        result.terms[()] = constant
    for indices, value in np.ndenumerate(one_body):
        coefficient = _coefficient(value)
        if coefficient != 0:
            p, q = indices
            result.terms[((p, 1), (q, 0))] = coefficient
    for indices, value in np.ndenumerate(two_body):
        coefficient = _coefficient(value)
        if coefficient != 0:
            p, q, r, s = indices
            result.terms[((p, 1), (q, 1), (r, 0), (s, 0))] = coefficient
    return result
