"""
Reference implementation for Trotter error coefficient estimation.

Based on: Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)

This file provides the original (slower) implementation for:
- Testing and validation
- Reference for understanding the algorithm
- Comparison against the optimized version

For production use, prefer trotter_coefficients_fast.py. Historical speedups
were measured on the earlier index domain; the corrected domain is not rebenchmarked here.
"""

import logging
import time
import numpy as np
import math

logger = logging.getLogger(__name__)

# --------------------------------------------------
# Helper functions for Pauli string arithmetic
# --------------------------------------------------
def pauli_dict(pauli_key):
    """Convert a tuple of (qubit, operator) pairs into a dictionary."""
    return {q: op for (q, op) in pauli_key}

def anticommute(key1, key2):
    """
    Return True if the Pauli strings (given by key1 and key2)
    anticommute. They anticommute if the number of qubits
    on which they both act nontrivially with different operators is odd.
    """
    d1 = pauli_dict(key1)
    d2 = pauli_dict(key2)
    count = 0
    for q in d1:
        if q in d2 and d1[q] != d2[q]:
            count += 1
    return (count % 2 == 1)

def pauli_product(p1, p2):
    """
    Multiply two single-qubit Pauli operators (ignoring phase).
    """
    if p1 == 'I':
        return p2
    if p2 == 'I':
        return p1
    if p1 == p2:
        return 'I'
    # For different Paulis, the unique remaining operator is returned.
    if {p1, p2} == {'X', 'Y'}:
        return 'Z'
    if {p1, p2} == {'X', 'Z'}:
        return 'Y'
    if {p1, p2} == {'Y', 'Z'}:
        return 'X'
    return 'I'

def pauli_product_key(key1, key2):
    """
    Given two Pauli strings (each represented as a tuple of (qubit, op)),
    return the Pauli string (as a sorted tuple) corresponding to their product.
    The product is computed qubit-by-qubit (ignoring any phase).
    """
    d1 = pauli_dict(key1)
    d2 = pauli_dict(key2)
    prod = {}
    for q in set(d1.keys()).union(d2.keys()):
        op1 = d1.get(q, 'I')
        op2 = d2.get(q, 'I')
        prod[q] = pauli_product(op1, op2)
    # Remove qubits where the product is identity.
    prod = {q: op for q, op in prod.items() if op != 'I'}
    return tuple(sorted(prod.items()))

# --------------------------------------------------
# Functions to compute (nested) commutator norms directly
# --------------------------------------------------
def compute_commutator(op1, op2):
    """
    Given two QubitOperator terms op1 and op2 (each with one Pauli string),
    compute their commutator [op1, op2].

    Returns:
      - (None, 0.0) if they commute;
      - ((prod_key, coeff), norm) if they anticommute,
        where prod_key is the resulting Pauli string (as a tuple)
        and coeff = 2*c1*c2 (with c1, c2 the term coefficients).
    The norm is computed as 2*|c1*c2|.
    """
    key1 = list(op1.terms.keys())[0]
    key2 = list(op2.terms.keys())[0]
    c1 = list(op1.terms.values())[0]
    c2 = list(op2.terms.values())[0]
    if not anticommute(key1, key2):
        return None, 0.0
    prod_key = pauli_product_key(key1, key2)
    coeff = 2 * c1 * c2
    norm_val = 2 * abs(c1 * c2)
    return (prod_key, coeff), norm_val

def compute_nested_commutator_norm(H_outer, inner_operator):
    """
    Given an outer QubitOperator H_outer and an inner operator (a tuple (prod_key, coeff))
    representing the inner commutator [H_j, H_k], compute the norm of the nested commutator:
        [H_outer, inner_operator].
    If the Pauli string of H_outer anticommutes with inner_operator's pauli key,
    the norm is 2*|c_outer * (inner coeff)|; otherwise, it is zero.
    """
    if inner_operator is None:
        return 0.0
    key_outer = list(H_outer.terms.keys())[0]
    c_outer = list(H_outer.terms.values())[0]
    if not anticommute(key_outer, inner_operator[0]):
        return 0.0
    return 2 * abs(c_outer * inner_operator[1])

# --------------------------------------------------
# Main Monte Carlo estimation
# --------------------------------------------------
def trotter_error_estimator(pauli_terms, time_limit, batch_size=100):
    """
    Monte Carlo estimation of Trotter error coefficients.

    Reference: Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)

    For first-order formulas (Equation 145):
      Error ≤ (t²/2) * C1, where C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||

    For second-order formulas (Equation 152):
      Error ≤ (t³/12) * C21 + (t³/24) * C22, where
      C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||  (all triples with k < i and k < j)
      C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||       (all pairs with k < j)

    Returns:
      (C1/2, C2) where C2 = C21/12 + C22/24. These are Monte Carlo
      estimates, not certified upper bounds. C21 includes ordered i,j and i=j.
    """
    N = len(pauli_terms)
    if N < 2:
        return 0.0, 0.0
    if time_limit <= 0 or batch_size < 1:
        raise ValueError("Positive time_limit and batch_size are required")

    def pair():
        i = np.random.randint(N)
        j = np.random.randint(N-1)
        if j >= i:
            j += 1
        return min(i, j), max(i, j)

    def nested():
        # Independent uniform rejection sampler, including equal outer/inner
        # indices. No fast-estimator helper or exact result is used here.
        while True:
            i, j, k = np.random.randint(N, size=3)
            if k < i and k < j:
                return i, j, k

    def estimate(draw, count):
        total = 0.0
        samples = 0
        start = time.monotonic()
        while time.monotonic() - start < time_limit/3:
            for _ in range(batch_size):
                total += draw()
                samples += 1
        return total * count/samples if samples else 0.0

    def first():
        i, j = pair()
        return compute_commutator(pauli_terms[i], pauli_terms[j])[1]

    def second21():
        i, j, k = nested()
        inner, _ = compute_commutator(pauli_terms[j], pauli_terms[k])
        return compute_nested_commutator_norm(pauli_terms[i], inner)

    def second22():
        k, j = pair()
        inner, _ = compute_commutator(pauli_terms[k], pauli_terms[j])
        return compute_nested_commutator_norm(pauli_terms[k], inner)

    c1 = estimate(first, N*(N-1)//2)
    c21 = estimate(second21, N*(N-1)*(2*N-1)//6)
    c22 = estimate(second22, N*(N-1)//2)
    return c1/2, c21/12+c22/24
