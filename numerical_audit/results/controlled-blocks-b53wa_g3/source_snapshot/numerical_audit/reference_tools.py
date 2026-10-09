"""Exact Pauli-action helpers adapted from Albert's pinned L-sweep mechanism audit.

Source: 89c701652c4c64918825e279e84a7bea5924822a,
codex_mechanism_experiment/mechanism_experiment.py.
The affine subspace is closed under every Pauli flip, so it retains leakage
outside the physical number/spin sector; it is not postselection.
"""
import numpy as np
from scipy import sparse

def affine_space(terms, hf_index):
    """Exact closure under every flip; coordinate bit k is basis[k]."""
    pivots = {}
    for _, flip, _, _ in terms:
        x = int(flip)
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    basis = [pivots[p] for p in sorted(pivots)]
    states = [int(hf_index)]
    for b in basis:
        states += [x ^ b for x in states]
    indices = np.asarray(states, dtype=np.int64)
    coordinates = {int(x) ^ int(hf_index): i for i, x in enumerate(indices)}
    reduced = []
    for c, x, z, y in terms:
        zz = sum(((b & int(z)).bit_count() % 2) << k for k, b in enumerate(basis))
        cc = c * (-1 if (int(hf_index) & int(z)).bit_count() % 2 else 1)
        reduced.append((cc, np.uint64(coordinates[int(x)]), np.uint64(zz), y))
    return indices, reduced, basis


def apply_p(state, term, parity):
    _, x, z, y = term
    src = np.arange(state.size, dtype=np.int64) ^ int(x)
    signs = 1 - 2 * parity[src & int(z)]
    return (1j ** int(y)) * signs * state[src]


def assemble_hamiltonian(terms, dimension):
    """Assemble exactly the retained Pauli sum; discard only exact zero entries."""
    coords = np.arange(dimension, dtype=np.int64)
    parity = np.array([int(x).bit_count() % 2 for x in coords], dtype=np.int8)
    grouped = {}
    for c, x, z, y in terms:
        x = int(x)
        if x not in grouped:
            grouped[x] = np.zeros(dimension, dtype=np.complex128)
        grouped[x] += c * (1j ** int(y)) * (1 - 2 * parity[coords & int(z)])
    rows, cols, data = [], [], []
    for x, values in grouped.items():
        keep = np.flatnonzero(values != 0)
        rows.append(keep ^ x)
        cols.append(keep)
        data.append(values[keep])
    H = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(dimension, dimension)).tocsr()
    if sparse.linalg.norm(H - H.getH()) > 1e-10:
        raise RuntimeError('Assembled Hamiltonian is not Hermitian')
    return H, parity


def leading_defect(terms, initial, H, parity):
    """Coefficient of dt^2 in S(dt)|HF>-exp(-iHdt)|HF>."""
    first = np.zeros_like(initial)
    second = np.zeros_like(initial)
    for term in terms:
        c = term[0]
        second -= 1j * c * apply_p(first, term, parity)
        second -= 0.5 * c * c * initial
        first -= 1j * c * apply_p(initial, term, parity)
    return second + 0.5 * (H @ (H @ initial))


def norm2(v):
    return float(np.vdot(v, v).real)


def defect_columns(terms,sector,H,parity):
    dim=H.shape[0];width=len(sector)
    E=np.zeros((dim,width),complex);E[sector,np.arange(width)]=1
    first=np.zeros_like(E);second=np.zeros_like(E)
    coords=np.arange(dim,dtype=np.int64)
    for c,x,z,y in terms:
        src=coords ^ int(x)
        phase=(1j**int(y))*(1-2*parity[src & int(z)])
        second-=1j*c*phase[:,None]*first[src]
        second[sector,np.arange(width)]-=.5*c*c
        destination=sector ^ int(x)
        initial_phase=(1j**int(y))*(1-2*parity[sector & int(z)])
        first[destination,np.arange(width)]-=1j*c*initial_phase
    return second+.5*(H@(H@E))

