"""Exact dyadic molecular parents and weighted Pauli-block combinatorics.

Python arbitrary-precision integers establish zero identities for the frozen
binary coefficients. This is independent of floating-point simulation and is
not a Lean certificate or a claim that floating molecular tensors obey every
analytic Coulomb symmetry exactly.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import math

from qhat.analysis.fermionic_coefficients import _unit_transform


def add_into(target, source, factor=1):
    for key, value in source.items():
        total = target.get(key, 0) + factor * value
        if total:
            target[key] = total
        else:
            target.pop(key, None)
    return target


@lru_cache(maxsize=32768)
def adjoint(key):
    rev = tuple((q, 1-a) for q, a in reversed(key))
    result = _unit_transform(rev, 'normal')
    if len(result) != 1 or result[0][1] not in (1, -1):
        raise ValueError('Canonical balanced component must have one signed adjoint.')
    return result[0][0], int(result[0][1].real)


def exact_parents(raw, cutoff=1e-12):
    """Hermitian part, whole-parent cutoff, exact JW; one new common target.

    Input complex parts must be identically zero. The Hermitian part uses an
    exact arithmetic average of paired canonical coefficients, not selection
    of a preferred coefficient. This correction and the subsequent whole-
    parent truncation are reported separately from the old/raw target.
    All maps use a single positive dyadic denominator and integer numerators.
    """
    pairs = {}
    for key, val in raw.terms.items():
        value = complex(val)
        if not math.isfinite(value.real) or value.imag != 0:
            raise ValueError('Expected finite exactly real molecular input.')
        pairs[key] = value.real.as_integer_ratio()
    denominator = max((d for _, d in pairs.values()), default=1) * 32
    normal = defaultdict(int)
    for key, (num, den) in pairs.items():
        integer = num * (denominator // den)
        for target, unit in _unit_transform(key, 'normal'):
            if complex(unit).imag or int(unit.real) != unit.real:
                raise ValueError('Unit normal ordering is not integral.')
            normal[target] += integer * int(unit.real)
    normal = {k: v for k, v in normal.items() if v}
    visited, parents = set(), []
    adjustment, discarded = 0, 0
    for key in sorted(normal):
        if key in visited:
            continue
        if sum(a for _, a in key)*2 != len(key):
            raise ValueError('Unbalanced canonical source.')
        partner, phase = adjoint(key)
        visited.update((key, partner))
        first, second = normal[key], phase*normal.get(partner, 0)
        if (first+second) % 2:
            raise ValueError('Insufficient dyadic scaling for exact averaging.')
        mean = (first+second)//2
        adjustment += abs(first-mean)
        if partner != key:
            adjustment += abs(second-mean)
        if not key:
            continue
        if abs(mean/denominator) <= cutoff:
            discarded += abs(mean)*(1 if key == partner else 2)
            continue
        components = {key: mean}
        if key != partner:
            components[partner] = phase*mean
        mapped, imaginary = defaultdict(int), defaultdict(int)
        for component, weight in components.items():
            for word, unit in _unit_transform(component, 'jw'):
                re, im = complex(unit).real.as_integer_ratio(), complex(unit).imag.as_integer_ratio()
                if weight % re[1] or weight % im[1]:
                    raise ValueError('Insufficient dyadic scaling for unit JW coefficients.')
                mapped[word] += weight//re[1]*re[0]
                imaginary[word] += weight//im[1]*im[0]
        if any(imaginary.values()):
            raise ValueError('Hermitian real parent has imaginary JW coefficient.')
        terms = {k: v for k,v in mapped.items() if v and k}
        if not terms:
            continue
        parents.append({'components': tuple(sorted(components)), 'weight': mean,
                        'coefficients': terms,
                        'scalar': mapped.get((), 0),
                        'four_distinct': len(key) == 4 and len({q for q,_ in key}) == 4})
    parents.sort(key=lambda p: (p['weight'], p['components']))
    total = {}
    for p in parents:
        add_into(total, p['coefficients'])
    return parents, total, denominator, {
        'canonical_hermitian_adjustment_fermionic_l1': adjustment/denominator,
        'discarded_parent_fermionic_l1': discarded/denominator,
        'input_constant': complex(raw.terms.get((), 0)).real,
        'removed_parent_pauli_scalar': sum(p['scalar'] for p in parents)/denominator,
        'parent_cutoff': cutoff,
        'policy': 'exact real-input CAR; exact adjoint average; whole-parent cutoff; exact JW consolidation; no Pauli cutoff',
    }


@lru_cache(maxsize=131072)
def symplectic(word):
    x = z = ny = 0
    for mode, label in word:
        bit = 1 << mode
        if label in ('X', 'Y'): x |= bit
        if label in ('Y', 'Z'): z |= bit
        if label == 'Y': ny += 1
    return x,z,ny


@lru_cache(maxsize=131072)
def from_symplectic(x,z):
    out = []
    remaining = x | z
    while remaining:
        bit = remaining & -remaining
        q = bit.bit_length()-1
        out.append((q, 'Y' if x&bit and z&bit else ('X' if x&bit else 'Z')))
        remaining ^= bit
    return tuple(out)


def pauli_commutator_term(left, right):
    """Return (word, real coefficient of i) for [left,right]."""
    x,z,y = symplectic(left); xx,zz,yy = symplectic(right)
    if ((x&zz).bit_count()+(z&xx).bit_count()) % 2 == 0:
        return None, 0
    dest = from_symplectic(x^xx,z^zz)
    phase = (y+yy-((x^xx)&(z^zz)).bit_count()+2*(z&xx).bit_count()) % 4
    if phase not in (1,3):
        raise ValueError('Anticommuting Hermitian Paulis must have imaginary product.')
    return dest, 2 if phase == 1 else -2


def commutator(left, right):
    """Integer coefficients of [A,B]/i; retain cancellation before norms."""
    values = defaultdict(int)
    for a, va in left.items():
        for b, vb in right.items():
            word, factor = pauli_commutator_term(a,b)
            if factor:
                values[word] += factor*va*vb
    return {word: value for word,value in values.items() if value}


def internally_commuting(block):
    return all(pauli_commutator_term(a,b)[1] == 0 for a,b in combinations(block,2))


def number_commutator(block):
    """Integer coefficients of [N,A]/i, denominator equals A's denominator."""
    values = defaultdict(int)
    for word, v in block.items():
        for i,(q,label) in enumerate(word):
            if label not in ('X','Y'): continue
            dest = word[:i]+((q,'Y' if label=='X' else 'X'),)+word[i+1:]
            values[dest] += -v if label == 'X' else v
    return {k:v for k,v in values.items() if v}


def decompose(parents, total):
    seen, owned = set(), []
    for p in parents:
        block = {k:total[k] for k in p['coefficients'] if k in total and k not in seen}
        seen.update(block)
        if block: owned.append(block)
    if seen != set(total):
        raise ValueError('Ownership not exhaustive.')
    groups = defaultdict(dict)
    for key,value in sorted(total.items()):
        groups[symplectic(key)[0]][key] = value
    masks = sorted(groups)
    return [p['coefficients'] for p in parents], owned, [groups[m] for m in masks], masks


def graph_and_color(blocks, keep_commutators=False):
    """Exact operator graph, deterministic largest-degree greedy coloring."""
    adjacent = [set() for _ in blocks]
    comms = {}
    for i,j in combinations(range(len(blocks)),2):
        c = commutator(blocks[i], blocks[j])
        if c:
            adjacent[i].add(j); adjacent[j].add(i)
            if keep_commutators: comms[i,j] = c
    visit = sorted(range(len(blocks)), key=lambda i:(-len(adjacent[i]), i))
    colors = {}
    for i in visit:
        unavailable = {colors[j] for j in adjacent[i] if j in colors}
        c = 0
        while c in unavailable: c += 1
        colors[i] = c
    groups = [[i for i in range(len(blocks)) if colors[i] == c] for c in range(max(colors.values(),default=-1)+1)]
    for group in groups:
        if any(j in adjacent[i] for i,j in combinations(group,2)):
            raise ValueError('Improper coloring.')
    info = {'vertex_count':len(blocks), 'edge_count':sum(map(len,adjacent))//2,
            'colors':groups,'color_count':len(groups),
            'edge_list':[[i,j] for i,ns in enumerate(adjacent) for j in sorted(ns) if i<j],
            'certificate':'exact arbitrary-precision integer Pauli commutators; no zero tolerance; not Lean'}
    return [i for group in groups for i in group], info, comms


def suffix_objective(order, comms):
    result = 0
    for at,i in enumerate(order):
        acc = {}
        for j in order[at+1:]:
            add_into(acc, comms.get((min(i,j),max(i,j)),{}), 1 if i<j else -1)
        result += sum(abs(v) for v in acc.values())
    return result


def suffix_greedy(count, comms):
    """Greedy L1 suffix-bound candidate, followed by strictly reducing swaps.

    Lambda_hat=sum ||Pauli coefficients of [G_i,suffix]||_1 bounds the
    sector operator-norm Lambda. Signs/cancellations within each commutator
    are retained. A decrease is only a decrease of a certified upper bound,
    not a prediction or guarantee that actual error decreases.
    """
    suffixes = [{} for _ in range(count)]
    for (i,j),c in comms.items():
        add_into(suffixes[i],c); add_into(suffixes[j],c,-1)
    remaining, order = set(range(count)), []
    while remaining:
        i = min(remaining,key=lambda k:(sum(abs(v) for v in suffixes[k].values()),k))
        order.append(i); remaining.remove(i)
        for j in remaining:
            add_into(suffixes[j],comms.get((min(i,j),max(i,j)),{}), -1 if j<i else 1)
    before = suffix_objective(order,comms)
    accepted = 0
    # Two deterministic passes; bounded candidate generation, not local-optimality claim.
    for _ in range(2):
        for k in range(len(order)-1):
            i,j = order[k:k+2]
            suffix = order[k+2:]
            ai,bj = {},{}
            for v in suffix:
                add_into(ai,comms.get((min(i,v),max(i,v)),{}),1 if i<v else -1)
                add_into(bj,comms.get((min(j,v),max(j,v)),{}),1 if j<v else -1)
            c = comms.get((min(i,j),max(i,j)),{})
            sign = 1 if i<j else -1
            aplus = dict(ai); add_into(aplus,c,sign)
            bminus = dict(bj); add_into(bminus,c,-sign)
            old = sum(map(abs,aplus.values()))+sum(map(abs,bj.values()))
            new = sum(map(abs,ai.values()))+sum(map(abs,bminus.values()))
            if new < old:
                order[k:k+2] = [j,i]; accepted += 1
    after = suffix_objective(order,comms)
    if after > before: raise ValueError('Objective increased.')
    return order, {'greedy_objective_integer':str(before),'final_objective_integer':str(after),
                   'accepted_adjacent_swaps':accepted,'adjacent_passes':2}
