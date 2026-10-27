"""Exact complete-unit actions, invariant sectors and ordered BCH coefficients.

Modeq is bit1<<q; coefficients have Hartree units, time has inverse-Hartree
units and hbar=1. No floating-point matrix entries are used by certificates.
"""
from fractions import Fraction as F
from intervals import require


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def add(a, b):
    return [[x+y for x,y in zip(u,v)] for u,v in zip(a,b)]


def times(a, scalar):
    return [[x*scalar for x in row] for row in a]


def multiply(a, b):
    n = len(a)
    out = zero(n)
    for i,row in enumerate(a):
        for k,x in enumerate(row):
            if x:
                for j,y in enumerate(b[k]):
                    if y:
                        out[i][j] += x*y
    return out


def commutator(a, b):
    return add(multiply(a,b), times(multiply(b,a),-1))


def entry_l1(a):
    """Conservative upper bound on the Euclidean induced operator norm."""
    return sum((abs(x) for row in a for x in row), F(0))


def encoded(a):
    return [[str(x) for x in row] for row in a]


def decode(rows, n):
    require(type(rows) is list and 0 < len(rows) <= 128, 'Bounded nonempty unit required')
    block = {}
    for word, coefficient in rows:
        w = tuple((q,p) for q,p in word)
        require(w and all(type(q) is int and 0 <= q < n and p in ('X','Y','Z') for q,p in w)
                and [q for q,p in w] == sorted({q for q,p in w}), 'Invalid Pauli word')
        require(type(coefficient) in (int,str), 'Exact integer coefficient required')
        c = int(coefficient)
        require(type(coefficient) is int or str(c) == coefficient, 'Canonical integer coefficient required')
        require(c and c.bit_length() <= 512 and w not in block, 'Invalid or duplicate coefficient')
        block[w] = c
    return block


def descriptor(block):
    xmask, rows = None, []
    for word,c in block.items():
        x = sum(1 << q for q,p in word if p in ('X','Y'))
        z = sum(1 << q for q,p in word if p in ('Z','Y'))
        ny = sum(p == 'Y' for q,p in word)
        require(ny % 2 == 0 and (xmask is None or xmask == x), 'Real equal-flip unit required')
        xmask = x
        rows.append((z, ny//2, c))
    for z,y,c in rows:
        for other_z,other_y,other_c in rows:
            require(((xmask & other_z).bit_count() + (xmask & z).bit_count()) % 2 == 0,
                    'Unit Pauli terms must commute')
    return xmask, rows


def action(desc, occupation):
    return sum(c * (-1 if ((z & occupation).bit_count()+nyhalf) % 2 else 1)
               for z,nyhalf,c in desc[1])


def spin_counts(v, n):
    return tuple(sum((v >> q) & 1 for q in range(parity,n,2)) for parity in (0,1))


def from_model(model, coefficient_scale=F(1)):
    n, electrons = model['n'], model['electrons']
    require(type(n) is int and 2 <= n <= 8 and n % 2 == 0 and
            type(electrons) is int and 0 < electrons <= n and electrons % 2 == 0,
            'Dimension/electron guard before allocation')
    require(type(model['blocks']) is list and 0 < len(model['blocks']) <= 128, 'Unit count guard')
    den = int(model['denominator'])
    require(0 < den and den.bit_length() <= 512 and type(coefficient_scale) in (int,F),
            'Bounded positive denominator and exact coefficient scale required')
    coefficient_scale = F(coefficient_scale)
    require(coefficient_scale > 0, 'Positive coefficient scale required')
    blocks = [decode(rows,n) for rows in model['blocks']]
    total, nonzero_masks = {}, set()
    desc = [descriptor(block) for block in blocks]
    for block,(mask,rows) in zip(blocks,desc):
        require(mask != 0 or len(block) == 1, 'Keep diagonal terms individually')
        require(mask == 0 or mask not in nonzero_masks, 'Complete unique flip classes required')
        if mask:
            nonzero_masks.add(mask)
        for word,c in block.items():
            require(word not in total, 'Overlapping unit partition')
            total[word] = c
        # Check all full-space columns before constructing a restricted matrix.
        for v in range(1 << n):
            amplitude = action((mask,rows),v)
            require(amplitude == action((mask,rows),v ^ mask), 'Hermitian unit action required')
            require(not amplitude or spin_counts(v,n) == spin_counts(v ^ mask,n),
                    'Complete unit fails full-space spin-sector conservation')
    require(total == decode(model['coefficients'],n), 'Unit partition differs from total Hamiltonian')
    hf = (1 << electrons)-1
    sector = [v for v in range(1 << n) if spin_counts(v,n) == spin_counts(hf,n)]
    require(len(sector) <= 36, 'Invariant-sector dimension ceiling')
    indices = {v:i for i,v in enumerate(sector)}
    matrices = []
    for mask,rows in desc:
        a = zero(len(sector))
        for j,v in enumerate(sector):
            amplitude = action((mask,rows),v)
            if amplitude:
                require(v ^ mask in indices, 'Sector action missing after full-space check')
                a[indices[v ^ mask]][j] = F(amplitude,den)*coefficient_scale
        require(a == [list(row) for row in zip(*a)], 'Projected matrix is not real symmetric')
        matrices.append(a)
    diag = sorted((i for i,d in enumerate(desc) if not d[0]),
                  key=lambda i:(-abs(next(iter(blocks[i].values()))),next(iter(blocks[i]))))
    off = sorted((i for i,d in enumerate(desc) if d[0]),
                 key=lambda i:(-sum(abs(c) for c in blocks[i].values()),desc[i][0]))
    require(diag and off, 'Molecular guard requires diagonal and off-diagonal units')
    return dict(n=n,electrons=electrons,den=den,blocks=blocks,descriptors=desc,sector=sector,
                hf=hf,hf_index=indices[hf],matrices=matrices,diag=diag,off=off)


def bch(data, order):
    """Return imaginary coefficient matrices of K1,K2: K=-i/2 sum[i<j][Gi,Gj]."""
    m, n = len(data['matrices']),len(data['sector'])
    require(sorted(order) == list(range(m)), 'Invalid unit permutation')
    k1,k2 = zero(n),zero(n)
    diagonal = set(data['diag'])
    for at,i in enumerate(order):
        for j in order[at+1:]:
            if i in diagonal and j in diagonal:
                continue
            c = times(commutator(data['matrices'][i],data['matrices'][j]),F(-1,2))
            if (i in diagonal) != (j in diagonal):
                k1 = add(k1,c)
            else:
                k2 = add(k2,c)
    for a in (k1,k2):
        require(all(a[i][j] == -a[j][i] for i in range(n) for j in range(n)), 'BCH skew symmetry')
    return k1,k2
