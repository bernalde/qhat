"""Finite polynomial propagation with separate outward and analytic enclosures."""
from fractions import Fraction as F
from intervals import Complex, Real, exp_i, require, ceil_div
from exact import zero, add, entry_l1


def point_zero(z):
    return z.real.lo == z.real.hi == z.imag.lo == z.imag.hi == 0


def identity(n, bits):
    return [[Complex.exact(int(i==j),bits=bits) for j in range(n)] for i in range(n)]


def matmul(a, b):
    n,bits = len(a),a[0][0].bits
    out = [[Complex.exact(bits=bits) for _ in range(n)] for _ in range(n)]
    for i,row in enumerate(a):
        for k,x in enumerate(row):
            if not point_zero(x):
                for j,y in enumerate(b[k]):
                    if not point_zero(y):
                        out[i][j] = out[i][j]+x*y
    return out


def matvec(a, v):
    out = []
    for row in a:
        value = Complex.exact(bits=v[0].bits)
        for x,y in zip(row,v):
            if not point_zero(x) and not point_zero(y):
                value = value+x*y
        out.append(value)
    return out


def power(a, exponent):
    require(type(exponent) is int and 0 <= exponent <= 4096, 'Finite matrix power bound')
    result = identity(len(a),a[0][0].bits)
    while exponent:
        if exponent & 1:
            result = matmul(result,a)
        exponent //= 2
        if exponent:
            a = matmul(a,a)
    return result


def radius_upper(value, bits):
    value = F(value)
    require(value >= 0, 'Nonnegative analytic radius')
    return F(ceil_div(value.numerator*(1 << bits),value.denominator),1 << bits)


def combine(a, b, bits):
    """Error of a product of approximate factors with contractive exact factors."""
    return radius_upper(a+b+a*b,bits)


def exp_matrix(h, tau, *, terms=96, bits=192):
    """Generic real-symmetric-H exponential, exact polynomial enclosed in intervals."""
    require(type(tau) in (int,F) and type(terms) is int and 2 <= terms <= 256, 'Exact time/term bound')
    n = len(h)
    require(0 < n <= 36 and all(len(row)==n for row in h), 'Matrix dimension bound')
    require(h == [list(row) for row in zip(*h)], 'Real symmetric generator required')
    norm = abs(tau)*entry_l1(h)
    scaling = 1
    while 2*norm > scaling:
        scaling *= 2
        require(scaling <= 4096, 'Matrix exponential scaling ceiling')
    if norm == 0:
        return identity(n,bits),F(0),dict(scaling=1,terms=terms,norm_upper='0')
    generator = [[Complex.exact(0,-tau*x/scaling,bits) for x in row] for row in h]
    term = result = identity(n,bits)
    for degree in range(1,terms):
        term = [[x.times(F(1,degree)) for x in row] for row in matmul(generator,term)]
        result = [[x+y for x,y in zip(u,v)] for u,v in zip(result,term)]
    result = power(result,scaling)
    local = F(1,1 << (terms-1))
    radius = radius_upper((1+local)**scaling-1,bits)
    return result,radius,dict(scaling=scaling,terms=terms,norm_upper=str(norm),analytic_radius=str(radius))


def exp_unit(h, tau, *, terms=96, bits=192):
    """Disjoint diagonal/weighted-flip blocks; no approximate state normalization."""
    n = len(h)
    result = identity(n,bits)
    used = set()
    radius = F(0)
    for j in range(n):
        if j in used:
            continue
        support = [i for i in range(n) if h[i][j]]
        require(len(support) <= 1, 'Complete real unit must have one output per column')
        if not support:
            used.add(j)
            continue
        i = support[0]
        a = h[i][j]
        require(h[j][i] == a, 'Hermitian complete unit required')
        value,tail,_ = exp_i(-tau*a,terms=terms,bits=bits)
        radius = max(radius,tail)
        if i == j:
            result[j][j] = value
            used.add(j)
        else:
            require(not h[j][j] and not h[i][i] and sum(bool(x) for x in h[i]) == 1,
                    'Disjoint real flip block required')
            result[i][i] = result[j][j] = Complex(value.real,Real.exact(0,bits))
            result[i][j] = result[j][i] = Complex(Real.exact(0,bits),value.imag)
            used.update((i,j))
    return result,radius_upper(radius,bits)


def leading_response(data, k1, k2, time, *, bits=192):
    n,hf = len(data['sector']),data['hf_index']
    diagonal,v = zero(n),zero(n)
    for i,h in enumerate(data['matrices']):
        if i in data['diag']:
            diagonal = add(diagonal,h)
        else:
            v = add(v,h)
    variance = Real.exact(0,bits)
    components = []
    for j in range(n):
        if j == hf or not k1[j][hf]:
            continue
        gap = diagonal[j][j]-diagonal[hf][hf]
        if gap:
            e,tail,_ = exp_i(time*gap,bits=bits)
            freal = e.imag.expanded(tail).divided_by(gap)
            fimag = (Real.exact(1,bits)-e.real.expanded(tail)).divided_by(gap)
        else:
            freal,fimag = Real.exact(time,bits),Real.exact(0,bits)
        absolute = (freal.square()+fimag.square()).times(k1[j][hf]**2)
        variance = variance+absolute
        components.append(dict(occupation=data['sector'][j],gap=str(gap),
                               imaginary_K1=str(k1[j][hf]),squared_norm=absolute.witness()))
    a = variance.sqrt()
    c = time*entry_l1(k2)+time**2*entry_l1(v)*(entry_l1(k1)+entry_l1(k2))
    return a,c,dict(a=a.witness(),C=str(c),variance=variance.witness(),components=components,
                    normV_upper=str(entry_l1(v)),normK1_upper=str(entry_l1(k1)),
                    normK2_upper=str(entry_l1(k2)),norm_method='sum_of_absolute_entries')


def response_guard(a, c, epsilon):
    interval = a.times(epsilon).expanded(c*epsilon**2)
    return Real(max(0,interval.lo),interval.hi,interval.bits)


def infidelity(reference, reference_radius, state, state_radius):
    overlap = Complex.exact(bits=state[0].bits)
    for x,y in zip(reference,state):
        overlap = overlap+x.conjugate()*y
    polynomial = Real.exact(1,state[0].bits)-overlap.abs_square()
    eta = combine(reference_radius,state_radius,state[0].bits)
    pad = eta*(2+eta)
    expanded = polynomial.expanded(pad)
    scale = 1 << expanded.bits
    lower,upper = max(0,expanded.lo),min(scale,expanded.hi)
    require(lower <= upper, 'Unitarity and error enclosure are inconsistent')
    return Real(lower,upper,expanded.bits),dict(overlap=overlap.witness(),
        polynomial_infidelity=polynomial.witness(),analytic_overlap_radius=str(eta),
        analytic_infidelity_radius=str(pad))
