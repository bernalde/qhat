"""Independent Fraction interval and generic-matrix certificate recheck.

Does not import producer arithmetic, unit actions, matrix routines or archived
Python. Uses160-bit rational intervals and a factorial Taylor remainder, so its
rigorous true-state/error enclosure can fit inside the producer's looser tail.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import factorial, isqrt
from pathlib import Path, PurePosixPath
import subprocess
import time

PRECISION = 160
SCALE = 1 << PRECISION
MODEL_SHA = 'ea38b9730eed2faf55fbd1fe22722d93b858b95fe4d804ecb72855d10176c134'
TIME = F(3,5)
TERMS = 96
ROOT = Path(__file__).resolve().parents[2]


def insist(ok, message):
    if not ok:
        raise ValueError(message)


def lower(x):
    return F(x.numerator*SCALE//x.denominator,SCALE)


def upper(x):
    return -lower(-x)


@dataclass(frozen=True, slots=True)
class Q:
    lo: F
    hi: F

    def __post_init__(self):
        insist(type(self.lo) is F and type(self.hi) is F and self.lo <= self.hi, 'Invalid rational interval')

    @classmethod
    def point(cls, x=0):
        insist(type(x) in (int,F), 'Exact rational input required')
        x = F(x)
        return cls(lower(x),upper(x))

    def __add__(self, other):
        return Q(lower(self.lo+other.lo),upper(self.hi+other.hi))

    def __neg__(self):
        return Q(-self.hi,-self.lo)

    def __sub__(self, other):
        return self+(-other)

    def __mul__(self, other):
        candidates = [a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Q(lower(min(candidates)),upper(max(candidates)))

    def scalar(self, x):
        x = F(x)
        a,b = (self.lo*x,self.hi*x) if x >= 0 else (self.hi*x,self.lo*x)
        return Q(lower(a),upper(b))

    def square(self):
        a = F(0) if self.lo <= 0 <= self.hi else min(self.lo**2,self.hi**2)
        return Q(lower(a),upper(max(self.lo**2,self.hi**2)))

    def sqrt(self):
        insist(self.lo >= 0,'Nonnegative square root interval')
        lo = isqrt(self.lo.numerator*SCALE*SCALE//self.lo.denominator)
        hi = isqrt(self.hi.numerator*SCALE*SCALE//self.hi.denominator)
        if F(hi*hi,SCALE*SCALE) < self.hi:
            hi += 1
        return Q(F(lo,SCALE),F(hi,SCALE))

    def expanded(self, radius):
        insist(radius >= 0,'Nonnegative analytic allowance')
        return Q(lower(self.lo-radius),upper(self.hi+radius))

    def witness(self):
        return dict(lo=str(self.lo),hi=str(self.hi),precision_bits=PRECISION)


# A complex rectangle is a pair(real_interval,imaginary_interval).
def zpoint(real=0, imag=0):
    return Q.point(real),Q.point(imag)


def zadd(a, b):
    return a[0]+b[0],a[1]+b[1]


def zmul(a, b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def zscale(a, x):
    return a[0].scalar(x),a[1].scalar(x)


def ziszero(a):
    return a[0].lo == a[0].hi == a[1].lo == a[1].hi == 0


def zidentity(n):
    return [[zpoint(int(i==j)) for j in range(n)] for i in range(n)]


def zmatrix_product(a, b):
    n = len(a)
    c = [[zpoint() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if not ziszero(a[i][k]):
                for j in range(n):
                    if not ziszero(b[k][j]):
                        c[i][j] = zadd(c[i][j],zmul(a[i][k],b[k][j]))
    return c


def zapply(a, vector):
    c = []
    for row in a:
        out = zpoint()
        for x,y in zip(row,vector):
            if not ziszero(x) and not ziszero(y):
                out = zadd(out,zmul(x,y))
        c.append(out)
    return c


def error_product(a, b):
    return upper(a+b+a*b)


def generic_exponential(h, tau):
    n = len(h)
    insist(0 < n <= 36 and all(len(row)==n for row in h),'Generic exponential size bound')
    insist(all(h[i][j] == h[j][i] for i in range(n) for j in range(n)), 'Symmetric H required')
    bound = abs(tau)*sum((abs(x) for row in h for x in row),F(0))
    repeat = 1
    while bound > F(repeat,2):
        repeat *= 2
        insist(repeat <= 4096,'Independent scaling ceiling')
    if bound == 0:
        return zidentity(n),F(0)
    generator = [[zpoint(0,-tau*x/repeat) for x in row] for row in h]
    term = zidentity(n)
    total = zidentity(n)
    for degree in range(1,TERMS):
        term = [[zscale(x,F(1,degree)) for x in row] for row in zmatrix_product(generator,term)]
        total = [[zadd(a,b) for a,b in zip(row,other)] for row,other in zip(total,term)]
    # Exact tail: the first omitted term is <=(1/2)^N/N!, and every
    # subsequent term ratio is <=1/(2(N+1)). This is tighter than2^-95.
    tail = F(1,(1 << TERMS)*factorial(TERMS))/(1-F(1,2*(TERMS+1)))
    radius = F(0)
    # Squaring of the interval polynomial, with a separate error recurrence.
    powers = repeat
    factor = total
    factor_radius = upper(tail)
    result = zidentity(n)
    while powers:
        if powers % 2:
            result = zmatrix_product(result,factor)
            radius = error_product(radius,factor_radius)
        powers //= 2
        if powers:
            factor = zmatrix_product(factor,factor)
            factor_radius = error_product(factor_radius,factor_radius)
    return result,radius


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rzero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def radd(a, b):
    return [[x+y for x,y in zip(row,other)] for row,other in zip(a,b)]


def rscale(a, scalar):
    return [[x*scalar for x in row] for row in a]


def rcommutator(a, b):
    n = len(a)
    return [[sum((a[i][k]*b[k][j]-b[i][k]*a[k][j] for k in range(n)),F(0))
             for j in range(n)] for i in range(n)]


def norm_upper(a):
    return sum((abs(x) for row in a for x in row),F(0))


def literal_matrices(model, scalar):
    """Local X/Y/Z actions over all64 columns; no producer mask descriptors."""
    insist((model['n'],model['electrons']) == (6,2),'Fixed independent molecular dimensions')
    n = model['n']
    full_n = 1 << n
    den = int(model['denominator'])
    insist(den > 0 and den.bit_length() <= 512 and len(model['blocks']) <= 128,'Exact input bounds')
    count = lambda v:tuple(sum((v >> q)&1 for q in range(parity,n,2)) for parity in (0,1))
    hf = (1 << model['electrons'])-1
    sector = [v for v in range(full_n) if count(v) == count(hf)]
    insist(len(sector) == 9,'Declared LiH sector dimension')
    matrices,all_terms,diagonal_ids = [],{},[]
    for index,block in enumerate(model['blocks']):
        insist(0 < len(block) <= 128,'Unit term bound')
        real = [[0 for _ in range(full_n)] for _ in range(full_n)]
        imag = [[0 for _ in range(full_n)] for _ in range(full_n)]
        is_diagonal = True
        for word,weight in block:
            key = tuple((q,p) for q,p in word)
            insist(key and key not in all_terms and type(weight) is str and str(int(weight)) == weight,
                   'Exact disjoint nonidentity input terms required')
            c = int(weight)
            insist(c and c.bit_length() <= 512,'Coefficient size bound')
            all_terms[key] = c
            for column in range(full_n):
                row,re,im = column,1,0
                for q,letter in word:
                    insist(type(q) is int and 0 <= q < n and letter in ('X','Y','Z'),'Local Pauli domain')
                    bit = (column >> q)&1
                    if letter == 'X':
                        row ^= 1 << q
                        is_diagonal = False
                    elif letter == 'Y':
                        row ^= 1 << q
                        re,im = (-im,re) if bit == 0 else (im,-re)
                        is_diagonal = False
                    elif bit:
                        re,im = -re,-im
                real[row][column] += c*re
                imag[row][column] += c*im
        for i in range(full_n):
            for j in range(full_n):
                insist(not imag[i][j] and real[i][j] == real[j][i],'Real Hermitian full-space unit')
                insist(not real[i][j] or count(i) == count(j),'Full-space sector conservation failed')
        matrices.append([[F(real[i][j],den)*scalar for j in sector] for i in sector])
        if is_diagonal:
            diagonal_ids.append(index)
    total = {tuple((q,p) for q,p in word):int(c) for word,c in model['coefficients']}
    insist(all_terms == total and len(total) == len(model['coefficients']),'Independent full unit partition check')
    return matrices,sector,sector.index(hf),diagonal_ids


def independent_bch(matrices, order, diagonal):
    n = len(matrices[0])
    k1,k2 = rzero(n),rzero(n)
    for i in range(len(order)):
        for j in range(i+1,len(order)):
            left,right = order[i],order[j]
            term = rscale(rcommutator(matrices[left],matrices[right]),F(-1,2))
            if left in diagonal and right in diagonal:
                insist(term == rzero(n),'Diagonal unit commutator must vanish')
            elif left in diagonal or right in diagonal:
                k1 = radd(k1,term)
            else:
                k2 = radd(k2,term)
    return k1,k2


def slope_and_constant(matrices, diagonal, hf, k1, k2):
    n = len(matrices[0])
    d,v = rzero(n),rzero(n)
    for i,matrix in enumerate(matrices):
        if i in diagonal:
            d = radd(d,matrix)
        else:
            v = radd(v,matrix)
    squared = Q.point()
    for j in range(n):
        if j == hf or not k1[j][hf]:
            continue
        gap = d[j][j]-d[hf][hf]
        if gap:
            exponential,tail = generic_exponential([[-gap]],TIME)
            real,imag = [x.expanded(tail) for x in exponential[0][0]]
            freal,fimag = imag.scalar(1/gap),(Q.point(1)-real).scalar(1/gap)
        else:
            freal,fimag = Q.point(TIME),Q.point()
        squared = squared+(freal.square()+fimag.square()).scalar(k1[j][hf]**2)
    constant = TIME*norm_upper(k2)+TIME**2*norm_upper(v)*(norm_upper(k1)+norm_upper(k2))
    return squared.sqrt(),constant


def producer_interval(raw):
    bits = raw['bits']
    insist(type(bits) is int and bits == 192,'Producer precision mismatch')
    lo,hi = int(raw['lo']),int(raw['hi'])
    insist(lo <= hi and max(abs(lo).bit_length(),abs(hi).bit_length()) <= 2048,'Producer interval domain')
    return F(lo,1 << bits),F(hi,1 << bits)


def contained(independent, raw, radius=F(0)):
    lo,hi = producer_interval(raw)
    insist(lo-radius <= independent.lo <= independent.hi <= hi+radius,
           'Independent rigorous enclosure is outside producer certificate')


def check_vector(vector, tail, raw, producer_tail):
    insist(len(vector) == len(raw),'State dimensions differ')
    for z,record in zip(vector,raw):
        for part,name in zip(z,('real','imag')):
            contained(part.expanded(tail),record[name],producer_tail)


def squared_error(reference, rt, state, st):
    overlap = zpoint()
    for a,b in zip(reference,state):
        overlap = zadd(overlap,zmul((a[0],-a[1]),b))
    value = Q.point(1)-overlap[0].square()-overlap[1].square()
    eta = error_product(rt,st)
    value = value.expanded(eta*(2+eta))
    return Q(max(F(0),value.lo),min(F(1),value.hi))


def verify(output, destination):
    insist(not destination.exists(),'Retain original verifier output')
    cpu,wall = time.process_time(),time.monotonic()
    frozen = read(output/'FROZEN.json')
    complete = read(output/'PRODUCER_COMPLETE.json')
    insist(digest(output/'FROZEN.json') == complete['freeze_sha256'],'Freeze digest changed')
    insist(digest(output/'model.json') == frozen['model_sha256'] == MODEL_SHA,'Fixed input digest changed')
    head = frozen['source_commit']
    insist(type(head) is str and len(head) == 40 and all(c in '0123456789abcdef' for c in head),'Invalid producer Git identity')
    listed = subprocess.check_output(['git','-C',str(ROOT),'ls-tree','-r','--name-only',head,'--',
                                     'numerical_audit/molecular_guard'],text=True).splitlines()
    expected_sources = {name for name in listed if Path(name).suffix in ('.py','.md')}
    expected_sources.update('numerical_audit/prospective_reference/'+name for name in ('engine.py','campaign.py','one_pass.py'))
    insist(set(frozen['source_sha256']) == expected_sources,'Incomplete producer source inventory')
    for name,expected in frozen['source_sha256'].items():
        path = PurePosixPath(name)
        insist(not path.is_absolute() and '..' not in path.parts and str(path) == name,'Unsafe archived source path')
        archived = output/'source_snapshot'/name
        insist(digest(archived) == expected,'Archived source digest changed')
        blob = subprocess.check_output(['git','-C',str(ROOT),'show',frozen['source_commit']+':'+name])
        insist(blob == archived.read_bytes(),'Archived producer source differs from its Git commit')
    preflight = frozen['mode'] == 'preflight'
    insist(frozen['mode'] in ('preflight','main'),'Unknown frozen panel')
    epsilons = [F(1),F(1,256)] if preflight else [F(1,1 << k) for k in range(0,21,2)]
    steps_panel = [2] if preflight else [1,2,4,8,16,32]
    scalar = F(7,5) if preflight else F(1)
    insist(frozen['epsilons'] == list(map(str,epsilons)) and frozen['steps'] == steps_panel and
           frozen['coefficient_scale'] == str(scalar) and frozen['time'] == str(TIME) and
           (frozen['bits'],frozen['terms']) == (192,TERMS),'Predeclared panel changed')
    matrices,sector,hf,diagonal = literal_matrices(read(output/'model.json'),scalar)
    insist(sector == frozen['sector'] and frozen['hf'] == sector[hf] and
           frozen['exact_unit_matrices'] == [[[str(x) for x in row] for row in a] for a in matrices],
           'Independent complete-unit reconstruction differs')
    insist(set(frozen['orders']) == {'current','diagonal_prefix'},'Wrong fixed schedule pair')
    for order in frozen['orders'].values():
        insist(sorted(order) == list(range(len(matrices))),'Invalid complete permutation')
    insist(frozen['orders']['current'] == frozen['legacy_current']['order'],'Frozen legacy order differs')
    model = read(output/'model.json')
    blocks = model['blocks']
    diagonal_priority = sorted(diagonal,key=lambda i:(-abs(int(blocks[i][0][1])),tuple(map(tuple,blocks[i][0][0]))))
    flip = lambda i:sum(1 << q for q,p in blocks[i][0][0] if p in ('X','Y'))
    off_priority = sorted((i for i in range(len(blocks)) if i not in diagonal),
                          key=lambda i:(-sum(abs(int(c)) for word,c in blocks[i]),flip(i)))
    insist(frozen['diagonal_ids'] == diagonal_priority and frozen['off_diagonal_ids'] == off_priority and
           frozen['orders']['diagonal_prefix'] == diagonal_priority+off_priority,'Frozen priorities/prefix differ')
    insist(frozen['full_space_invariance_columns_checked'] == 64*len(matrices),'Incorrect invariance scope')
    leading = read(output/'leading_response.json')
    verified_slopes = {}
    for label,order in frozen['orders'].items():
        k1,k2 = independent_bch(matrices,order,diagonal)
        record = leading[label]
        insist(record['imaginary_K1'] == [[str(x) for x in row] for row in k1] and
               record['imaginary_K2'] == [[str(x) for x in row] for row in k2],'BCH decomposition differs')
        a,c = slope_and_constant(matrices,diagonal,hf,k1,k2)
        contained(a,record['a'])
        insist(c == F(record['C']),'Response remainder constant differs')
        verified_slopes[label] = a,c
    initial = [zpoint(int(j==hf)) for j in range(len(sector))]
    rows,comparisons,guards = [],[],[]
    for at,epsilon in enumerate(epsilons):
        actual = [rscale(h,1 if i in diagonal else epsilon) for i,h in enumerate(matrices)]
        h = rzero(len(sector))
        for a in actual:
            h = radd(h,a)
        reference_operator,rt = generic_exponential(h,TIME)
        reference = zapply(reference_operator,initial)
        rawref = read(output/f'references/epsilon-{at:02d}.json')
        insist(rawref['epsilon'] == str(epsilon),'Reference identity differs')
        check_vector(reference,rt,rawref['state'],F(rawref['analytic_radius']))
        guard_record = complete['response_guards'][at]
        insist(guard_record['epsilon'] == str(epsilon),'Guard identity differs')
        guard_intervals = {}
        for label,(slope,c) in verified_slopes.items():
            a = slope.scalar(epsilon).expanded(c*epsilon**2)
            guard_intervals[label] = Q(max(F(0),a.lo),a.hi)
            contained(guard_intervals[label],guard_record['intervals'][label])
        ga,gb = [producer_interval(guard_record['intervals'][key]) for key in ('current','diagonal_prefix')]
        guard_winner = 'current' if ga[1] < gb[0] else ('diagonal_prefix' if gb[1] < ga[0] else 'unresolved')
        insist(guard_record['strict_smaller_response'] == guard_winner,'Guard decision differs from intervals')
        guards.append(dict(epsilon=str(epsilon),verified_producer_decision=guard_winner,
                           independent_response_intervals={k:v.witness() for k,v in guard_intervals.items()}))
        for steps in steps_panel:
            operators = [generic_exponential(a,TIME/steps) for a in actual]
            errors = {}
            for label,order in frozen['orders'].items():
                state,st = initial,F(0)
                for _ in range(steps):
                    for index in reversed(order):
                        operator,tail = operators[index]
                        state = zapply(operator,state)
                        st = error_product(st,tail)
                name = f'rows/epsilon-{at:02d}-r{steps:02d}-{label}.json'
                raw = read(output/name)
                insist((raw['epsilon'],raw['steps'],raw['order']) == (str(epsilon),steps,label),'Circuit identity differs')
                check_vector(state,st,raw['state'],F(raw['state_analytic_radius']))
                error = squared_error(reference,rt,state,st)
                contained(error,raw['squared_infidelity'])
                errors[label] = error
                rows.append(dict(path=name,independent_squared_infidelity=error.witness()))
            a,b = errors['current'],errors['diagonal_prefix']
            winner = 'current' if a.hi < b.lo else ('diagonal_prefix' if b.hi < a.lo else 'unresolved')
            producer = next(x for x in complete['finite_rankings'] if x['epsilon'] == str(epsilon) and x['steps'] == steps)
            pa,pb = [producer_interval(read(output/f'rows/epsilon-{at:02d}-r{steps:02d}-{key}.json')['squared_infidelity'])
                     for key in ('current','diagonal_prefix')]
            producer_winner = 'current' if pa[1] < pb[0] else ('diagonal_prefix' if pb[1] < pa[0] else 'unresolved')
            insist(producer['strict_smaller_error'] == producer_winner and producer_winner in ('unresolved',winner),
                   'Producer strict ranking contradicted')
            comparisons.append(dict(epsilon=str(epsilon),steps=steps,independent_strict_smaller_error=winner))
            print('INDEPENDENTLY_ENCLOSED',str(epsilon),'r',steps,winner,flush=True)
    expected = len(epsilons)*len(steps_panel)*2
    insist(len(rows) == complete['circuits'] == expected and
           len(epsilons) == complete['reference_evolutions'] and
           complete['rows'] == [x['path'] for x in rows],'Incomplete declared panel')
    result = dict(status='passed',circuits=len(rows),reference_evolutions=len(epsilons),
        producer_freeze_sha256=digest(output/'FROZEN.json'),verifier_sha256=digest(Path(__file__)),
        precision_bits=PRECISION,method='Independent literal local Pauli actions,Fraction intervals,generic matrix exponentials,factorial tail',
        rows=rows,response_guards=guards,finite_rankings=comparisons,CPU_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,
        scope='Software interval certificate;no archived Python executed;no concrete Lean instantiation')
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    args = parser.parse_args()
    verify(args.output,args.destination)
