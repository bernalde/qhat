"""Deterministic rational extended-Hubbard scaling family in a fixed basis."""
from fractions import Fraction
from functools import reduce
import math


def synthetic(n):
    if n not in (4,6,8,10,12):raise ValueError('Declared scaling range is 4,6,8,10,12 modes')
    coefficients={};off=[];scalar=Fraction()
    def add(word,c):coefficients[word]=coefficients.get(word,Fraction())+c
    # epsilon_p n_p; omit its scalar consistently.
    for p in range(n):
        epsilon=Fraction(2*(p//2+1)+(p%2),20)
        scalar+=epsilon/2
        add(((p,'Z'),),-epsilon/2)
    # Spin-preserving hopping, complete XX/YY partners with exact JW strings.
    for p in range(n):
        for q in range(p+1,n):
            if p%2!=q%2:continue
            t=-Fraction(3,10*(1+(q-p)//2))
            middle=tuple((j,'Z') for j in range(p+1,q))
            block={((p,'X'),)+middle+((q,'X'),):t/2,
                   ((p,'Y'),)+middle+((q,'Y'),):t/2}
            off.append(block)
            for w,c in block.items():add(w,c)
    # Long-range density coupling, with a stronger on-site alpha/beta term.
    for p in range(n):
        for q in range(p+1,n):
            U=Fraction(3,5) if p//2==q//2 else Fraction(1,10*(1+(q-p)))
            scalar+=U/4
            add(((p,'Z'),),-U/4);add(((q,'Z'),),-U/4)
            add(((p,'Z'),(q,'Z')),U/4)
    coefficients={w:c for w,c in coefficients.items() if c}
    den=reduce(math.lcm,(c.denominator for c in coefficients.values()),1)
    integer={w:int(c*den) for w,c in coefficients.items()}
    diagonals=[{w:c} for w,c in sorted(integer.items()) if all(p=='Z' for q,p in w)]
    blocks=diagonals+[{w:int(c*den) for w,c in b.items()} for b in off]
    encode=lambda b:[[[[q,p] for q,p in w],str(c)] for w,c in sorted(b.items())]
    return dict(case='extended_hubbard_'+str(n),family='synthetic extended Hubbard',n=n,
        electrons=n//2,denominator=str(den),coefficients=encode(integer),blocks=[encode(b) for b in blocks],
        coefficient_policy='exact rational deterministic model; scalar identity removed',
        removed_scalar_exact=str(scalar),
        formula='epsilon_p=(2*(floor(p/2)+1)+(p%2))/20; t_pq=-3/[10*(1+(q-p)/2)] for equal parity; U_pq=3/5 on-site, else 1/[10*(1+q-p)]')
