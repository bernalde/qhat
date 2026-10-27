"""Outward dyadic intervals for exact-rational complex certificate calculations.

Bounds are integer multiples of2^-bits. No binary floating-point arithmetic is
used. Arithmetic rounding is enclosed here; analytic Taylor error is separate.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


def require(ok, message):
    if not ok:
        raise ValueError(message)


def ceil_div(a, b):
    require(type(b) is int and b > 0, 'Positive interval divisor required')
    return -((-a) // b)


@dataclass(frozen=True, slots=True)
class Real:
    lo: int
    hi: int
    bits: int = 192

    def __post_init__(self):
        require(type(self.lo) is int and type(self.hi) is int and self.lo <= self.hi and
                type(self.bits) is int and 16 <= self.bits <= 1024, 'Invalid dyadic interval')

    @property
    def scale(self):
        return 1 << self.bits

    @classmethod
    def exact(cls, value, bits=192):
        require(type(value) in (int, Fraction), 'Use an explicit exact rational, not a float')
        require(type(bits) is int and 16 <= bits <= 1024, 'Invalid dyadic precision')
        value = Fraction(value)
        scale = 1 << bits
        return cls(value.numerator * scale // value.denominator,
                   ceil_div(value.numerator * scale, value.denominator), bits)

    def match(self, other):
        require(isinstance(other, Real) and self.bits == other.bits, 'Mixed interval precisions')

    def __add__(self, other):
        self.match(other)
        return Real(self.lo + other.lo, self.hi + other.hi, self.bits)

    def __neg__(self):
        return Real(-self.hi, -self.lo, self.bits)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        self.match(other)
        products = [self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi]
        return Real(min(products) // self.scale, ceil_div(max(products), self.scale), self.bits)

    def times(self, value):
        require(type(value) in (int, Fraction), 'Exact scalar required')
        value = Fraction(value)
        low, high = (self.lo, self.hi) if value >= 0 else (self.hi, self.lo)
        return Real(low * value.numerator // value.denominator,
                    ceil_div(high * value.numerator, value.denominator), self.bits)

    def divided_by(self, value):
        require(type(value) in (int, Fraction) and value != 0, 'Nonzero exact divisor required')
        return self.times(1 / Fraction(value))

    def square(self):
        low = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        high = max(self.lo**2, self.hi**2)
        return Real(low // self.scale, ceil_div(high, self.scale), self.bits)

    def sqrt(self):
        require(self.lo >= 0, 'Nonnegative interval required for square root')
        low = isqrt(self.lo * self.scale)
        upper_square = self.hi * self.scale
        high = isqrt(upper_square)
        if high * high < upper_square:
            high += 1
        return Real(low, high, self.bits)

    def expanded(self, radius):
        require(type(radius) in (int, Fraction) and radius >= 0, 'Nonnegative exact radius required')
        amount = ceil_div(Fraction(radius).numerator * self.scale, Fraction(radius).denominator)
        return Real(self.lo - amount, self.hi + amount, self.bits)

    def contains(self, value):
        value = Fraction(value)
        return self.lo * value.denominator <= value.numerator * self.scale <= self.hi * value.denominator

    def bounds(self):
        return Fraction(self.lo, self.scale), Fraction(self.hi, self.scale)

    def witness(self):
        return dict(lo=str(self.lo), hi=str(self.hi), bits=self.bits)


@dataclass(frozen=True, slots=True)
class Complex:
    real: Real
    imag: Real

    def __post_init__(self):
        self.real.match(self.imag)

    @property
    def bits(self):
        return self.real.bits

    @classmethod
    def exact(cls, real=0, imag=0, bits=192):
        return cls(Real.exact(real, bits), Real.exact(imag, bits))

    def __add__(self, other):
        return Complex(self.real + other.real, self.imag + other.imag)

    def __neg__(self):
        return Complex(-self.real, -self.imag)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Complex(self.real * other.real - self.imag * other.imag,
                       self.real * other.imag + self.imag * other.real)

    def times(self, value):
        return Complex(self.real.times(value), self.imag.times(value))

    def conjugate(self):
        return Complex(self.real, -self.imag)

    def abs_square(self):
        return self.real.square() + self.imag.square()

    def power(self, exponent):
        require(type(exponent) is int and 0 <= exponent <= 4096, 'Finite power bound')
        result = Complex.exact(1, bits=self.bits)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def witness(self):
        return dict(real=self.real.witness(), imag=self.imag.witness())


def exp_i(angle, *, terms=96, bits=192):
    """Polynomial enclosure for exp(i*angle), plus a separate analytic radius.

    Scaling ensures norm(generator)<=1/2. The96-term tail2^-95 is conservative;
    finite powering uses(1+delta)^s-1. The interval covers arithmetic rounding.
    """
    require(type(angle) in (int, Fraction) and type(terms) is int and 2 <= terms <= 256,
            'Bounded exact angle and Taylor term count required')
    angle = Fraction(angle)
    if not angle:
        return Complex.exact(1, bits=bits), Fraction(0), dict(terms=terms, scaling=1)
    scaling = 1
    while 2 * abs(angle) > scaling:
        scaling *= 2
        require(scaling <= 4096, 'Exponential scaling ceiling')
    generator = Complex.exact(0, angle / scaling, bits)
    term = result = Complex.exact(1, bits=bits)
    for degree in range(1, terms):
        term = (generator * term).times(Fraction(1, degree))
        result = result + term
    result = result.power(scaling)
    local = Fraction(1, 1 << (terms - 1))
    radius = (1 + local)**scaling - 1
    return result, radius, dict(terms=terms, scaling=scaling, local_analytic_radius=str(local),
                               total_analytic_radius=str(radius))
