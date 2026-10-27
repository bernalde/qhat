"""Exact Fraction containment and independent high-precision series checks."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pytest
from intervals import Real, Complex, exp_i


@pytest.mark.parametrize('a,b', [(F(-7, 3), F(11, 5)), (F(2, 7), F(1, 13)),
                                (F(-1, 17), F(-19, 23)), (F(0), F(-1, 19))])
def test_exact_rational_point_arithmetic_is_enclosed(a, b):
    x, y = Real.exact(a, 32), Real.exact(b, 32)
    assert x.contains(a) and y.contains(b)
    assert (x + y).contains(a + b)
    assert (x - y).contains(a - b)
    assert (x * y).contains(a * b)
    assert x.times(b).contains(a * b)
    if b:
        assert x.divided_by(b).contains(a / b)


def test_interval_extrema_include_all_sign_combinations():
    x, y = Real(-7, 11, 16), Real(-19, 23, 16)
    result = x * y
    for a in [F(-7, 1 << 16), F(0), F(11, 1 << 16)]:
        for b in [F(-19, 1 << 16), F(0), F(23, 1 << 16)]:
            assert result.contains(a * b)
    assert x.square().lo == 0
    assert x.square().contains(F(121, 1 << 32))


def test_square_root_has_exact_endpoint_inequalities():
    x = Real.exact(F(2, 3), 96)
    root = x.sqrt()
    low, high = root.bounds()
    assert low**2 <= F(2, 3) <= high**2
    assert (root.hi - root.lo) <= 2
    with pytest.raises(ValueError, match='Nonnegative'):
        Real(-1, 2, 32).sqrt()


def test_complex_product_matches_rational_algebra():
    x = Complex.exact(F(2, 3), F(-5, 7), 64)
    y = Complex.exact(F(-11, 13), F(17, 19), 64)
    z = x * y
    assert z.real.contains(F(2, 3) * F(-11, 13) - F(-5, 7) * F(17, 19))
    assert z.imag.contains(F(2, 3) * F(17, 19) + F(-5, 7) * F(-11, 13))
    assert x.abs_square().contains(F(4, 9) + F(25, 49))
    assert (x * x.conjugate()).imag.contains(0)


@pytest.mark.parametrize('angle', [F(0), F(3, 5), F(-7, 3), F(19, 7)])
def test_scaled_exponential_matches_independent_decimal_series(angle):
    value, analytic, details = exp_i(angle)
    with localcontext() as context:
        context.prec = 110
        x = Decimal(angle.numerator) / Decimal(angle.denominator)
        cosine = ct = Decimal(1)
        sine = st = x
        for k in range(1, 180):
            ct *= -x * x / Decimal((2 * k - 1) * (2 * k))
            st *= -x * x / Decimal((2 * k) * (2 * k + 1))
            cosine += ct
            sine += st
    # The110-digit reference error is far below the separately added analytic
    # tail; the zero-angle case is exact on both paths.
    assert value.real.expanded(analytic).contains(F(cosine))
    assert value.imag.expanded(analytic).contains(F(sine))
    assert abs(angle / details['scaling']) <= F(1, 2)


def test_invalid_precision_and_unbounded_scaling_are_refused():
    with pytest.raises(ValueError, match='not a float'):
        Real.exact(.1)
    with pytest.raises(ValueError, match='Mixed'):
        Real.exact(F(1, 3), 32) + Real.exact(F(1, 3), 64)
    with pytest.raises(ValueError, match='Nonzero'):
        Real.exact(1).divided_by(0)
    with pytest.raises(ValueError, match='scaling ceiling'):
        exp_i(F(10000))


@pytest.mark.parametrize('bits', [-1, 15, 1025, 10**9, True])
def test_precision_is_checked_before_allocating_scale(bits):
    with pytest.raises(ValueError, match='precision'):
        Real.exact(1, bits)


def test_nonzero_analytic_tail_is_separate_from_polynomial_enclosure():
    value, analytic, _ = exp_i(F(1, 2), terms=2, bits=64)
    assert value.real.lo == value.real.hi == 1 << 64
    assert value.imag.lo == value.imag.hi == 1 << 63
    assert analytic == F(1, 2)
    # sin(1/2) is strictly below1/2; a polynomial interval alone is not an
    # exponential enclosure. The analytic allowance must remain in the caller.
    sine_low, sine_high = F(1, 2) - F(1, 48), F(1, 2) - F(1, 48) + F(1, 3840)
    assert not value.imag.contains(sine_high)
    assert value.imag.expanded(analytic).contains(sine_low)
    assert value.imag.expanded(analytic).contains(sine_high)
