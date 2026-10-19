"""Tests for common/hermiticity.py."""

import numpy as np
import pytest

from qhat.common.hermiticity import (
    check_hermitian_one_body,
    check_hermitian_two_body,
    max_relative_imaginary,
    real_pauli_coefficients,
    require_real_scalar,
)


def random_hermitian_one_body(n, rng):
    h = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (h + h.conj().T) / 2


def random_hermitian_two_body(n, rng):
    """Two-body tensor of a Hermitian operator, symmetrized so h[p,q,r,s] = conj(h[s,r,q,p])."""
    h = rng.normal(size=(n,) * 4) + 1j * rng.normal(size=(n,) * 4)
    return (h + h.transpose(3, 2, 1, 0).conj()) / 2


class TestRequireRealScalar:

    @pytest.mark.parametrize("x", [1.5, -2, np.float64(3.0), np.int64(4), 2.0 + 0j,
                                   np.complex128(-5.0 + 1e-19j), np.array(6.0)])
    def test_accepts_real_values(self, x):
        result = require_real_scalar("x", x)
        assert type(result) is float
        assert result == pytest.approx(float(np.real(x)))

    @pytest.mark.parametrize("x", [1.0 + 1e-6j, np.complex128(0.5j), 1j])
    def test_rejects_complex_values(self, x):
        with pytest.raises(ValueError, match="x must be real"):
            require_real_scalar("x", x)

    @pytest.mark.parametrize("x", [np.inf, -np.inf, np.nan, complex(np.nan, 0)])
    def test_rejects_non_finite_values(self, x):
        with pytest.raises(ValueError, match="must be finite"):
            require_real_scalar("x", x)

    @pytest.mark.parametrize("x", ["1.0", None, True, [1.0]])
    def test_rejects_non_numbers(self, x):
        with pytest.raises(TypeError):
            require_real_scalar("x", x)

    def test_tolerance_is_relative(self):
        assert require_real_scalar("x", 1e6 + 1e-7j) == 1e6
        with pytest.raises(ValueError):
            require_real_scalar("x", 1e-6 + 1e-7j)


class TestCheckHermitianOneBody:

    def test_accepts_hermitian(self):
        check_hermitian_one_body("h1", random_hermitian_one_body(5, np.random.default_rng(0)))

    def test_accepts_real_symmetric_and_zero(self):
        h = np.random.default_rng(1).normal(size=(4, 4))
        check_hermitian_one_body("h1", h + h.T)
        check_hermitian_one_body("h1", np.zeros((3, 3)))

    def test_rejects_non_hermitian(self):
        h = random_hermitian_one_body(4, np.random.default_rng(2))
        h[1, 3] += 0.1j
        with pytest.raises(ValueError, match=r"h1 is not Hermitian.*index \((1, 3|3, 1)\)"):
            check_hermitian_one_body("h1", h)

    def test_rejects_round_off_above_tolerance(self):
        h = random_hermitian_one_body(4, np.random.default_rng(3))
        h[0, 1] += 1e-9
        with pytest.raises(ValueError):
            check_hermitian_one_body("h1", h)
        check_hermitian_one_body("h1", h, rtol=1e-6)

    def test_rejects_non_finite(self):
        h = np.eye(2)
        h[0, 0] = np.nan
        with pytest.raises(ValueError, match="non-finite"):
            check_hermitian_one_body("h1", h)


class TestCheckHermitianTwoBody:

    def test_accepts_hermitian(self):
        check_hermitian_two_body("h2", random_hermitian_two_body(3, np.random.default_rng(4)))

    def test_accepts_hermitian_operator_with_unsymmetrized_tensor(self):
        """Only the part antisymmetric in (p, q) and (r, s) defines the operator."""
        n = 3
        rng = np.random.default_rng(5)
        h = random_hermitian_two_body(n, rng)
        # A term symmetric under p <-> q contributes nothing to a+_p a+_q a_r a_s.
        junk = rng.normal(size=(n,) * 4) + 1j * rng.normal(size=(n,) * 4)
        junk = junk + junk.transpose(1, 0, 2, 3)
        unsymmetrized = h + junk
        assert np.abs(unsymmetrized - unsymmetrized.transpose(3, 2, 1, 0).conj()).max() > 0.1
        check_hermitian_two_body("h2", unsymmetrized)

    def test_rejects_non_hermitian(self):
        h = random_hermitian_two_body(3, np.random.default_rng(6))
        h[0, 1, 2, 0] += 0.05j
        with pytest.raises(ValueError, match="h2 is not Hermitian"):
            check_hermitian_two_body("h2", h)


class TestMaxRelativeImaginary:

    def test_values(self):
        assert max_relative_imaginary(np.ones(3)) == 0.0
        assert max_relative_imaginary(np.zeros(3, dtype=complex)) == 0.0
        assert max_relative_imaginary(np.array([2.0, 0.5j])) == pytest.approx(0.25)


class TestRealPauliCoefficients:

    def test_returns_floats(self):
        terms = {(): np.complex128(-1e4 + 1e-19j), ((0, 'Z'),): 0.5 + 0j, ((1, 'X'),): -0.25}
        result = real_pauli_coefficients(terms)
        assert result == {(): -1e4, ((0, 'Z'),): 0.5, ((1, 'X'),): -0.25}
        assert all(type(c) is float for c in result.values())

    def test_rejects_imaginary_parts(self):
        terms = {((0, 'Z'),): 0.5 + 0.01j, ((1, 'X'),): -0.25}
        with pytest.raises(ValueError, match="not real.*not Hermitian"):
            real_pauli_coefficients(terms)

    def test_large_identity_does_not_mask_imaginary_parts(self):
        """Scaling by the identity term would hide an imaginary part of 1e-9 here."""
        terms = {(): 1e5, ((0, 'Z'),): 0.5 + 1e-9j, ((1, 'X'),): -0.25}
        with pytest.raises(ValueError):
            real_pauli_coefficients(terms)

    def test_dense_identity_key(self):
        terms = {"II": 1e5, "ZI": 0.5 + 1e-9j}
        with pytest.raises(ValueError):
            real_pauli_coefficients(terms)

    def test_identity_only(self):
        assert real_pauli_coefficients({(): 2.0 + 0j}) == {(): 2.0}
        with pytest.raises(ValueError):
            real_pauli_coefficients({(): 2.0 + 1j})

    def test_rejects_non_finite(self):
        with pytest.raises(ValueError, match="not finite"):
            real_pauli_coefficients({((0, 'Z'),): np.nan})
