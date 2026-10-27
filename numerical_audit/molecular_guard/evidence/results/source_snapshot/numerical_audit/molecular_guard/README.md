# Fixed molecular guard and finite-step enclosures

Status: implementation in progress. The protocol was committed before any guard
or physical error. No independently audited molecular certificate or new Lean
kernel result is available yet.

`exact.py` checks every complete unit on every full-space occupation column
before restricting to a conserved spin sector. It reconstructs exact rational
unit matrices and the physical BCH coefficients
`K(epsilon)=epsilon*K1+epsilon^2*K2`, with `K=-i/2 sum[i<j][Gi,Gj]`.
The written product applies its units in reverse chronological order.

`intervals.py` implements outward dyadic rectangles with integer arithmetic;
`propagation.py` encloses exact exponential Taylor polynomials and separately
tracks their analytic truncation errors. Scaling places each generator below
norm1/2. The geometric96-term bound is2^-95. Products of approximate operators
whose exact operators are unitary have error at most`(1+a)(1+b)-1`; each positive
radius is rounded upward to prevent unbounded rational denominator growth.
The Euclidean operator-norm bounds use the conservative sum of absolute entries.

Complete real equal-flip units have disjoint weighted2x2 blocks, permitting
scalar polynomial evaluation. The reference full-Hamiltonian exponential uses a
generic matrix polynomial and finite powering. Arithmetic interval widths enclose
coefficient/operation rounding, while analytic radii cover the difference between
those polynomials and the exact exponentials. Approximate states are never
renormalized. Overlap uncertainty`eta` expands squared infidelity by
`eta*(2+eta)`; only then is the true unitary error interval intersected with[0,1].

The leading response uses exact diagonal energy gaps, the analytic zero-gap
limit and coherent complex transport. The bound
`C=T||K2||+T^2||V||(||K1||+||K2||)` gives the standard-deviation enclosure
`epsilon*a +/- C*epsilon^2` for0<=epsilon<=1. Strictly separated intervals are
required for rankings; shared endpoints remain unresolved.

`study.py freeze` authenticates the fixed model/source bytes, constructs both
predeclared permutations and archives them before calculating guards or errors.
`study.py evaluate` consumes that freeze and retains the entire declared panel.
Producer completion is not independent verification. The independent generic
rational interval checker, delivery/replay validation, paired reviews and actual
publication/merge remain required. See[PROTOCOL.md](PROTOCOL.md) for the exact
input,panel,resource ceilings and physical-versus-scaled-Hamiltonian scope.

Tests compare integer Pauli actions with literal local matrices, BCH signs with
hand-computable commutators, matrix polynomials with rational coefficients,
chronological conventions with independent SciPy propagation, and scalar
intervals with independent high-precision series. These software tests do not
establish a Lean molecular proof.
