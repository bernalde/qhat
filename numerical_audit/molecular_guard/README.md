# Fixed molecular guard and finite-step enclosures

Status: the complete frozen LiH panel now has independently verified software
interval certificates for all 132 circuits and 11 reference evolutions. The
current order has strictly lower squared infidelity in all 66 tested pairs,
including six at physical epsilon=1. The leading-response guard remains
unresolved at epsilon=1, 1/4 and 1/16; it is strict on the predeclared scaled
grid from 1/64 downward. See [the complete report](evidence/report/README.md).

The protocol and both orders were frozen before errors. All payload bytes,
archive extraction and staged Git blobs were checked. Input provenance is
explicitly developmental until the actual issue #17 merge. Paired Astra/Argon
source reviews, publication and merge CI remain pending. The separate generic Lean
modules passed a fresh full-root build of 80 modules and audit of all 2,527
project constants, including 20 new public endpoints. Their exact source and
hypotheses are retained with the paper proof inventory. Independent semantic
assessment remains required; these generic proofs do not instantiate the molecular certificate.

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
The independent generic rational interval checker and complete delivery/replay
validation have passed. Paired source reviews and actual publication/merge remain
required. See[PROTOCOL.md](PROTOCOL.md) for the exact
input,panel,resource ceilings and physical-versus-scaled-Hamiltonian scope.

Tests compare integer Pauli actions with literal local matrices, BCH signs with
hand-computable commutators, matrix polynomials with rational coefficients,
chronological conventions with independent SciPy propagation, and scalar
intervals with independent high-precision series. These software tests do not
establish a Lean molecular proof.


Report replays may execute at a later commit. The reporter authenticates both
verifier files against their executed Git revision and the maintained source,
then requires those source hashes to agree across reports. Only execution
revision fields may differ: all scientific tables, producer identity, reporting
source hash, resource snapshot and exact molecular input merge stay bound.
The original development report is retained when the merged input pin is added.


The numerical input is now pinned to verified QHAT PR #8 merge
`1eeae416db1cac0de3e5617ca30ecf95729df04d`, with identical supplied model bytes.
`evidence/HISTORICAL_DEVELOPMENT_STORAGE.json` authenticates all 192 original
evidence files, including inert copies of the earlier development report and
provenance. The maintained report uses the actual merged input and unchanged
scientific CSVs. `FINAL_INPUT_BINDING.json` distinguishes this reporting refresh
from the original physical computation and independent interval audit.


Fresh freezes authenticate the direct maintained Python/Markdown source files
and the three explicit legacy helpers. Nested historical snapshots and report
README files remain data and cannot enlarge that producer-source inventory.
The verifier still rejects omitted maintained sources and authenticates each
frozen source byte against the producing Git revision. Source-contract fixtures
stop before molecular reconstruction; fresh numerical replay is recorded
separately when the executed verifier revision changes.
