# Exact response bounds: useful restricted certificates, full-order refusal

The new separator/conditioning bounds are useful on the two retained eight-label
windows. The archived best order has a certified objective ratio of1.00009015
for CH4 and1.00494983 for CO: at most0.0091% and0.495% above the restricted
optimum, respectively. The complete full-order projection failed its resource
preflight, so this work establishes no full-molecule approximation ratio.

| Frozen response-table domain | Retained parent / lower bound | Archived best order / lower bound |
|---|---:|---:|
| CH4, positions157:164 |1.00041838|1.00009015|
| CO, positions168:175 |1.32168549|1.00494983|

Every outside label and both endpoints stay fixed. Original word counts give
r=15 for these tables, so the same ratios apply to Q and J=Q/r^4. The bound
relaxes the coupled local ordering choices, proposes rational separators, and
checks every comparison branch with exact integer/rational arithmetic. It does
not enumerate all orders. A separately implemented arithmetic checker accepted
both certificates; the table bit digests match the archived exhaustive study,
and exact checks establish L<=archived optimum<=retained parent. The archived
best orders are reused witnesses; no fresh enumeration or new physical circuit
was run. The two cases share one cold process, so their timings do not establish
a cross-case speedup or a warm benchmark.

The target is the supplied frozen binary64 response table interpreted exactly
as dyadic rationals. These certificates do not cover exact transcendental
transport, the new normal optimizer's separately aggregated tables, different
outside orders, arbitrary endpoints, or actual Trotter error. The experiment
supports certifying a small conditioned window; it does not make the global
ordering problem efficiently approximable.

## Full-order molecular feasibility

Before examining bound quality, the frozen resource-only preflight selected
the larger retained development input, CH4:12 modes,351 original units and1722
original words. A distinct input scales every exact coefficient by7/5 and
freezes a new T=1 transport kernel. The full same-model incumbent response
completed. The complete eight-target projection then returned `scope_limit`
because its scope inventory exceeded50,000. This is an incomplete projection,
not a completed zero lower bound. Its guarded job used13.744 CPU seconds and
57.178 wall seconds, with an owned-cgroup peak around313.4MiB.

The predeclared feasibility gate therefore failed. The two unaltered full
development cases and four prospective fourteen-mode full-domain bound cases
were not run; they remain explicitly unperformed. No easier input, reduced
target set, increased scope cap or fresh budget replaced this failed attempt.
The original refusal checker recorded the input/kernel/result digests and checked
the bounded cap declaration and executable sources; it explicitly records
`scientific_verification_performed=false` and verifies no positive certificate.
It did not enforce equality of the producer's input/kernel digests with the
supplied artifacts. The paired review independently compared the retained
digests and found those particular links agree. This passive agreement must not
be confused with enforcement by that historical checker.

This failure concerns the current explicit scope representation under the
frozen limits. It is not an impossibility theorem for molecular approximation.
An implicit support/separation oracle that avoids materializing the whole
scope inventory is a concrete next research direction. Its computational cost
and useful lower-bound strength remain unestablished.

## Validation and evidence

The Julia backend passed570 assertions across10 test sets, covering exact dense
polynomial identities, sign/normalization conventions, complete projection
coverage, branch arithmetic, zero relaxations, caps and refusal receipts.
Two parser failures and one malformed-fixture failure preceded that pass; all
raw attempts and their costs are retained. The malformed fixture had converted
its intended Float64 back to an integer before validation; the fixture was fixed
without weakening cap validation.

Producer source is55cc7ee9e053ec47b9432d603c0a08b2af300d91. The execution freeze
binds all backend bytes. `analyze.py` derives the exact ratios from hash-bound
raw results and oracle checks. `evidence/ANALYSIS.a1.json` retains the exact
rational ratios and cumulative costs through the window checks:116.16 charged
CPU seconds and481.411 guarded wall seconds, including all failed attempts.
The derivation's own later cost and packaging costs are separate receipts.
These costs are well inside the unchanged2160-CPU/10800-wall campaign pool.

## Driver hardening after the retained study

The review fixes add handled input/kernel byte admission before user-artifact
hashing, enforce available producer digests during refusal replay, and label
pre-digest or missing-artifact refusals explicitly `recorded_refusal_unbound`.
These unbound receipts authenticate no complete input/kernel link and verify
no scientific or positive certificate. A completed proof still requires both
artifacts and their matching digests. Final output-cap failures now attempt a
small incomplete receipt under the same cap so a batch can continue; all other
I/O errors and non-overwrite checks remain failures. The revised suite passed624 assertions across15 test sets. The standalone
admission probe failed on the immutable old checker (1 pass,4 expected failures)
and passed on the revised checker (11 passes). The new small-cap and sparse-file
checks exercise actual guarded driver paths. All three attempts, including the
expected failure, are retained in [the Address validation](evidence/address_2026-10-28/VALIDATION.json).

The existing independent scientific replay retains stricter cumulative
accounting than the producer's separate component budgets. Its reconstructed
projection stays live while uncached rational Taylor replay shares the same
monotone Budget. A valid near-cap producer proof may therefore be refused by
the checker under the same Caps. Improving that replay efficiency requires
separate future work with explicit live-storage accounting; no Budget or cap
has been weakened or reset in this driver fix.

Historical evidence, producer source hashes, numerical results and execution
freezes are unchanged. Current-source regression fixtures are newly generated.
Because executable-source identity is enforced, the revised active checker
does not retroactively accept old producer records as if the new code had
produced or checked them.

No new molecular performance winner, universal approximation algorithm or
physical-error guarantee follows from these results.
