# Sparse ordering on larger molecular cases

This directory provides the frozen twelve-mode comparison for the reusable
`analysis.sparse_response_ordering` library. The [completed evidence guide](evidence/README.md)
and [full six-arm report](evidence/report/REPORT.md) retain 24 orders, 96 actual
circuits, eight references and 240 independently verified comparisons. Read
[PROTOCOL.md](PROTOCOL.md) and [REGISTRY.json](REGISTRY.json) before reproducing it.

Cubic response has the lowest or tied state error in all 16 case/time/budget
settings. Its advantage over first-response search is strongest for CO; stretched
N2 ties that cheaper search, and CH4 improves only about 1% over cheap current
placement. Actual construction CPU is 3.4–30 seconds for cubic. These are
dependent observations from four active-space cases in three families, with no
universal winner, literature-novelty or hardware claim. All circuits use the
symmetric formula; first-response and individual-Pauli controls use a
first-response search objective.

The frozen design contains four previously unrecorded geometries from three
molecular families, six ordering arms, two validation times and two rotation caps. It
retains all 24 orders, 96 actual primitive circuits, eight reference evolutions
and 240 paired comparisons. Exact reference propagation is an offline checker;
training reads only the Hamiltonian, HF reference and declared training budget.

The sparse trainer allocates reached configurations rather than Hilbert-space
arrays. This benchmark's independent validator intentionally uses a bounded
4096-state calculation. The separately versioned line compiler preserves the
earlier routing/cleanup rules while extending its validated register domain to
twelve modes. Actual primitive circuits are checked through an exact Clifford
frame and independent Pauli rotations, with literal per-gate small-case checks.

Use the pinned Python 3.12 numerical environment and run every command through
the task's shared operating-system resource guard. Commit all scientific source
before chemistry generation or order freeze. Use fresh directories and retain
every failed or partial attempt:

```sh
python -m pytest numerical_audit/scalable_ordering/test_pipeline.py -q
python numerical_audit/scalable_ordering/generate_inputs.py --output FRESH_INPUTS
python numerical_audit/scalable_ordering/preflight.py --inputs FRESH_INPUTS --output FRESH_PREFLIGHT
python numerical_audit/scalable_ordering/study.py prepare --inputs FRESH_INPUTS --output FRESH_STUDY
python numerical_audit/scalable_ordering/study.py run --output FRESH_STUDY
python numerical_audit/scalable_ordering/study.py summarize --output FRESH_STUDY
python numerical_audit/scalable_ordering/verify.py --output FRESH_STUDY --destination FRESH_VERIFICATION
```

The method is a bounded heuristic, not a globally optimal order or a finite-error
certificate. The library's separate conservative step selector provides an
ideal-product-formula bound under its stated assumptions. Actual setup, search,
compilation, validation, refusals and output charges must remain distinct.
No chemistry-model/basis error or efficient correlated-state preparation claim
is established by these HF experiments.

The representative preflight selects the generated case with the largest
Pauli-word count (lexical case tie), scales every supplied coefficient by 7/5,
and tests all six arms at time 1.5 and cap 30U. It freezes those orders before
preflight errors. This resource-only selection and altered Hamiltonian do not
expose any main-case error; all preflight outputs and failures remain retained.
