# Prospective reference correction and ordering comparison

Read PROTOCOL.md for the predeclared inputs, controls, resource limits and
decision rule. This module is numerical research code; it does not alter
QHAT's default algorithms, the manuscript, or the separate JW/BK project.

Use the existing Python 3.12 numerical environment. The run records exact
package versions. Dependencies are NumPy, SciPy, Numba, PySCF, OpenFermion,
OpenFermion-PySCF, h5py and pytest (tests only). Do not install QHAT's optional
Julia/resource-estimation stack solely to execute this comparison.

From this directory, run the following phases through guarded.py, one at a
time. Each label must be new. Supply absolute interpreter, QHAT and output
paths; use a new output directory after changing scientific source.

```
python guarded.py --label tests --logs ABS_LOGS -- ABS_PYTHON -m pytest test_engine.py -q
python guarded.py --label generate --logs ABS_LOGS -- ABS_PYTHON campaign.py generate --qhat-root ABS_QHAT --output ABS_OUTPUT
python guarded.py --label prepare --logs ABS_LOGS -- ABS_PYTHON campaign.py prepare --qhat-root ABS_QHAT --output ABS_OUTPUT
python guarded.py --label run --logs ABS_LOGS -- ABS_PYTHON campaign.py run --qhat-root ABS_QHAT --output ABS_OUTPUT
python guarded.py --label verify --logs ABS_LOGS -- ABS_PYTHON verify.py --qhat-root ABS_QHAT --output ABS_OUTPUT
python guarded.py --label summarize --logs ABS_LOGS -- ABS_PYTHON campaign.py summarize --qhat-root ABS_QHAT --output ABS_OUTPUT
```

Generation freezes source/protocol/environment before SCF, then authenticates
all three inputs. Preparation freezes every ordering before propagation.
The physical grid has 216 rows (three geometries, two determinants, nine
arms, four budgets); all circuits retain the entire Hilbert space. Verification
uses literal tensor-Pauli phases, CAR matrices, completed-unit Taylor products,
full reference integration, independent exponential actions and a separate
rotation implementation. Source review and numerical checking have separate
scope; neither is a new Lean proof or a certificate of optimal physical error.

one_pass.py is the unchanged exact supplied-weight placement rule from the
preceding diagonal-rule experiment. Its original SHA256 is
3cfaced88f188d418ad5fac1be585478ad7894d6d4da53142489f82c8e82fcb6.
The protocol manifest records its actual bytes and all active dependencies.
The corrected rule uses a different reference with complex response features;
it does not inherit the current rule's integer-score guarantee.

The sparse feature construction avoids full-H diagonalization and exponential
state allocation during search. Its pair/destination tables and complex
arithmetic can nevertheless be expensive: report measured setup and search
costs, not an assumed speed advantage. Full matrices used for independent
validation are deliberately limited to eight qubits.

After the completed campaign, start with [the decision report](report/README.md).
The compact frozen inputs, designs, states and verification ledgers are retained
in `evidence/results/`. The focused CI lane replays a scratch copy, preserving
the committed evidence and checking the independently derived decision summary.
