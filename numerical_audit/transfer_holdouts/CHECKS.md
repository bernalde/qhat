# Transfer implementation checks

All 78 source checks passed at implementation commit
`7b99593eb95e86ef3ea40a6c8621d6f34c08823f`. The retained worker receipt is
`transfer-full-tests-a1.resources.json` and its JUnit artifact is
`runtime/transfer-full-tests-a1/junit.xml` in the owned campaign area.

The suite includes the original response/state/training checks, literal coherent
objectives and signed cubic columns, exact integer enumeration of every legacy
initializer slot, both formula orders and circuit orientations, real/imaginary
states, full-space number leakage, independent primitive and unfused replay,
DAG resources, maximal affordability, compilation censors, frozen-source/identity
and allocation guards, all 1,440 logical identities and exact work sharing,
reviewed-interface gating and complete outer caller cost boundaries.

Run the six files with the existing pinned interpreter:

```sh
python -m pytest -q numerical_audit/transfer_holdouts/test_cubic_response.py \
  numerical_audit/transfer_holdouts/test_states.py \
  numerical_audit/transfer_holdouts/test_compilation.py \
  numerical_audit/transfer_holdouts/test_panel.py \
  numerical_audit/transfer_holdouts/test_study.py \
  numerical_audit/transfer_holdouts/test_objective_reference.py
```

The worker used 2.21 measured CPU seconds and 9.064973325 wall seconds, with
MemoryMax=512 MiB, no swap, quarter-core CPU and one thread. Cumulative unit-check
usage is 5.54 CPU seconds and 30.223323790 wall seconds against the unchanged
original 120/900 allowance. Earlier attempts remain in the campaign ledger.

These checks validate the implementation on bounded toy fixtures and the exact
panel definitions. They do not establish transfer on the six real molecular
inputs. The independently checked 30 reference states are separate evidence.
The full representative preflight, real training/validation panel, independent
complete replay, archive/Git delivery and paper checks remain required after
issue #16's reviewed interface is available. Paired source reviews, publication,
separate eligible merges and their actual CI are also pending.

Current development adds the whole-panel independent verifier and an explicit
resource-only partial-trace audit after the tested commit above. Their added
integration checks are queued behind the active issue #16 replay. The retained
78-check result remains valid historical evidence; it does not claim acceptance
of those later executable changes. No real transfer validation has run.

Reporting/resource-ledger and deterministic delivery code were added later as
well. Their corruption/refusal checks are queued with the updated suite; complete
end-to-end delivery still needs the declared representative preflight.
