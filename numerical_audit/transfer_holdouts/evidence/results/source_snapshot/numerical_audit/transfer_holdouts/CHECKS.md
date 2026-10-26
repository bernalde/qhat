# Transfer implementation checks

All **103 source checks passed** at implementation commit
`3b8485d39e26d2b5eeb7bb7bc3f6d711f3704d52`. The retained receipt is
`transfer-full-tests-a3.resources.json`, with JUnit under
`runtime/transfer-full-tests-a3/junit.xml` in the owned campaign area.
This subsequent documentation update changes no tested Python source.

The nine-file suite covers coherent response and signed cubic columns, literal
integer legacy placement, both formula orders/orientations, full-space leakage,
primitive/unfused replay and DAG counts, maximal affordability and censors,
frozen identities/allocation guards, all 1,440 logical identities and exact work
sharing, reviewed-interface gating, complete caller-cost boundaries, whole-panel
independent verification, raw resource ledgers and deterministic archive,
extraction, inventory and Git-delivery refusal checks.

```sh
python -m pytest -q numerical_audit/transfer_holdouts/test_cubic_response.py \
  numerical_audit/transfer_holdouts/test_states.py \
  numerical_audit/transfer_holdouts/test_compilation.py \
  numerical_audit/transfer_holdouts/test_panel.py \
  numerical_audit/transfer_holdouts/test_study.py \
  numerical_audit/transfer_holdouts/test_objective_reference.py \
  numerical_audit/transfer_holdouts/test_verify.py \
  numerical_audit/transfer_holdouts/test_reporting.py \
  numerical_audit/transfer_holdouts/test_delivery.py
```

The successful worker used 2.21 CPU seconds and 9.048795731 wall seconds,
MemoryMax=512 MiB, no swap, quarter-core CPU and one thread. Cumulative unit-check
usage is 9.70 CPU seconds and 48.332141866 wall seconds against the unchanged
original 120/900 allowance. Every earlier attempt remains charged and retained.
In particular, attempt a2 had 61 passes and 42 fixture setup errors because the
coordinator omitted the custom temporary directory's parent. Creating that
parent in a fresh a3 output fixed setup; the source was unchanged. The earlier
78-check pass remains historical evidence at `7b99593`.

These are bounded implementation and corruption/refusal checks, not evidence of
transfer on the six real molecular inputs. All 30 reference states were checked
independently as separate evidence. No real transfer training or validation has
run. The representative full-path preflight, actual complete panel, independent
replay, archive/Git delivery, paper checks, paired reviews, publication, separate
eligible merges and actual-merge CI remain required. Real execution waits for
issue #16's legitimately reviewed interface; synthetic unit fixtures do not
satisfy that gate.
