# Recovered optimized and Strang benchmark campaign

This directory contains the durable recovery of the optimized-order comparison campaign after the WSL reboot erased the previous temporary staging area. No temporary numerical results were treated as surviving evidence. Existing reviewed input designs and exact target states are read from QHAT's committed artifacts; neither Hamiltonians nor reference states are regenerated here.

Read `PROTOCOL.md` for the frozen scientific comparison. There are 14 cases, four optimized first-order methods at seven rotation budgets and four Strang controls at six budgets: **728 propagated states**. Strang is omitted at budget U because at least one complete step generally costs more than U; its first tested budget is 2U. All reported work counts are actual Pauli rotations after exact commuting endpoint-block fusion. They are not compiled gate counts.

The supported commands, using the environment pinned by the enclosing benchmark extension, are:

```sh
python checked_prepare.py
python check_oracles.py
python -O check_resume_guard.py
python check_floating_guard.py
python execute.py freeze
python execute.py run
python execute.py finish
```

Use `--cases CASE ...` to select cases for `checked_prepare.py` and `execute.py run`. Resume with these same commands; do not bypass the provenance wrapper by calling low-level preparation functions. `QHAT_ROOT` can point to a relocated QHAT checkout. The surrounding `common.py` must remain next to this directory. Read `CHECKPOINT_RECOVERY.md` for the explicit adoption record of the twelve small cases prepared before the recovery guard was introduced.

The floating-point search correction is documented by `floating_guard_diagnostic.json`: the old incremental-score consistency guard aborts on F2 even though the proposed insertion agrees with a positive exact-amplitude rescore. The recovered optimizer uses exact integer destination amplitudes, an engineering roundoff allowance for screening, and nonnegative quadratic rescores before accepting candidate improvements. This is not a formally verified floating-point error bound. `oracle_checks.json` checks every insertion position for two starts on a 64-term molecular subset, exact equality of the specialized Pauli pair cache with the original class algebra, and fused versus unfused Strang products. All eight development JD, J0 and unrestricted Pauli JD schedules are also required to reproduce the archived schedules exactly.

`designs/` stores plans and optimizer diagnostics; `checkpoints/` stores per-start accepted-sweep checkpoints and input/source bindings. `ALL_SCHEDULES_FROZEN.json` and `EXECUTION_MANIFEST.json` are written before propagation. `rows/` stores one measurement record per budget, while `states/` stores the actual state and, for Strang, integer pulse indices and duration numerators. The shared denominator is 2r and each word's summed duration numerator is 2r. Rows and states are saved atomically, then hashed. `results.json` and `completion.json` are written only after all 728 rows exist.

Search budgets count all insertion destinations examined, including commuting swaps. The matched-budget Pauli control receives the conserving JD search's total destination count, split between its two starts. This is an evaluation-count control, not a wall-clock-time claim; setup and search times are separate. Independent circuit/objective verification and cross-method result analysis are maintained by the enclosing benchmark campaign.

Read `timing_scope.json` before comparing wall times. In this campaign, row `seconds` includes Strang pulse construction, state compression/output, and one design hash as well as evolution; fixed-campaign rows time evolution only. Pair setup is measured once and shared by JD/J0 or unrestricted/budgeted Pauli searches, even though its time is repeated in both arm records. Search counts are comparable algorithmic measures; these heterogeneous wall-clock intervals are descriptive and must not be presented as identical timing protocols.

## Completed checks and analysis

`archived_schedule_checks.json` verifies all 72 relevant literal archived schedules, and `archived_state_checks.json` verifies all 30 repeated old r=100 state vectors bitwise. `REPORT.md` and `summary.json` summarize the eight optimized/symmetric controls; `RUNTIME_REPORT.md` gives search, shared setup, propagation interval and memory measurements. `SOURCE_PROOF_CORRESPONDENCE.md` and `EMPIRICAL_CLAIM_REVIEW.md` record independent read-only reviews of the current written mathematics and empirical claims.

`boundary_sensitivity.py` is a separate cost-only analysis. It computes a uniform maximal pairwise-commuting prefix and suffix for each symmetric schedule, proves and checks overlap-safe exact fusion, and reaccounts the work of **existing** states at their **existing** step counts. `BOUNDARY_SENSITIVITY.md` and `boundary_sensitivity.json` report least-tested-work comparisons and explicit coverage of the separately running published V2 campaign. No original row or frozen schedule changes, and no unmeasured error at a newly feasible larger step count is inferred. Run `python boundary_sensitivity.py` after all V2 rows exist to refresh its coverage. The dense checks include overlapping endpoints and fully commuting schedules; complete-class-label and refined-label Strang formulas have identical flattened lists in all fourteen cases and therefore identical unitaries at the same step count.
