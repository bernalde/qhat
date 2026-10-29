# Prospective normal-window comparison

This study tests the opt-in normal Julia insertion backend on four larger 14-mode/four-electron instances in two previously seen families. The existing twelve-mode development results remain a separate completed study. This README and PROTOCOL.md were prepared before molecular generation; execution is established only by retained freezes, receipts and the eventual report.

The registry is fixed before generation; final source/configuration/resource freezes and a bounded end-to-end altered-input preflight are required before the main comparison. Physical errors are computed only after every constructor output is frozen. Full-space circuits, strong complete-unit and Pauli controls, total construction accounting, and retained losses/refusals are required. See PROTOCOL.md for the current preparation contract; preparation does not assert execution or completion.

Use the existing pinned scientific Python environment and Julia1.10.11, with
one numerical thread. No dependency installation or chemistry download is part
of these commands. Every numerical/build/test command is root-dispatched through
the frozen quarter-core/512MiB guard; the following describes the stage order:

1. Commit all scientific source. Run `generate_inputs.py --prepare --output DATA`,
   each of the four declared `--case` jobs, then `--finalize`. `DATA` must be a
   fresh child of the common campaign artifact directory so guard logs and
   finite resource ledgers can already exist outside it.
2. Run `preflight.py DATA NEW_PREFLIGHT_DATA`; this selects the largest retained
   word count, lexical tie, and encodes exact coefficient scaling by7/5.
3. Run `construct_campaign.py plan --root CAMPAIGN --inputs DATA_OR_PREFLIGHT
   --mode main|preflight --output NEW_PLAN --python PYTHON --julia JULIA`.
   This writes command arrays; it never launches them. Dispatch each fresh
   command once, preserving its named receipt and result or honest missing result.
4. Run `construct_campaign.py freeze --root CAMPAIGN --plan PLAN --output NEW_FREEZE`.
   All64 main calls (or8 preflight calls) need retained receipts before any
   per-case physical index is published. Missing orders remain censored.
5. Run `panel.py PANEL_INVENTORY NEW_PANEL_DIRECTORY`. It emits all88 main
   candidates per case, or44 preflight candidates, only for the complete valid
   roster. An incomplete inventory emits an audit only and exits2.
6. Run `objective_check.jl INPUT_TOML CANDIDATES_TOML NEW_CHECK`, then
   `physical.jl INPUT_TOML CANDIDATES_TOML NEW_PHYSICAL_DIRECTORY` and
   `independent.jl INPUT_TOML PHYSICAL_MANIFEST_TOML NEW_CHECK` under their
   own finite reservations. Main execution additionally requires the frozen
   resource-only preflight gate to pass; an easier case is never substituted.

The campaign planner preserves the selected Python virtual-environment path
without resolving its interpreter symlink. The computational kernels are Julia;
the data/identity glue and unchanged Python control methods remain explicit.
