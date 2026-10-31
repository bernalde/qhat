# Replaying the archived molecular evidence

These instructions use the positional interfaces from the completed archive replay. The path-variable examples below document those interfaces; they have not been run as an additional experiment.

The authoritative execution input is `archive/main-a1.data.tar.gz`, SHA-256 `4d1639d4c8c9ae29226d0b7e8bf1112897ad5c439327966baff8512f67da1892`. Loose `data/` files are selected convenience copies and are not a complete execution root. Extracting the archive yields `data/`, its byte inventory, original states, inputs, panels, witnesses and full `source_snapshot/` closure. `package.py extract` checks every regular member against the sealed inventory and writes `EXTRACTION_VERIFIED.json`. Use a fresh destination; do not rewrite original or extracted DATA.

The numerical producer is `aea354e287f2277d8b11cc51695a6d49f7e84537`, with freeze `f3ed635ac7831edcb3feb51d96e6116c2c95f7ad9e2b169b006ca862006776e3`. Use Python 3.12.12 and Julia 1.12.6, the frozen `analysis/julia_window` project, one Julia/GC/BLAS thread and the recorded dependency environment. The guard requires the documented Linux/systemd user-service environment. Retain actual runtime hashes and all receipts for a new verification run; a relocated executable path is not the original execution identity. This recipe installs nothing and does not run chemistry or repeat a constructor.

Set these variables to absolute local paths. `MG_BOOT` is the study directory in a trusted full checkout with byte-matched `package.py`, `common.py`, guard and referenced source closure. `MG_REPLAY` must not already exist. The bootstrap code is needed to authenticate/extract the archive; numerical commands subsequently use the **extracted** snapshot.

```bash
MG_PYTHON=/absolute/path/to/python3.12
MG_JULIA=/absolute/path/to/julia-1.12.6
MG_BOOT=/absolute/path/to/checkout/numerical_audit/molecular_generalization
MG_BUNDLE=/absolute/path/to/evidence/main_2026-10-30
MG_REPLAY=/absolute/path/to/fresh-verification-run
MG_ARCHIVE="$MG_BUNDLE/archive/main-a1.data.tar.gz"

sha256sum "$MG_ARCHIVE"
mkdir "$MG_REPLAY"
mkdir "$MG_REPLAY/logs"
mkdir "$MG_REPLAY/independent" "$MG_REPLAY/objective"
```

Compare the printed SHA-256 to the value above and the delivery manifest before extraction. The following guard function uses the executed dispatcher's actual guard CLI and numerical profile. Record prospective caps and a fresh immutable plan before dispatch; keep every failed cost. The original replay reserved 900 CPU/4510 wall seconds in a 1000/5000 package pool for seven stages. A new verification attempt has its own plans/receipts/costs and does not overwrite or reset the original campaign. Keep a finite cumulative allowance; these individual caps are not unlimited retries. The guard's numerical profile preserves quarter-core quota, 448/512 MiB RAM and no swap.

```bash
mg_guard() {
    local mg_label="${1//_/-}" mg_cpu="$2" mg_wall="$3"
    shift 3
    "$MG_PYTHON" "$MG_BOOT/tools/resource_guard.py" \
        --label "$mg_label" --logs "$MG_REPLAY/logs" \
        --cpu-seconds "$mg_cpu" --wall-seconds "$mg_wall" \
        --disk-mib 1024 --memory-profile numerical --output "$MG_REPLAY" \
        -- "$@"
}

mg_guard portable-extract 10 60 \
    "$MG_PYTHON" "$MG_BOOT/package.py" extract \
    "$MG_ARCHIVE" "$MG_REPLAY/extracted"

MG_DATA="$MG_REPLAY/extracted/data"
MG_SNAPSHOT="$MG_DATA/source_snapshot"
MG_STUDY="$MG_SNAPSHOT/numerical_audit/molecular_generalization"
```

Require a successful extraction and its authenticated inventory before continuing. The guard function converts underscores to hyphens only in service labels; scientific case identifiers and paths remain unchanged. Run the following for **each** fixed case, first `HF_r1p00_631g_12m8e`, then `H2S_r1p40_a96_631g_12m8e`, sequentially. Each output directory is fresh and outside extracted DATA. These are the exact positional interfaces used in the executed dispatcher's `commands` map: `verify.jl MODEL CANDIDATES PRIMARY_DIRECTORY OUTPUT` and `objective.jl MODEL CANDIDATES OUTPUT WITNESSES`.

```bash
MG_CASE=HF_r1p00_631g_12m8e

mg_guard "portable-verify-$MG_CASE" 240 1200 \
    "$MG_JULIA" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
    --heap-size-hint=256M --project="$MG_SNAPSHOT/analysis/julia_window" \
    "$MG_STUDY/verify.jl" \
    "$MG_DATA/inputs/$MG_CASE/model.toml" \
    "$MG_DATA/panels/$MG_CASE.toml" \
    "$MG_DATA/physical/$MG_CASE" \
    "$MG_REPLAY/independent/$MG_CASE"

mg_guard "portable-objective-$MG_CASE" 180 900 \
    "$MG_JULIA" --startup-file=no --history-file=no --threads=1 --gcthreads=1 \
    --heap-size-hint=256M --project="$MG_SNAPSHOT/analysis/julia_window" \
    "$MG_STUDY/objective.jl" \
    "$MG_DATA/inputs/$MG_CASE/model.toml" \
    "$MG_DATA/panels/$MG_CASE.toml" \
    "$MG_REPLAY/objective/$MG_CASE" \
    "$MG_DATA/witnesses/$MG_CASE.toml"
```

Then set `MG_CASE=H2S_r1p40_a96_631g_12m8e` and execute the same two commands once. Do not treat every nonzero exit as the known censor or rerun until an exit changes. Expected outcomes are:

| Fresh numerical command | Required payload outcome |
|---|---|
| HF `verify.jl` | Exit 0; 72 logical rows, 52 distinct unfused circuits and two references checked |
| H2S `verify.jl` | Exit 0; same denominators |
| HF `objective.jl` | Exit 2 / `partial_or_censored`; 72 selected scores and seven proposed values agree, 14 unique J checks; exactly original `rep2.window.B64.raw1` remains unpromoted/native-refused |
| H2S `objective.jl` | Exit 0 / `complete`; 72 selected scores and seven raw claims verified, 13 unique J checks |

Inspect each actual report and terminal resource receipt. HF's expected exit 2 must have zero mismatches/floor censorings/failures and the **exact** original native refusal `SearchDeadline` / “whole-call window allowance exhausted,” with no original native score. A guard termination, cleanup/resource failure or additional refusal is not an expected success. The witness file must not be omitted, and an expected HF payload exit must not be hidden by a generic shell “ignore errors” command.

These commands replay numerical claims; they do not alone certify all historical metadata, original receipt joins, CSV readback or publication acceptance. The executed full consumer sources are retained in `replay/source/`, together with `replay/REPLAY_FINAL.json`, `replay/consumer/REPLAY_CHECKS.json`, all plans/receipts and the source bindings. Consult the bundle manifest for the exact storage map.

`replay_dispatch.py` is a host-bound historical operator, not a portable CLI: it embeds the original ROOT/RUN/OUT paths and requires the original producer/visual gate. Do not execute it blindly on another machine or rewrite its archived bytes. `metadata_consumer.py` has the real argument set `--producer-data`, `--data`, `--consumer-root`, `--numerical-root`, `--producer-analysis`, `--producer-visual-gate`, `--metadata-plan`; it also validates exact original owner/gate/source relationships. Reusing it elsewhere requires a separately labelled, prospectively bound later-consumer plan and verified path mapping. These instructions do not invent a generic portable “all gates pass” runner.

Original and replayed scientific outcomes remain `ANALYSIS.complete=false`, `scientific_gate_complete=false` and both primary decisions `incomplete`. The recorded `delivery_complete=true` / `replay_complete=true` mean that all required processing reproduced the retained limitation; they do not supply the missing original native check. No constructor retry, native backfill, promotion or scientific-rule change is part of reproduction.
