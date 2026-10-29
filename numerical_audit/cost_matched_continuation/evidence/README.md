# Evidence for the matched-allowance continuation study

[RESULTS.json](results/RESULTS.json) is a byte-identical copy of the complete main analysis, SHA-256 `0855a33c981b379d40d92853be38761b81e056dbca2dece9d62a8214eb3b899f`. It contains all twenty order records, eighty logical physical rows, thirty-two pairs, independent reports and the scientific-cutoff ledger. [SCIENCE_ADMISSION.json](results/SCIENCE_ADMISSION.json) records the coordinator's final scientific inspection. The [study report](../REPORT.md) explains the outcome and its limits.

The five archives under `archives/` retain complete main and preflight DATA_ROOTs, their source snapshots and binary states, operator/validation snapshots, and every completed OS job at the packaging snapshot. [COPY_MANIFEST.a1.json](COPY_MANIFEST.a1.json) gives each archive's SHA-256 and every member's path, size and SHA-256. Source paths in that manifest are recorded provenance; the extractor never opens them. Historical configuration and source files remain inside archive wrappers.

The main scientific ledger is the immutable 31-job `main-a1/data/campaign/LEDGER_CUTOFF.json`. The outer main archive also contains the later analysis job. The packaging job could not contain its own future completion receipt; subsequent archive/extraction/replay checks are recorded separately in the delivery evidence. Original failures and earlier source/accounting versions remain retained. No task/handoff records, reviewer packets/reports, credentials or raw model streams are part of this payload.

From the repository root, use a fresh output directory whose parent already exists:

```sh
python numerical_audit/cost_matched_continuation/evidence/tools/extract_evidence.py \
  --evidence numerical_audit/cost_matched_continuation/evidence \
  --output /tmp/qhat-cmc-evidence
python /tmp/qhat-cmc-evidence/main-a1/data/source_snapshot/numerical_audit/cost_matched_continuation/analyze.py \
  --root /tmp/qhat-cmc-evidence/main-a1/data \
  --output /tmp/qhat-cmc-frozen-analysis.json
```

The extracted frozen analyzer is SHA-256 `ca1b9f5007f98581d52e4f4f20afd19d43902b4770d2e59dd6559ca10cd69686`; it is the consumer used for the original complete-output identity claim. To inspect the retained data with the maintained consumer, write a separate output:

```sh
python numerical_audit/cost_matched_continuation/analyze.py \
  --root /tmp/qhat-cmc-evidence/main-a1/data \
  --output /tmp/qhat-cmc-maintained-analysis.json
```

The maintained consumer inventories optional unpinned legacy logs but uses only authenticated evidence for terminal CPU and memory accounting. Its output records the executing analyzer version, including the recursively replayed preflight version. Keep this later output separate from the original `RESULTS.json`, source freezes and archive manifest.

[The post-review replay record](results/CONSUMER_REPLAY.json) authenticates all 2,174 extracted members before and after checking both consumers. The frozen consumer again reproduced the entire original result byte-for-byte. The maintained main result differs only at `analyzer_sha256` and `preflight_prerequisite.recomputed_analyzer_sha256`; a separate A3 replay differs only at `analyzer_sha256`. All other JSON values, including scientific outcomes, costs and admission gates, match. The record retains exact output hashes and separate regression/replay receipts; these later checks do not add numerical propagations to the study.

Extraction validates archive hashes, complete regular-member inventories and all member bytes before writing its success record. Outputs are fresh and partial failure evidence is preserved. The extractor requires POSIX filesystem operations. The consumer uses Python's standard library (`tomllib` requires Python 3.11 or newer); the demonstrated environment is the recorded Python 3.12 research interpreter. Original absolute command paths remain provenance in the data.

This is archive extraction and analyzer-consumption replay. Fresh numerical regeneration additionally requires the recorded scientific Python environment, Julia 1.10.11, archived producer source and the declared resource guard. Consumer portability does not imply a new physical propagation or a cross-platform bitwise numerical guarantee. The prospective README and protocol remain unchanged; completed results are in REPORT.md.
