# Exact source used for this run

This file matches the source hash recorded in the run manifest before execution.
Subsequent review added three provenance checks/fields to the live runner:
verification of every fixture file, an explicit first-order/chronology label,
and the discarded scalar coefficient. Numerical evolution was unchanged.
The updated runner passed a separate three-ordering HCN smoke run.

For byte-identical source reconstruction, overlay this archived file at
`numerical_audit/run_explicit_precision.py` in a separate checkout. It is an
archive, not a standalone script to execute from this source_snapshot directory.
Other source hashes match the current audit implementation.
