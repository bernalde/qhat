# Read-only Lean bridge appendix correspondence

Passed: all 19 public theorem names appear exactly once, in source order, in the new manuscript appendix. The mathematical steps follow the frozen proof guide and actual Lean source. Source hashes and the complete statement list are recorded in lean_bridge_correspondence_audit.json. No Lean recompilation was performed or claimed in this transcription check.

The support-to-weighted-shift derivation, literal matrix multiplication, pair commutator, molecular-mask specialization, actual orbital-sign conjugation, diagonal energy-gap identity and list BCH accumulation all preserve their Lean hypotheses and orientation. In particular, the scalar diagonal identity requires no support condition, whereas its commutator-entry specialization does; the appendix correctly distinguishes these statements. Repeated labels and repeated masks are allowed, and the label type need not be finite.

The proof of the BCH-entry identity expands exactly the imported head-tail recursion `comm head tail.sum + bch tail`. The prose invokes the earlier entry specialization for readability where Lean simplifies the same shift conditional directly; this is an equivalent mathematical step, not an alternative hypothesis or result.

The appendix correctly labels these as supporting formalizations rather than 19 new contributions and explicitly leaves the sinc integral, variance identity, complete response objective, ranking condition and numerical search outside the Lean bridge. No mathematical or source-correspondence defect was found. The external Ryabinkin locator and the prior 42-constant audit claim were not independently re-reviewed in this source-transcription check.
