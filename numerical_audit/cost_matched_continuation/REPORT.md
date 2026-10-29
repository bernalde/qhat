# Cubic continuation under equal total CPU allowances

Full-order cubic continuation produced lower physical error than Julia width-eight window refinement in all 32 paired observations on the two declared new 14-mode/four-active-electron geometries. The predeclared primary rule selects `full_advantage`: both case medians favor full continuation by more than the fixed `1e-12` ambiguity, every primary repetition favors full continuation, and all source, execution, cost, objective and independent physical gates pass. This is a result for the declared methods and allowances on these two geometries; it does not establish a universal ordering method.

The [complete results](evidence/results/RESULTS.json) retain all 20 frozen orders, 80 logical physical rows, 32 paired observations, work/stopping records and accounting. Their authoritative main analysis is `ANALYSIS.json`, SHA-256 `0855a33c981b379d40d92853be38761b81e056dbca2dece9d62a8214eb3b899f`, with `status=complete` and no issues.

The [protocol](PROTOCOL.json) and [registry](REGISTRY.json) fixed both geometries, training objective, allowances and physical conditions before observing their physical errors. LiH has bond length 1.70 Angstrom, no frozen core and active spatial orbitals 0–6. Linear symmetric BeH2 has hydrogen coordinates ±1.60 Angstrom, frozen spatial orbital 0 and active orbitals 1–7. Both use RHF/6-31g with 14 modes and four active electrons. These are fresh geometries in previously studied families, with no new-family claim. The geometry audit is limited to recorded numerical-audit JSON at baseline `646c029d71fa7abd2c06ae64693843f4aeb5f1e8`. Energies are in Hartrees and times in inverse Hartrees, with hbar=1.

Each case has two deterministic timing repetitions. Within a repetition, both methods continue exactly the same measured frozen parent at total CPU allowances B=32 and B=64 seconds. Full continuation uses four deterministic starts prepared inside the measured parent process; its insertion search may change every position, including global endpoints. The Julia arm revalidates the supplied parent and searches strictly interior windows of width eight, with at most 32 windows and 100 sweeps. Native evaluation/move ceilings are one million/1,000; all declared storage and support limits remain fixed. The two search domains and amounts of arithmetic differ. Restart seeds depend on case and repetition, so the repetitions include fixed restart variation as well as timing variation.

Both methods train and are evaluated from the same Hartree–Fock determinant, with occupied active modes 0–3. Training uses the cubic response objective at T=1 and a 30U ideal original-word rotation cap. Orders and panels were frozen before physical evaluation at T=0.65 and 1.35, with 8U and 24U caps and the symmetric formula. U is the original Pauli-word count of the corresponding input. Circuits propagate the full 16,384-amplitude register; only the numerical exact-H reference uses a checked number sector. The affordable Strang step counts are four at 8U and twelve at 24U for every retained row. Physical error is the stable evaluation of one minus the absolute overlap of the normalized exact and Trotter states. Norm drift and number/spin leakage are retained separately.

The primary condition is B=64, T=1.35, 8U. Lower projective physical error is better; Δ below is window error minus full error. These displayed values reproduce the retained analysis.

| Case | Repetition | Full error | Window error | Δ |
| --- | ---: | ---: | ---: | ---: |
| LiH, 1.70 Angstrom | 1 | 9.761260607563768e-8 | 1.1202312396998652e-7 | 1.4410517894348847e-8 |
| LiH, 1.70 Angstrom | 2 | 1.0248259266941073e-7 | 1.1202312396998652e-7 | 9.540531300575794e-9 |
| BeH2, ±1.60 Angstrom | 1 | 1.2250890185845893e-7 | 1.9628895487952813e-7 | 7.37800530210692e-8 |
| BeH2, ±1.60 Angstrom | 2 | 1.4211340811467015e-7 | 1.9628895487952813e-7 | 5.417554676485798e-8 |

The median of the two already-paired differences is `1.197552459746232e-8` for LiH and `6.397779989296359e-8` for BeH2. Both are positive and exceed `1e-12`; none of the four primary pairs opposes that direction. These are paired medians, not differences of method medians.

Across both times, rotation caps and repetitions, each CPU allowance contributes 16 pairs:

| Total CPU allowance | Full wins | Ties | Window wins | Unavailable | Median paired Δ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32 seconds | 16 | 0 | 0 | 0 | 6.452569566568908e-10 |
| 64 seconds | 16 | 0 | 0 | 0 | 6.690282429172729e-10 |

These aggregate medians mix the declared physical conditions and are descriptive. All 32 pairs are eligible; none is omitted. No physical error is at or below the `1e-12` ratio floor: the smallest of all 80 retained errors is `2.1127083732513202e-11`. Ratios are withheld at the declared floor when applicable; the rule never discards the paired errors or losses. The two timing repetitions are not independent molecular samples, and no confidence intervals are reported.

Equal allowances did not produce equal actual spending. All eight full continuations stopped at `cpu_limit` after scoring all four declared starts. Six window continuations stopped at the child CPU limit and retained verified checkpoints. The two BeH2 B=64 window runs stopped at the declared 200,000-atom limit, at total conservative charges of 50.463 and 48.698 CPU seconds, retaining nine checked updates each. Those support stops are valid interrupted checkpoints, not convergence. The four parent searches reached their intentional 3,000-evaluation ceiling and retained their measured restart plans. No support cap, geometry, allowance or search setting was retuned and no fresh-budget constructor retry was made.

Complete construction costs include native Python worker startup, imports/loading, restart preparation, search, scoring, serialization and teardown, together with owned Julia 1.10.11 startup/JIT and child work. The pinned Python scientific environment and exact source bytes are archived. Runs used one numerical thread, a quarter-core quota, 448/512 MiB MemoryHigh/MemoryMax and no swap. Each continuation's total charge includes its frozen parent prefix and passes the B CPU and 5B wall ceilings.

| Construction accounting boundary | Raw CPU seconds | Conservative CPU seconds | Guarded wall seconds |
| --- | ---: | ---: | ---: |
| Actual study executions: four shared parents once + sixteen continuations | 688.008 | 688.028 | 2782.919 |
| Attributed comparisons: parent prefix charged to each arm/budget continuation | 732.144 | 732.176 | 2982.081 |

Wall displays are rounded to three decimals; full precision remains in the results. For each parent P and its four continuations C, actual construction cost is P+ΣC and attributed comparison cost is 4P+ΣC. Physical rows never multiply constructor charges. The four unique parents cost 14.716 conservative CPU seconds in total; repeated attribution adds 44.148 CPU seconds. These are constructor costs, not total project CPU or quantum execution costs. The separate 31-job main scientific cutoff reconciles inputs, constructors, panel assembly, physical evaluation, objective replay and independent checking at 816.058 conservative CPU seconds and 3310.320 guarded wall seconds. Prior preflight/validation attempts and later analysis, reporting, packaging and review costs remain separate.

Independent checking covers all 80 logical physical rows through 44 distinct unfused circuit replays: 20 for LiH and 24 for BeH2. Four independent reference checks cover the two times for each case. Objective replay verifies all 80 scored logical rows with zero unscored rows or mismatches; exact deduplication yields five unique objectives for LiH and six for BeH2. Deduplication reduces numerical tasks without removing logical observations. Original state bytes, source/input hashes, objective witnesses, resource receipts and guard joins are retained with the results.

The practical conclusion is that, from these common parents and under the declared 32/64-second CPU allowances and fixed capacity limits, full-order continuation gave better physical errors on both new geometries. In particular, the BeH2 B=64 result compares full continuation against a window method that stopped early at its declared support limit; it is not an experiment forcing both methods to spend exactly 64 seconds. No local/global optimality, molecular approximation guarantee, hardware-performance conclusion or automatic default-optimizer replacement follows. This compares cubic insertion searches; the separate exact first-order dynamic-programming results are a different study.

The [evidence guide](evidence/README.md), [copy manifest](evidence/COPY_MANIFEST.a1.json) and [main archive](evidence/archives/main-a1.tar.gz) document the retained payload and its source bindings. All five archives passed independent extraction with 2,174 member-byte checks. Replaying the frozen analyzer against the relocated main DATA_ROOT from outside the checkout reproduced the entire original analysis byte-for-byte, including its SHA-256 above. This consumer check validates retained evidence joins; independent numerical propagation was performed during the original study. [Delivery evidence](evidence/results/DELIVERY_EVIDENCE.json) records the check sources and receipts. The frozen prospective README and protocol retain their original chronology.
