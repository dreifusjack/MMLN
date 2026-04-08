## Summary

End-to-end integration test suite for verifying fork optimizations against the
upstream `eaegerber/MMLN` package. Compares FMLN, MMLN, and MDres output
chain-by-chain across 10 synthetic scenarios and 3 proposal types.

**Scope change from ticket:** Instead of the 4 original datasets (small sim,
medium sim, MLB, pollen), we built 10 targeted synthetic scenarios that
stress-test specific numerical edge cases. Tests FMLN, MMLN, and MDres
independently with `FUNC=` and `SCENARIO=` filtering.

### What's included

- `optimization-testing/` directory with Makefile-driven workflow
- `make generate-data` / `make reference` / `make test` / `make compare` / `make clean`
- Resume support for reference generation (skips existing fixtures)
- Divergence trajectory plots for normbeta chains
- README updated with `Rcpp::compileAttributes()` build instructions

### Test results

| Function | Proposal  | Result |
|----------|-----------|--------|
| FMLN     | norm      | 40/40 PASS (bitwise identical) |
| FMLN     | beta      | 40/40 PASS (bitwise identical) |
| FMLN     | normbeta  | 36/40 PASS (1 scenario fails) |
| MMLN     | norm      | 60/60 PASS (bitwise identical) |
| MMLN     | beta      | all FAIL |
| MMLN     | normbeta  | all FAIL |

## Known issues found

### 1. Upstream bug: `m` variable shadowing in MMLN (`MMLN_VERBOSE_BUG.md`)

The MMLN progress bar code (`verbose=TRUE`) overwrites `m` (number of
random-effect groups) with the ETA minutes value. This corrupts the random
intercept loop for all subsequent iterations, producing incorrect posterior
samples. The bug exists in both upstream and fork, but produces different
results each run because the corrupted value depends on wall-clock time.

This is the root cause of **all MMLN beta and normbeta failures** — the fork
code is correct, but reference and fork get different corrupted `m` values.

**Fix:** Rename `m` → `mn` in the ETA formatting block. Should be reported
upstream.

### 2. Normbeta alpha overflow on `large_coefs` (`NORMBETA_DIVERGENCE.md`)

The alpha formula `(1 + exp(mu)) / sigma - exp(mu) / (1 + exp(mu))` overflows
to NaN when mu is large, desynchronizing the RNG stream. Only affects the
`large_coefs` scenario (beta in [-4,4]).

**Fix:** Rewrite using stable sigmoid form `1 / (1 + exp(-mu))`. Also remove
the `alpha < 1e-3` safety floor in C++ that doesn't exist in upstream R.
