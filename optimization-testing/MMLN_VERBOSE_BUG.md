# MMLN `verbose=TRUE` Bug: Variable Shadowing Corrupts Random Intercepts

## Summary

The MMLN function's progress bar code overwrites the variable `m` (number of
random-effect groups) with the ETA minutes value. This causes the random
intercept update loop to run the wrong number of iterations for the remainder
of the MCMC chain, producing **incorrect posterior samples** whenever
`verbose=TRUE`.

This bug exists in both the upstream `eaegerber/MMLN` and our fork.

## Location

`R/mln_functions.R`, inside the MMLN function's verbose progress block
(approximately line 378):

```r
# m is initialized as ncol(Z) — the number of groups (e.g. 50)
m <- ncol(Z)

for(it in seq_len(n_iter)) {
    # ... psi update uses m ...
    for(j in seq_len(m)) {     # <-- loops over random-effect groups
      # ...
    }

    # ... later in the same loop body ...
    if (verbose && ...) {
      # ...
      h <- floor(eta_sec / 3600)
      m <- floor((eta_sec %% 3600) / 60)   # <-- OVERWRITES m with ETA minutes!
      s <- round(eta_sec %% 60)
      # ...
    }
}
```

After the first verbose update fires, `m` changes from the number of groups to
the minutes component of the ETA. For a run taking 3 minutes, `m` becomes `3`
instead of `50`. For a fast run (<60s), `m` becomes `0`, and the psi loop is
skipped entirely for all subsequent iterations.

## Impact

- **MMLN results with `verbose=TRUE` are non-reproducible** across runs because
  the corrupted `m` value depends on wall-clock time
- The random intercepts (`psi`) stop updating after the first progress bar tick,
  so the posterior samples for `psi`, `Phi`, and downstream chains are incorrect
- The Phi update (`df = nu_P + m`) also uses the corrupted `m`, compounding the
  error
- **FMLN is not affected** because it does not use `m` in its loop body (no
  random intercepts)

## How it manifests in test results

Our test suite runs both reference (upstream) and fork with `verbose=TRUE`. Both
have the bug, but `m` gets a different corrupted value each time because
wall-clock speed differs:

- **norm proposal** (~2-8s): ETA is always <60s, so `m = 0` in both runs →
  same wrong value → chains match (bitwise identical)
- **beta proposal** (~130-265s): ETA varies between runs, so `m` gets different
  minute values → chains diverge
- **normbeta proposal**: same as beta, plus the separate C++ vs R divergence

This is why MMLN beta fails on ALL scenarios despite the fork not changing any
beta code — it's not a fork bug, it's a verbose bug.

## Fix

Rename the ETA variable to avoid shadowing. In both FMLN and MMLN:

```r
# Before (buggy)
h <- floor(eta_sec / 3600)
m <- floor((eta_sec %% 3600) / 60)
s <- round(eta_sec %% 60)
eta_str <- sprintf("%02d:%02d:%02d", h, m, s)

# After (fixed)
h <- floor(eta_sec / 3600)
mn <- floor((eta_sec %% 3600) / 60)
s <- round(eta_sec %% 60)
eta_str <- sprintf("%02d:%02d:%02d", h, mn, s)
```

This should be reported upstream to `eaegerber/MMLN`.

## Verification

```r
# This should print 0 after the fix:
set.seed(42)
res1 <- MMLN(..., verbose = FALSE)
set.seed(42)
res2 <- MMLN(..., verbose = TRUE)
max(abs(res1$w_chain[[N]] - res2$w_chain[[N]]))
```
