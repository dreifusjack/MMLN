# MMLN `verbose=TRUE` Bug: Variable Shadowing Corrupts Random Intercepts

## Summary

The MMLN function's progress bar ETA code reuses the variable name `m`, which
already holds the number of random-effect groups. After the first progress
update, every subsequent MCMC iteration uses the wrong value of `m`, producing
**incorrect posterior samples** whenever `verbose=TRUE`.

This bug exists in both the upstream `eaegerber/MMLN` and our fork.

## What `m` is and where it's used

`m` is defined at **line 235** of `R/mln_functions.R`:

```r
N <- nrow(Y); d <- ncol(Y) - 1; p <- ncol(X); m <- ncol(Z)
```

`m` is the number of columns in the random-effects design matrix `Z` — i.e.,
the number of groups (e.g., 50 teams in the baseline scenario). It is used in
two critical places inside the MCMC loop:

1. **Line 309** — the random intercept (`psi`) update loops over all `m` groups:
   ```r
   for(j in seq_len(m)) {
       idx <- which(Z[,j] == 1)
       # ... update psi[j,] using observations in group j ...
   }
   ```

2. **Line 330** — the Wishart degrees of freedom for the `Phi` (group covariance) update:
   ```r
   df = prior_settings$nu_P + m
   ```

## The overwrite

At **line 378**, inside the verbose progress block, the ETA formatting code
assigns the number of **minutes** remaining to a variable also called `m`:

```r
if (verbose && (it %% max(1, floor(n_iter / 100)) == 0 || it == n_iter)) {
    setTxtProgressBar(pb, it)
    elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
    eta_sec     <- elapsed_sec / it * (n_iter - it)

    h <- floor(eta_sec / 3600)
    m <- floor((eta_sec %% 3600) / 60)   # <-- OVERWRITES m (group count) with ETA minutes
    s <- round(eta_sec %% 60)

    eta_str <- sprintf("%02d:%02d:%02d", h, m, s)
    cat(sprintf("\r ETA: %s", eta_str))
    flush.console()
}
```

## When the overwrite happens

With `n_iter = 1000`, the verbose block fires every `max(1, floor(1000/100)) = 10`
iterations. So at **iteration 10**, `m` is overwritten for the first time.

The psi update and Phi update run *before* the verbose block in each iteration,
so iteration 10 itself is still correct. But from **iteration 11 onward**, `m`
holds the ETA minutes value instead of the group count.

## What goes wrong after the overwrite

| Run duration | `m` becomes | Effect on psi loop (`seq_len(m)`) | Effect on Wishart df |
|-------------|-------------|----------------------------------|---------------------|
| < 1 minute  | `0`         | Loop runs 0 times — **psi never updates again** | df too small |
| 1-2 minutes | `1`         | Only 1 group updated per iteration | df off by ~49 |
| 3+ minutes  | `3+`        | Only 3+ groups updated | df still wrong |

For a typical run taking ~2 minutes, `m` becomes `1` or `2` instead of `50`.
The random intercepts for groups 2-50 freeze at their iteration-10 values for
the remaining 990 iterations.

## Why results differ between runs

The corrupted value of `m` depends on **wall-clock time** (`eta_sec`), which
varies with CPU load, system speed, and background processes. Two runs with the
same seed produce different `m` values → different posterior samples →
non-reproducible results.

This is why our test suite shows:
- **MMLN norm passes**: Runs take ~2-8 seconds. `m` becomes `0` on both
  reference and fork runs → same wrong value → chains match by accident.
- **MMLN beta/normbeta fail**: Runs take 2-4 minutes. `m` varies between
  reference and fork runs → different wrong values → chains diverge.

## FMLN is not affected

FMLN has the same ETA code at line 160, but FMLN has no random effects — it
never defines or uses `m`. The overwrite is harmless.

## The fix

Rename the ETA minutes variable from `m` to `mn` on lines 378 and 381:

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

After the fix, verbose and quiet runs should produce identical results:

```r
library(MMLN)
set.seed(42)
res1 <- MMLN(Y, X, Z, n_iter=100, burn_in=10, thin=2, mh_scale=1, verbose=FALSE, proposal="beta")
set.seed(42)
res2 <- MMLN(Y, X, Z, n_iter=100, burn_in=10, thin=2, mh_scale=1, verbose=TRUE, proposal="beta")
max(abs(res1$w_chain[[45]] - res2$w_chain[[45]]))
# Should be 0
```
