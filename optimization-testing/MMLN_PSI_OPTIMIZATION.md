# MMLN Random Intercept Update: Optimization Notes

## What was optimized

The per-group ψ (psi) update is the dominant bottleneck in `MMLN()`, particularly for
large-m scenarios (e.g. `many_groups`: m=500, n_i=4, N=2000). The upstream loop ran each
of the following m times per MCMC iteration:

```r
for(j in seq_len(m)) {
  idx <- which(Z[, j] == 1)                                      # O(N) scan
  R_j <- R_tot[idx, , drop = FALSE]
  V_j <- chol2inv(chol(chol2inv(chol(Phi)) + length(idx) * Sigma_inv))  # 2 Choleskys
  M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
  psi[j, ] <- mvnfast::rmvn(1, mu = as.vector(M_j), sigma = V_j)
}
```

Three redundancies eliminated:

1. **`which(Z[, j] == 1)` called m×n_iter times.** Z is constant; group membership indices
   are now pre-computed once before the sampler loop as `group_indices`.

2. **`chol2inv(chol(Phi))` recomputed m times per iteration.** `Phi_inv` is now computed
   once per iteration. Additionally, `V_j = (Phi^{-1} + n_j Σ^{-1})^{-1}` depends only on
   the group size n_j, not the group itself, so it is computed once per unique group size
   (often just one value) and cached.

3. **Per-group M_j matmul replaced by a single batched matmul per unique size.**
   `group_sums %*% (Sigma_inv %*% V_j)` computes all m posterior means at once for groups
   sharing the same V_j.

## Runtime results (1000 iterations, many_groups scenario: m=500)

| Version | many_groups time | vs upstream |
|---|---|---|
| Upstream | ~36s | — |
| Optimized (current, sequential RNG) | ~13s | **~64% faster** |
| Alternative (batched rmvn, different RNG) | ~6.7s | **~81% faster** |

## The RNG tradeoff

The remaining per-group loop calls `mvnfast::rmvn(1, ...)` m times, which is the dominant
remaining cost (~7s of the 13s total). Replacing it with a single batch call
`rmvn(m, mu=rep(0,d), sigma=chol(V_j), isChol=TRUE)` + mean addition achieves the 6.7s
runtime.

However, `mvnfast::rmvn` uses the **sitmo Threefry PRNG** (a counter-based cryptographic
RNG), seeded once per call from R's base RNG via `runif(1)`. Calling `rmvn(1,...)` m times
produces a different random sequence than calling `rmvn(m,...)` once, even under the same
`set.seed()`. This causes the entire MCMC chain to diverge from the upstream reference
because the different ψ draws propagate into all subsequent updates (Φ, Σ, β, W).

Example divergence under the batched approach (baseline scenario, 1000 iterations):

```
baseline / normbeta / MMLN / beta_chain     max|diff| = 4.67e-01
baseline / normbeta / MMLN / sigma_chain    max|diff| = 1.67e-01
baseline / normbeta / MMLN / w_chain        max|diff| = 5.49e+00
baseline / normbeta / MMLN / mhaccept_chain max|diff| = 5.20e-02
baseline / normbeta / MMLN / phi_chain      max|diff| = 7.21e-01
baseline / normbeta / MMLN / psi_chain      max|diff| = 1.30e+00
```

These are not algorithmic errors — the chains target the correct posterior — but they are
different sample paths and cannot be compared against the upstream fixtures without
regenerating them.

The current implementation accepts the ~64% speedup and preserves the upstream RNG stream.
To get the additional ~17% (reaching ~81% total), regenerate fixtures from the fork code.
