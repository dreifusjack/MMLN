# Cross-Platform Chain Divergence: BLAS/LAPACK Non-Determinism

## Summary

MMLN chains produced on Linux and macOS diverge even with identical seeds and
code. This is **not a bug** — it is an inherent consequence of different
BLAS/LAPACK implementations across operating systems.

## Root Cause

Both `mvnfast::rmvn` and base R's `rWishart` call LAPACK's `dpotrf` (Cholesky
decomposition) internally. Linux typically uses **OpenBLAS** while macOS uses
**Apple Accelerate**. These produce Cholesky factors that differ at machine
epsilon (~1e-16).

The RNG streams themselves are cross-platform deterministic:
- `mvnfast` uses sitmo threefry (counter-based, integer-only arithmetic)
- `rWishart` uses R's Mersenne-Twister

The divergence comes from the LAPACK layer, not the random number generators.

## Why MMLN Amplifies It

MMLN computes Cholesky decompositions on **stochastic matrices** (`Sigma`,
`Phi`, `V_j`) that change every MCMC iteration. Each iteration's ~1e-16
rounding difference feeds into the next iteration's inputs, compounding over
1000 iterations into meaningful divergence.

| MCMC step | LAPACK-dependent call |
|-----------|----------------------|
| W proposal | `mvnfast::rmvn` → `arma::chol(sigma)` |
| psi update | `chol(Phi)`, `chol2inv`, `mvnfast::rmvn(1, ..., sigma=V_j)` |
| Phi update | `chol(S1)`, `rWishart` → `dpotrf` |
| beta update | `mvnfast::rmvn(1, ..., sigma=cov_b)` |
| Sigma update | `chol(S2)`, `rWishart` → `dpotrf` |

Every step feeds a fresh stochastic matrix through a platform-dependent
Cholesky, so the drift accumulates multiplicatively.

## Implication for Testing

Reference fixtures and fork results **must be generated on the same machine**.
Cross-machine comparisons require tolerances, not exact equality. The test
suite already enforces this by running both `make reference` and `make test`
locally.
