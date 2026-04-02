# Normbeta Numerical Divergence: `large_coefs` Scenario

## Summary

The C++ `mh_update_cpp` produces divergent chains compared to the upstream R
implementation (`normbetapropdist` / `normbetaloglike`) on the `large_coefs`
scenario only. All other scenarios (9/10) are bitwise identical or within
tolerance.

| Chain          | max\|diff\| |
|----------------|-------------|
| w_chain        | 1.97e+01    |
| sigma_chain    | 6.48e+00    |
| beta_chain     | 8.78e-01    |
| mhaccept_chain | 3.33e-02    |

## Root Cause

The alpha parameter in the normbeta proposal is computed as:

```
alpha = (1 + exp(mu)) / sigma_ii - exp(mu) / (1 + exp(mu))
```

The `large_coefs` scenario has beta coefficients in [-4, 4] with `Sigma = 2I`,
so `mu = X * beta` can reach values where `exp(mu)` overflows to `Inf`. When
that happens:

- `(1 + Inf) / sigma_ii` = `Inf`
- `Inf / (1 + Inf)` = `Inf / Inf` = **`NaN`**
- `alpha` = `Inf - NaN` = **`NaN`**

This NaN propagates through `digamma`/`trigamma` into the proposal, which gets
rejected (`ratio[is.na(ratio)] <- -Inf`). If the R and C++ code paths hit NaN
at even slightly different thresholds (due to intermediate rounding), one
accepts a proposal while the other rejects it, **permanently desynchronizing the
RNG stream**. Every subsequent iteration diverges — not from accumulated
floating-point drift, but from a single flipped accept/reject decision.

### Additional difference: safety floor

The C++ code adds a floor not present in upstream R:

```cpp
if (alpha < 1e-3) alpha = 1e-3;  // C++ only — not in R normbetapropdist
```

This was borrowed from `betapropdist` (the exact Beta proposal) but the R
`normbetapropdist` never had it. It can mask NaN and produce different proposal
values even when overflow hasn't occurred.

## Proposed Fix: Numerically Stable Alpha Formula

Rewrite the alpha computation to avoid `exp(mu)` overflow entirely. The
problematic term `exp(mu) / (1 + exp(mu))` is mathematically the sigmoid
function, which has a stable form:

```
sigmoid(mu) = 1 / (1 + exp(-mu))
```

And `(1 + exp(mu)) / sigma_ii` can be rewritten using the log-sum-exp identity:

```
1 + exp(mu) = exp(mu) * (1 + exp(-mu))
```

So the full alpha formula becomes:

```
alpha = exp(mu) * (1 + exp(-mu)) / sigma_ii  -  1 / (1 + exp(-mu))
      = exp(mu) / sigma_ii  +  1/sigma_ii  -  sigmoid(mu)
```

Or more directly, factoring through `sigmoid(mu) = exp(mu) / (1 + exp(mu))`:

```
alpha = (1 + exp(mu)) / sigma_ii - sigmoid(mu)
      = sigmoid(mu) * ((1 + exp(mu))^2 / (exp(mu) * sigma_ii) - 1)
```

A practical stable implementation:

```
sig   = 1 / (1 + exp(-mu))          # sigmoid, always in (0,1)
alpha = sig * (1/sig + exp(mu)) / sigma_ii - sig
      = (1 + exp(mu)) / sigma_ii - sig
```

The key insight: replace `exp(mu) / (1 + exp(mu))` with `1 / (1 + exp(-mu))`.
For large positive mu, `exp(-mu) ≈ 0` so this evaluates cleanly to ~1. The
`Inf/Inf` indeterminate form is eliminated.

The same rewrite should be applied to:
- `mh_update.cpp` (C++ — the fork)
- `normbetapropdist` in `mln_helpers.R` (R — should be proposed upstream)
- `normbetaloglike` in `mln_helpers.R` (R — should be proposed upstream)

With this fix, the safety floor becomes unnecessary and should be removed.

## Alternatives

**Accept the divergence.** The `large_coefs` scenario deliberately targets
overflow edge cases. We could increase the tolerance for this scenario or
exclude it from the normbeta comparison. This is the simplest option but
leaves a latent numerical bug in production.

**Clamp mu.** Cap `mu` at ±500 before computing `exp(mu)`. Both R and C++
would agree in the overflow region. Simpler than rewriting the formula but
less principled — it changes the statistical behavior for extreme inputs.

## Recommendation

Apply the stable formula rewrite (option A) to both C++ and R. This:
- Eliminates the overflow that causes NaN
- Removes the need for the safety floor
- Makes C++ and R agree for all inputs, not just non-extreme ones
- Is a correctness improvement to the math, not a workaround
