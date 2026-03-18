library(MMLN)
library(mvnfast)

# Tests for update_psi_cpp (random_intercepts.cpp)
# Acceptance criteria:
#   1. C++ posterior mean M_j converges to R reference within 0.05 (distributional)
#   2. Group index pre-computation correct for unequal group sizes
#   3. Handles a group with only one observation

# R reference: compute posterior M_j and V_j without drawing (deterministic)
compute_MV_R <- function(Z, R_tot, Phi, Sigma_inv) {
  m <- ncol(Z)
  d <- ncol(R_tot)
  Phi_inv <- chol2inv(chol(Phi))
  M_list <- vector("list", m)
  V_list <- vector("list", m)
  for (j in seq_len(m)) {
    idx <- which(Z[, j] == 1)
    R_j <- R_tot[idx, , drop = FALSE]
    V_j <- solve(Phi_inv + length(idx) * Sigma_inv)
    M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
    M_list[[j]] <- as.vector(M_j)
    V_list[[j]] <- V_j
  }
  list(M = M_list, V = V_list)
}

# Draw n_draws samples from C++ and return empirical mean.
# group_idx pre-computed once and reused across all draws.
compute_psi_empirical_mean <- function(group_idx, R_tot, Phi, Sigma_inv, n_draws = 10000) {
  m <- length(group_idx)
  d <- ncol(R_tot)
  psi_accum <- matrix(0, m, d)
  for (i in seq_len(n_draws)) {
    psi_accum <- psi_accum + update_psi_cpp(group_idx, R_tot, Phi, Sigma_inv)
  }
  psi_accum / n_draws
}

# ---- Shared inputs ----
set.seed(42)
m <- 5
d <- 2
group_sizes <- c(10, 10, 10, 10, 10)
N <- sum(group_sizes)
Z <- model.matrix(~ factor(rep(seq_len(m), group_sizes)) - 1)
R_tot <- matrix(rnorm(N * d), N, d)
Phi <- diag(d) * 0.5
Sigma_inv <- diag(d) * 2.0
group_idx <- build_group_idx(Z)
ref <- compute_MV_R(Z, R_tot, Phi, Sigma_inv)

# ---- Test 0: Core math -- M_j and V_j match R to 1e-8 (deterministic) ----
# This directly verifies the posterior mean and covariance formulas are correct.
# No randomness involved. A failure here means wrong math, not just sampling noise.
cat("Test 0: Core math -- M_j and V_j match R reference to 1e-8\n")
cpp_moments <- compute_psi_moments_cpp(group_idx, R_tot, Phi, Sigma_inv)
for (j in seq_len(m)) {
  M_diff <- max(abs(as.vector(cpp_moments$M[[j]]) - ref$M[[j]]))
  V_diff <- max(abs(as.matrix(cpp_moments$V[[j]]) - ref$V[[j]]))
  cat(sprintf("  Group %d: M_j max diff = %.2e, V_j max diff = %.2e\n", j, M_diff, V_diff))
  stopifnot("FAIL: C++ M_j does not match R to 1e-8" = M_diff < 1e-8)
  stopifnot("FAIL: C++ V_j does not match R to 1e-8" = V_diff < 1e-8)
}
cat("  PASS\n\n")

# ---- Test 1: Equal groups ----
cat("Test 1: Posterior M_j convergence (equal groups)\n")
set.seed(1)
emp_mean <- compute_psi_empirical_mean(group_idx, R_tot, Phi, Sigma_inv)
for (j in seq_len(m)) {
  diff <- max(abs(emp_mean[j, ] - ref$M[[j]]))
  cat(sprintf("  Group %d: max diff = %.4f\n", j, diff))
  stopifnot("FAIL: empirical mean deviates from M_j by more than 0.05" = diff < 0.05)
}
cat("  PASS\n\n")

# ---- Test 2: Unequal group sizes ----
cat("Test 2: Unequal group sizes\n")
set.seed(42)
group_sizes2 <- c(3, 7, 15, 1, 24)
N2 <- sum(group_sizes2)
Z2 <- model.matrix(~ factor(rep(seq_len(m), group_sizes2)) - 1)
R_tot2 <- matrix(rnorm(N2 * d), N2, d)
group_idx2 <- build_group_idx(Z2)
ref2 <- compute_MV_R(Z2, R_tot2, Phi, Sigma_inv)
set.seed(2)
emp_mean2 <- compute_psi_empirical_mean(group_idx2, R_tot2, Phi, Sigma_inv)
for (j in seq_len(m)) {
  diff <- max(abs(emp_mean2[j, ] - ref2$M[[j]]))
  cat(sprintf("  Group %d (n=%d): max diff = %.4f\n", j, group_sizes2[j], diff))
  stopifnot("FAIL: unequal group empirical mean deviates from M_j" = diff < 0.05)
}
cat("  PASS\n\n")

# ---- Test 3: Single-observation group edge case ----
cat("Test 3: Single-observation group\n")
set.seed(42)
group_sizes3 <- c(1, 9, 5, 3, 2)
N3 <- sum(group_sizes3)
Z3 <- model.matrix(~ factor(rep(seq_len(m), group_sizes3)) - 1)
R_tot3 <- matrix(rnorm(N3 * d), N3, d)
group_idx3 <- build_group_idx(Z3)
ref3 <- compute_MV_R(Z3, R_tot3, Phi, Sigma_inv)
set.seed(3)
emp_mean3 <- compute_psi_empirical_mean(group_idx3, R_tot3, Phi, Sigma_inv)
for (j in seq_len(m)) {
  diff <- max(abs(emp_mean3[j, ] - ref3$M[[j]]))
  cat(sprintf("  Group %d (n=%d): max diff = %.4f\n", j, group_sizes3[j], diff))
  stopifnot("FAIL: single-obs group empirical mean deviates from M_j" = diff < 0.05)
}
cat("  PASS\n\n")

cat("All tests passed.\n")
