library(MMLN)
library(mvnfast)

# Tests for update_psi_cpp (random_intercepts.cpp)
# Acceptance criteria:
#   1. Deterministic: C++ posterior means M_j and covariances V_j match R to 1e-8
#   2. Group index pre-computation correct for unequal group sizes (no mismatch / no reordering)
#   3. Stochastic: empirical mean of C++ psi draws converges to M_j (no outlier group beyond tolerance)

# ---- Reference computations (deterministic) ----

# R reference (solve-form): compute posterior M_j and V_j without drawing.
compute_MV_R <- function(Z, R_tot, Phi, Sigma_inv) {
  m <- ncol(Z)
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

# Exact upstream (main repo) reference for V_j:
#   V_j <- chol2inv(chol(chol2inv(chol(Phi)) + n_j * Sigma_inv))
compute_MV_R_gerber <- function(Z, R_tot, Phi, Sigma_inv) {
  m <- ncol(Z)
  M_list <- vector("list", m)
  V_list <- vector("list", m)
  for (j in seq_len(m)) {
    idx <- which(Z[, j] == 1)
    R_j <- R_tot[idx, , drop = FALSE]
    V_j <- chol2inv(chol(chol2inv(chol(Phi)) + length(idx) * Sigma_inv))
    M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
    M_list[[j]] <- as.vector(M_j)
    V_list[[j]] <- V_j
  }
  list(M = M_list, V = V_list)
}

idx_check_group_idx <- function(Z, group_idx) {
  m <- ncol(Z)
  ok <- logical(m)
  n_j <- integer(m)
  for (j in seq_len(m)) {
    # Expected 0-based indices for equality with build_group_idx.
    idx_r <- sort(as.integer(which(Z[, j] == 1) - 1L))
    idx_c <- sort(as.integer(group_idx[[j]]))
    ok[j] <- length(idx_r) == length(idx_c) && all(idx_r == idx_c)
    n_j[j] <- length(idx_r)
  }
  list(ok = ok, n_j = n_j)
}

det_diff_table <- function(cpp_moments, ref_moments) {
  m <- length(ref_moments$M)
  M_max <- numeric(m)
  V_max <- numeric(m)
  V_frob <- numeric(m)
  d <- length(ref_moments$M[[1]])
  for (j in seq_len(m)) {
    M_max[j] <- max(abs(as.vector(cpp_moments$M[[j]]) - ref_moments$M[[j]]))
    Vdiff <- as.matrix(cpp_moments$V[[j]]) - ref_moments$V[[j]]
    V_max[j] <- max(abs(Vdiff))
    V_frob[j] <- sqrt(sum(Vdiff^2))
  }
  data.frame(
    group = seq_len(m),
    d = d,
    M_max_diff = M_max,
    V_max_diff = V_max,
    V_frob_norm = V_frob
  )
}

# ---- Stochastic helper (C++ draws) ----

# Draw n_draws samples from C++ and return empirical mean.
compute_psi_empirical_mean <- function(group_idx, R_tot, Phi, Sigma_inv, n_draws = 10000) {
  m <- length(group_idx)
  d <- ncol(R_tot)
  psi_accum <- matrix(0, m, d)
  for (i in seq_len(n_draws)) {
    psi_accum <- psi_accum + update_psi_cpp(group_idx, R_tot, Phi, Sigma_inv)
  }
  psi_accum / n_draws
}

# ---- Controls ----

n_draws <- as.integer(Sys.getenv("MMLN_NDRAWS", "10000"))
emp_tol <- 0.05
det_tol <- 1e-8

# ---- Scenario 1: Equal group sizes ----

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

group_idx_check <- idx_check_group_idx(Z, group_idx)
cat("Group idx pre-computation (equal groups):", all(group_idx_check$ok), "\n")
if (!all(group_idx_check$ok)) {
  stop("FAIL: group_idx pre-computation mismatch (equal groups)")
}

ref <- compute_MV_R(Z, R_tot, Phi, Sigma_inv)
ref_gerber <- compute_MV_R_gerber(Z, R_tot, Phi, Sigma_inv)

# Sanity: the two reference implementations should agree (deterministically).
for (j in seq_len(m)) {
  Mdiff <- max(abs(ref$M[[j]] - ref_gerber$M[[j]]))
  Vdiff <- max(abs(ref$V[[j]] - ref_gerber$V[[j]]))
  stopifnot("FAIL: R reference mismatch between solve-form and Gerber-form M_j" = Mdiff < 1e-10)
  stopifnot("FAIL: R reference mismatch between solve-form and Gerber-form V_j" = Vdiff < 1e-10)
}

cat("\nTest 0 (equal groups): Deterministic diffs (C++ vs Gerber-form reference)\n")
cpp_moments <- compute_psi_moments_cpp(group_idx, R_tot, Phi, Sigma_inv)
det_tbl <- det_diff_table(cpp_moments, ref_gerber)
print(det_tbl, row.names = FALSE)
cat(sprintf(
  "\nOverall max M_j diff = %.2e; overall max V_j diff = %.2e\n",
  max(det_tbl$M_max_diff), max(det_tbl$V_max_diff)
))
stopifnot("FAIL: deterministic M_j outlier exceeds 1e-8" = max(det_tbl$M_max_diff) < det_tol)
stopifnot("FAIL: deterministic V_j outlier exceeds 1e-8" = max(det_tbl$V_max_diff) < det_tol)

cat("  PASS\n\n")

cat(sprintf("Test 1 (equal groups): Empirical mean diffs (n_draws=%d)\n", n_draws))
set.seed(1)
emp_mean <- compute_psi_empirical_mean(group_idx, R_tot, Phi, Sigma_inv, n_draws = n_draws)
emp_diffs <- numeric(m)
for (j in seq_len(m)) {
  emp_diffs[j] <- max(abs(emp_mean[j, ] - ref$M[[j]]))
}
emp_tbl <- data.frame(group = seq_len(m), n = group_sizes, max_abs_mean_diff = emp_diffs)
print(emp_tbl, row.names = FALSE)
cat(sprintf("Overall max empirical mean diff = %.4f\n", max(emp_diffs)))
stopifnot("FAIL: empirical mean outlier exceeds tolerance" = max(emp_diffs) < emp_tol)
cat("  PASS\n\n")

# ---- Scenario 2: Unequal group sizes ----

cat("Test 2: Unequal group sizes\n")
set.seed(42)
group_sizes2 <- c(3, 7, 15, 1, 24)
N2 <- sum(group_sizes2)
Z2 <- model.matrix(~ factor(rep(seq_len(m), group_sizes2)) - 1)
R_tot2 <- matrix(rnorm(N2 * d), N2, d)
group_idx2 <- build_group_idx(Z2)

group_idx_check2 <- idx_check_group_idx(Z2, group_idx2)
cat("Group idx pre-computation (unequal groups):", all(group_idx_check2$ok), "\n")
if (!all(group_idx_check2$ok)) {
  stop("FAIL: group_idx pre-computation mismatch (unequal groups)")
}

ref2 <- compute_MV_R(Z2, R_tot2, Phi, Sigma_inv)
ref2_gerber <- compute_MV_R_gerber(Z2, R_tot2, Phi, Sigma_inv)

cat("\nDeterministic diffs (C++ vs Gerber-form reference) - unequal groups:\n")
cpp_moments2 <- compute_psi_moments_cpp(group_idx2, R_tot2, Phi, Sigma_inv)
det_tbl2 <- det_diff_table(cpp_moments2, ref2_gerber)
print(det_tbl2, row.names = FALSE)
stopifnot("FAIL: deterministic M_j outlier exceeds 1e-8 (unequal groups)" = max(det_tbl2$M_max_diff) < det_tol)
stopifnot("FAIL: deterministic V_j outlier exceeds 1e-8 (unequal groups)" = max(det_tbl2$V_max_diff) < det_tol)

set.seed(2)
emp_mean2 <- compute_psi_empirical_mean(group_idx2, R_tot2, Phi, Sigma_inv, n_draws = n_draws)
emp_diffs2 <- numeric(m)
for (j in seq_len(m)) {
  emp_diffs2[j] <- max(abs(emp_mean2[j, ] - ref2$M[[j]]))
}
emp_tbl2 <- data.frame(group = seq_len(m), n = group_sizes2, max_abs_mean_diff = emp_diffs2)
cat(sprintf("\nEmpirical mean diffs (n_draws=%d) - unequal groups:\n", n_draws))
print(emp_tbl2, row.names = FALSE)
cat(sprintf("Overall max empirical mean diff = %.4f\n", max(emp_diffs2)))
stopifnot("FAIL: empirical mean outlier exceeds tolerance (unequal groups)" = max(emp_diffs2) < emp_tol)

cat("  PASS\n\n")

# ---- Scenario 3: Single-observation edge case ----

cat("Test 3: Single-observation group edge case\n")
set.seed(42)
group_sizes3 <- c(1, 9, 5, 3, 2)
N3 <- sum(group_sizes3)
Z3 <- model.matrix(~ factor(rep(seq_len(m), group_sizes3)) - 1)
R_tot3 <- matrix(rnorm(N3 * d), N3, d)
group_idx3 <- build_group_idx(Z3)

group_idx_check3 <- idx_check_group_idx(Z3, group_idx3)
cat("Group idx pre-computation (single-observation edge case):", all(group_idx_check3$ok), "\n")
if (!all(group_idx_check3$ok)) {
  stop("FAIL: group_idx pre-computation mismatch (single-observation edge case)")
}

ref3 <- compute_MV_R(Z3, R_tot3, Phi, Sigma_inv)
ref3_gerber <- compute_MV_R_gerber(Z3, R_tot3, Phi, Sigma_inv)

cat("\nDeterministic diffs (C++ vs Gerber-form reference) - single-observation:\n")
cpp_moments3 <- compute_psi_moments_cpp(group_idx3, R_tot3, Phi, Sigma_inv)
det_tbl3 <- det_diff_table(cpp_moments3, ref3_gerber)
print(det_tbl3, row.names = FALSE)
stopifnot("FAIL: deterministic M_j outlier exceeds 1e-8 (single-observation)" = max(det_tbl3$M_max_diff) < det_tol)
stopifnot("FAIL: deterministic V_j outlier exceeds 1e-8 (single-observation)" = max(det_tbl3$V_max_diff) < det_tol)

set.seed(3)
emp_mean3 <- compute_psi_empirical_mean(group_idx3, R_tot3, Phi, Sigma_inv, n_draws = n_draws)
emp_diffs3 <- numeric(m)
for (j in seq_len(m)) {
  emp_diffs3[j] <- max(abs(emp_mean3[j, ] - ref3$M[[j]]))
}
emp_tbl3 <- data.frame(group = seq_len(m), n = group_sizes3, max_abs_mean_diff = emp_diffs3)
cat(sprintf("\nEmpirical mean diffs (n_draws=%d) - single-observation:\n", n_draws))
print(emp_tbl3, row.names = FALSE)
cat(sprintf("Overall max empirical mean diff = %.4f\n", max(emp_diffs3)))
stopifnot("FAIL: empirical mean outlier exceeds tolerance (single-observation)" = max(emp_diffs3) < emp_tol)

cat("  PASS\n\n")

cat("All tests passed.\n")
