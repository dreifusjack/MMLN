# test_update_psi.R
# Acceptance tests for the C++ random-intercept update.
# Run after devtools::load_all() to ensure the package is compiled.
#
# Tests:
#   0 - Core math: compute_psi_moments_cpp M_j / V_j match R reference to 1e-8
#   1 - Equal groups: empirical mean over 10,000 draws converges to M_j (tol 0.05)
#   2 - Unequal groups: same convergence test with group sizes c(3,7,15,1,24)
#   3 - Single-obs group: group with n=1 handled correctly

library(MMLN)

cat("=== test_update_psi.R ===\n\n")

# ---- helpers ----------------------------------------------------------------

# R reference: compute posterior mean and covariance for one group
r_psi_moments <- function(R_j, Phi, Sigma_inv) {
  Phi_inv <- chol2inv(chol(Phi))
  n_j     <- nrow(R_j)
  A_j     <- Phi_inv + n_j * Sigma_inv
  V_j     <- chol2inv(chol(A_j))
  s_j     <- colSums(R_j)
  M_j     <- V_j %*% (Sigma_inv %*% s_j)
  list(M = as.vector(M_j), V = V_j)
}

pass <- function(name) cat(sprintf("  PASS  %s\n", name))
fail <- function(name, msg) { cat(sprintf("  FAIL  %s: %s\n", name, msg)); stop(msg) }

check <- function(name, cond, msg = "condition not met") {
  if (isTRUE(cond)) pass(name) else fail(name, msg)
}

# ---- shared fixtures --------------------------------------------------------

set.seed(42)
d <- 3   # ALR dimension
N <- 50
m <- 4   # number of groups (equal size = 50/4 → actually use explicit Z below)

# ---- TEST 0: core math -------------------------------------------------------

cat("Test 0: core math (compute_psi_moments_cpp matches R to 1e-8)\n")

# Build a balanced Z (equal groups of 10 each, m=5, N=50)
m0 <- 5; n_per <- 10; N0 <- m0 * n_per
Z0 <- matrix(0, N0, m0)
for (j in seq_len(m0)) Z0[(j-1)*n_per + seq_len(n_per), j] <- 1

Phi0      <- diag(d) + 0.3 * matrix(1, d, d); Phi0 <- (Phi0 + t(Phi0)) / 2
Sigma0    <- diag(d) * 0.5
Sigma_inv0 <- chol2inv(chol(Sigma0))
R_tot0    <- matrix(rnorm(N0 * d), N0, d)

group_idx0 <- build_group_idx(Z0)
cpp_moments <- compute_psi_moments_cpp(group_idx0, R_tot0, Phi0, Sigma_inv0)

for (j in seq_len(m0)) {
  idx_r <- which(Z0[, j] == 1)
  ref   <- r_psi_moments(R_tot0[idx_r, , drop = FALSE], Phi0, Sigma_inv0)
  cpp_M <- as.vector(cpp_moments$M[[j]])
  cpp_V <- cpp_moments$V[[j]]
  check(sprintf("Test0 group%d M", j), max(abs(cpp_M - ref$M)) < 1e-8,
        sprintf("max|M diff|=%.2e", max(abs(cpp_M - ref$M))))
  check(sprintf("Test0 group%d V", j), max(abs(cpp_V - ref$V)) < 1e-8,
        sprintf("max|V diff|=%.2e", max(abs(cpp_V - ref$V))))
}

# ---- TEST 1: equal groups, convergence of empirical mean -------------------

cat("\nTest 1: equal groups — empirical mean converges to M_j (tol 0.05)\n")

n_draws <- 10000
# Reuse fixtures from Test 0
psi_draws <- array(0, dim = c(n_draws, m0, d))
for (s in seq_len(n_draws)) {
  psi_draws[s, , ] <- update_psi_cpp(group_idx0, R_tot0, Phi0, Sigma_inv0)
}

emp_means <- apply(psi_draws, c(2, 3), mean)  # m0 x d

for (j in seq_len(m0)) {
  idx_r <- which(Z0[, j] == 1)
  ref   <- r_psi_moments(R_tot0[idx_r, , drop = FALSE], Phi0, Sigma_inv0)
  err   <- max(abs(emp_means[j, ] - ref$M))
  check(sprintf("Test1 group%d emp_mean", j), err < 0.05,
        sprintf("max|emp_mean - M_j|=%.4f", err))
}

# ---- TEST 2: unequal groups -------------------------------------------------

cat("\nTest 2: unequal groups — empirical mean converges to M_j (tol 0.05)\n")

grp_sizes <- c(3, 7, 15, 1, 24)
m2 <- length(grp_sizes); N2 <- sum(grp_sizes)
Z2 <- matrix(0, N2, m2)
row_start <- 1
for (j in seq_len(m2)) {
  Z2[row_start:(row_start + grp_sizes[j] - 1), j] <- 1
  row_start <- row_start + grp_sizes[j]
}

set.seed(7)
Phi2      <- diag(d) * 2
Sigma2    <- diag(d) * 0.8; Sigma2[1,2] <- Sigma2[2,1] <- 0.2
Sigma_inv2 <- chol2inv(chol(Sigma2))
R_tot2    <- matrix(rnorm(N2 * d), N2, d)
group_idx2 <- build_group_idx(Z2)

psi_draws2 <- array(0, dim = c(n_draws, m2, d))
for (s in seq_len(n_draws)) {
  psi_draws2[s, , ] <- update_psi_cpp(group_idx2, R_tot2, Phi2, Sigma_inv2)
}
emp_means2 <- apply(psi_draws2, c(2, 3), mean)

for (j in seq_len(m2)) {
  idx_r <- which(Z2[, j] == 1)
  ref   <- r_psi_moments(R_tot2[idx_r, , drop = FALSE], Phi2, Sigma_inv2)
  err   <- max(abs(emp_means2[j, ] - ref$M))
  check(sprintf("Test2 group%d (n=%d) emp_mean", j, grp_sizes[j]), err < 0.05,
        sprintf("max|emp_mean - M_j|=%.4f", err))
}

# ---- TEST 3: single-obs group -----------------------------------------------

cat("\nTest 3: single-observation group handled correctly\n")

# group 4 in Test 2 has n=1; verify moments match R reference
j_single <- 4
idx_single <- which(Z2[, j_single] == 1)  # length 1
ref_single <- r_psi_moments(R_tot2[idx_single, , drop = FALSE], Phi2, Sigma_inv2)
cpp_m3 <- compute_psi_moments_cpp(group_idx2, R_tot2, Phi2, Sigma_inv2)

M_single <- as.vector(cpp_m3$M[[j_single]])
V_single <- cpp_m3$V[[j_single]]
check("Test3 single-obs M", max(abs(M_single - ref_single$M)) < 1e-8,
      sprintf("max|M diff|=%.2e", max(abs(M_single - ref_single$M))))
check("Test3 single-obs V", max(abs(V_single - ref_single$V)) < 1e-8,
      sprintf("max|V diff|=%.2e", max(abs(V_single - ref_single$V))))

# also check draw does not error or produce NaN/Inf
psi_single_draws <- replicate(200, update_psi_cpp(group_idx2, R_tot2, Phi2, Sigma_inv2)[j_single, ])
check("Test3 single-obs no NaN", !any(is.nan(psi_single_draws)), "NaN in draws")
check("Test3 single-obs no Inf", !any(is.infinite(psi_single_draws)), "Inf in draws")

cat("\n=== All tests passed ===\n")
