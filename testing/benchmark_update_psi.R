library(MMLN)
library(mvnfast)
library(microbenchmark)

# Benchmark: R loop vs C++ for random intercept update.
# group_idx pre-computed once per dataset before timing -- mirrors real MCMC usage.

update_psi_R <- function(Z, R_tot, Phi, Sigma_inv) {
  m <- ncol(Z)
  d <- ncol(R_tot)
  psi <- matrix(0, m, d)
  Phi_inv <- chol2inv(chol(Phi))
  for (j in seq_len(m)) {
    idx <- which(Z[, j] == 1)
    R_j <- R_tot[idx, , drop = FALSE]
    V_j <- solve(Phi_inv + length(idx) * Sigma_inv)
    M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
    psi[j, ] <- mvnfast::rmvn(1, mu = as.vector(M_j), sigma = V_j)
  }
  psi
}

# ---- Dataset 1: Simulated scale (m=20, N=400, d=2) ----
set.seed(42)
m1 <- 20; d1 <- 2
group_sizes1 <- rep(20, m1)
N1 <- sum(group_sizes1)
Z1 <- model.matrix(~ factor(rep(seq_len(m1), group_sizes1)) - 1)
R_tot1 <- matrix(rnorm(N1 * d1), N1, d1)
Phi1 <- diag(d1) * 0.5
Sigma_inv1 <- diag(d1) * 2.0
group_idx1 <- build_group_idx(Z1)

cat("Dataset 1: Simulated scale (m=20, N=400, d=2)\n")
bm1 <- microbenchmark(
  R   = update_psi_R(Z1, R_tot1, Phi1, Sigma_inv1),
  Cpp = update_psi_cpp(group_idx1, R_tot1, Phi1, Sigma_inv1),
  times = 500
)
print(bm1)
cat(sprintf("Speedup: %.1fx\n\n",
  median(bm1$time[bm1$expr == "R"]) / median(bm1$time[bm1$expr == "Cpp"])))

# ---- Dataset 2: Baseball scale (m=500, N~6250, d=3) ----
set.seed(42)
m2 <- 500; d2 <- 3
group_sizes2 <- sample(5:20, m2, replace = TRUE)
N2 <- sum(group_sizes2)
Z2 <- model.matrix(~ factor(rep(seq_len(m2), group_sizes2)) - 1)
R_tot2 <- matrix(rnorm(N2 * d2), N2, d2)
Phi2 <- diag(d2) * 0.5
Sigma_inv2 <- diag(d2) * 2.0
group_idx2 <- build_group_idx(Z2)

cat(sprintf("Dataset 2: Baseball scale (m=500, N=%d, d=3)\n", N2))
bm2 <- microbenchmark(
  R   = update_psi_R(Z2, R_tot2, Phi2, Sigma_inv2),
  Cpp = update_psi_cpp(group_idx2, R_tot2, Phi2, Sigma_inv2),
  times = 100
)
print(bm2)
cat(sprintf("Speedup: %.1fx\n\n",
  median(bm2$time[bm2$expr == "R"]) / median(bm2$time[bm2$expr == "Cpp"])))

# ---- Dataset 3: Unequal groups stress test (m=100, mixed sizes, d=2) ----
set.seed(42)
m3 <- 100; d3 <- 2
group_sizes3 <- c(rep(1, 10), rep(5, 40), rep(50, 50))
N3 <- sum(group_sizes3)
Z3 <- model.matrix(~ factor(rep(seq_len(m3), group_sizes3)) - 1)
R_tot3 <- matrix(rnorm(N3 * d3), N3, d3)
Phi3 <- diag(d3) * 0.5
Sigma_inv3 <- diag(d3) * 2.0
group_idx3 <- build_group_idx(Z3)

cat(sprintf("Dataset 3: Unequal groups (m=100, N=%d, mix of n=1/5/50, d=2)\n", N3))
bm3 <- microbenchmark(
  R   = update_psi_R(Z3, R_tot3, Phi3, Sigma_inv3),
  Cpp = update_psi_cpp(group_idx3, R_tot3, Phi3, Sigma_inv3),
  times = 200
)
print(bm3)
cat(sprintf("Speedup: %.1fx\n",
  median(bm3$time[bm3$expr == "R"]) / median(bm3$time[bm3$expr == "Cpp"])))
