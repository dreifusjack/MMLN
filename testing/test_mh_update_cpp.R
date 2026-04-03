library(MMLN)
library(mvnfast)

# Test numerical equivalence of mh_update_cpp against the pure-R normbeta logic.
#
# This is the exact code path used when proposal == "normbeta" in FMLN().
# We compare:
#   - proposed W_new
#   - log_q_old and log_q_new
#
# For deterministic comparison we set the same RNG seed and rely on mh_update_cpp()
# using R's RNG for its proposal draws.

set.seed(1)
N <- 5
d <- 2
mh_scale <- 1.0

# Use diagonal Sigma (both R and C++ only access Sigma[i,i] here)
Sigma <- diag(d)
Sigma_scaled <- mh_scale * Sigma

Y <- matrix(sample(1:20, N * (d + 1), replace = TRUE), nrow = N, ncol = d + 1)
Mu <- matrix(rnorm(N * d), nrow = N, ncol = d)
W <- matrix(rnorm(N * d), nrow = N, ncol = d)

ymu <- cbind(Y, Mu) # N x (d+1+d)

cat("Test mh_update_cpp vs R normbeta proposal/loglike\n")

# 1) R reference
set.seed(123)
W_new_R <- t(apply(ymu, 1, normbetapropdist, Sigma = Sigma_scaled))

log_q_old_R <- apply(cbind(ymu, W), 1, normbetaloglike, Sigma = Sigma_scaled)
log_q_new_R <- apply(cbind(ymu, W_new_R), 1, normbetaloglike, Sigma = Sigma_scaled)

# 2) C++ implementation (same RNG seed => same z draws)
set.seed(123)
out_cpp <- MMLN:::mh_update_cpp(W, Y, Mu, Sigma_scaled)
W_new_C <- out_cpp$W_new
log_q_old_C <- out_cpp$log_q_old
log_q_new_C <- out_cpp$log_q_new

# 3) Compare
diff_W <- max(abs(W_new_C - W_new_R))
diff_old <- max(abs(log_q_old_C - log_q_old_R))
diff_new <- max(abs(log_q_new_C - log_q_new_R))

cat(sprintf("  max|W_new_cpp - W_new_R|   = %.3e\n", diff_W))
cat(sprintf("  max|log_q_old_cpp - R|     = %.3e\n", diff_old))
cat(sprintf("  max|log_q_new_cpp - R|     = %.3e\n", diff_new))

stopifnot(diff_W < 1e-10)
stopifnot(diff_old < 1e-10)
stopifnot(diff_new < 1e-10)

cat("PASS\n")

