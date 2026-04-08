# mmln_utils.R - Mixed Effects MLN Utilities

# Required libraries
library(mvnfast)

# Simulate mixed MLN data with one random intercept per group
simulate_mixed_mln_data <- function(m,                # number of groups
                                    n_i,              # observations per group (scalar or vector)
                                    p,                # number of fixed covariates (including intercept)
                                    d,                # number of non-baseline logits (categories = d+1)
                                    beta,             # p x d fixed-effect coefficients
                                    Sigma,            # d x d covariance for within-group errors
                                    Phi,              # d x d covariance for random intercepts
                                    PA_mean = 200     # avg total counts per obs
) {
  if(length(n_i) == 1) n_i <- rep(n_i, m)
  id <- rep(seq_len(m), times = n_i)
  N <- length(id)
  X <- cbind(1, matrix(rnorm(N * (p - 1)), nrow = N, ncol = p - 1))
  Z <- model.matrix(~ factor(id) - 1)
  psi_mat <- mvnfast::rmvn(m, mu = rep(0, d), sigma = Phi)
  PA <- round(rpois(N, PA_mean))
  W <- matrix(0, nrow = N, ncol = d + 1)
  for(j in seq_len(N)) {
    eps <- mvnfast::rmvn(1, mu = rep(0, d), sigma = Sigma)
    y_j <- X[j, , drop=FALSE] %*% beta + eps + psi_mat[id[j], ]
    exp_y <- exp(y_j)
    probs <- c(exp_y / (1 + sum(exp_y)), 1 / (1 + sum(exp_y)))
    W[j, ] <- as.vector(rmultinom(1, size = PA[j], prob = probs))
  }
  list(
    W = W,
    X = X,
    Z = Z,
    id = id,
    PA = PA,
    beta = beta,
    Sigma = Sigma,
    Phi = Phi
  )
}

# Gibbs sampler for mixed-effects MLN model with progress timer
run_mixed_gibbs_sampler <- function(W, X, Z, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = c("norm", "beta", "normbeta")) {
  proposal <- match.arg(proposal)
  N <- nrow(W); d <- ncol(W) - 1; p <- ncol(X); m <- ncol(Z)
  if(is.null(prior_settings)) prior_settings <- list(
    beta_var = 10,
    nu_S     = d + 1, Lambda_S = diag(d),
    nu_P     = d + 1, Lambda_P = diag(d)
  )
  keep <- seq(burn_in + 1, n_iter, by = thin)
  n_save <- length(keep)
  beta_chain  <- vector("list", n_save)
  sigma_chain <- vector("list", n_save)
  phi_chain   <- vector("list", n_save)
  psi_chain   <- vector("list", n_save)
  y_chain     <- vector("list", n_save)

  # initialize
  Y         <- alr(compress_counts(W))   # from mln_helpers.R
  beta      <- matrix(0, p, d)
  Sigma     <- diag(d)
  psi       <- matrix(0, m, d)
  Phi       <- diag(d)
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv  <- chol2inv(chol(crossprod(X)))

  warned_na_ratio <- FALSE
  if(verbose) {
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    start_time <- Sys.time()
  }
  save_i <- 1

  for(it in seq_len(n_iter)) {
    # current mean
    Mu  <- X %*% beta + Z %*% psi   # N x d
    wmu <- cbind(W, Mu)

    # MH update of latent Y
    if(proposal == "norm") {
      Y_prop    <- Y + mvnfast::rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- log_q_new <- rep(0, N)
    } else if(proposal == "beta") {
      P_old     <- t(apply(Y, 1, ytopstar))
      P_new     <- t(apply(wmu, 1, betapropdist, Sigma = mh_scale * Sigma))
      Y_prop    <- t(apply(P_new, 1, pstartoy))
      log_q_old <- apply(cbind(wmu, P_old), 1, betaloglike, Sigma = mh_scale * Sigma)
      log_q_new <- apply(cbind(wmu, P_new), 1, betaloglike, Sigma = mh_scale * Sigma)
    } else {
      Y_prop    <- t(apply(wmu, 1, normbetapropdist, Sigma = mh_scale * Sigma))
      log_q_old <- apply(cbind(wmu, Y), 1, normbetaloglike, Sigma = mh_scale * Sigma)
      log_q_new <- apply(cbind(wmu, Y_prop), 1, normbetaloglike, Sigma = mh_scale * Sigma)
    }
    expY   <- rowSums(exp(Y)); expYp <- rowSums(exp(Y_prop))
    ll_old <- rowSums(W[,1:d] * Y[,1:d]) - rowSums(W * log1p(expY)) -
      0.5 * rowSums((Y - Mu) %*% Sigma_inv * (Y - Mu))
    ll_new <- rowSums(W[,1:d] * Y_prop[,1:d]) - rowSums(W * log1p(expYp)) -
      0.5 * rowSums((Y_prop - Mu) %*% Sigma_inv * (Y_prop - Mu))
    ratio  <- ll_new - ll_old + log_q_new - log_q_old
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accepted <- log(runif(N)) < ratio
    Y[accepted, ] <- Y_prop[accepted, ]

    # update random intercepts psi_j
    R_tot <- Y - X %*% beta
    for(j in seq_len(m)) {
      idx <- which(Z[, j] == 1)
      R_j <- R_tot[idx, , drop = FALSE]
      V_j <- solve(chol2inv(chol(Phi)) + length(idx) * Sigma_inv)
      M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
      psi[j, ] <- mvnfast::rmvn(1, mu = as.vector(M_j), sigma = V_j)
    }

    # update Phi
    S_psi <- t(psi) %*% psi
    Phi   <- solve(rWishart(1,
                            df    = prior_settings$nu_P + m,
                            Sigma = solve(prior_settings$Lambda_P + S_psi))[,,1])

    # update beta
    R     <- Y - Z %*% psi
    beta0 <- S_xx_inv %*% (t(X) %*% R)
    cov_b <- kronecker(S_xx_inv, Sigma)
    beta  <- matrix(mvnfast::rmvn(1,
                                  mu    = as.vector(beta0),
                                  sigma = cov_b),
                    nrow = p)

    # update Sigma
    Eps   <- R - X %*% beta
    S_mat <- t(Eps) %*% Eps
    Sigma <- solve(rWishart(1,
                            df    = prior_settings$nu_S + N,
                            Sigma = solve(prior_settings$Lambda_S + S_mat))[,,1])
    Sigma_inv <- chol2inv(chol(Sigma))

    # save samples
    if(it %in% keep) {
      beta_chain[[save_i]]  <- beta
      sigma_chain[[save_i]] <- Sigma
      phi_chain[[save_i]]   <- Phi
      psi_chain[[save_i]]   <- psi
      y_chain[[save_i]]     <- Y
      save_i <- save_i + 1
    }

    if (verbose && (it %% max(1, floor(n_iter / 100)) == 0 || it == n_iter)) {
      setTxtProgressBar(pb, it)
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / it) * (n_iter - it)
      cat(sprintf("\r ETA: %s", format(.POSIXct(eta, tz = "GMT"), "%M:%S")))
      flush.console()
    }
  }

  if(verbose) close(pb)

  list(
    beta_chain  = beta_chain,
    sigma_chain = sigma_chain,
    phi_chain   = phi_chain,
    psi_chain   = psi_chain,
    y_chain     = y_chain
  )
}

# Example test for mixed-effects Gibbs sampler (unchanged)
sim <- simulate_mixed_mln_data(
  m       = 30,
  n_i     = 10,
  p       = 2,
  d       = 2,
  beta    = matrix(c(-2, 0.5, 1, -0.3), nrow = 2, ncol = 2),
  Sigma   = diag(c(1.2, .7)),
  Phi     = diag(c(0.5, 0.3)),
  PA_mean = 100
)
res <- run_mixed_gibbs_sampler(
  W       = sim$W,
  X       = sim$X,
  Z       = sim$Z,
  n_iter  = 1000,
  burn_in = 180,
  thin    = 2,
  proposal= "normbeta",
  verbose = TRUE
)
str(res, max.level = 1)


beta_samples3 <- matrix(0, nrow = length(res$beta_chain), ncol = length(sim$beta))
for(i in 1:length(res$beta_chain)){
  beta_samples3[i,] <- c(res$beta_chain[[i]])
}
plot(beta_samples3[,2])
abline(h = sim$beta[2,1])
head(beta_samples3)
sim$beta

w_preds3 <- lapply(seq_along(res$y_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res$beta_chain[[i]],
                              Sigma = res$sigma_chain[[i]],
                              PA = sim$PA,
                              Z = sim$Z,
                              psi = res$psi_chain[[i]],
                              mixed = TRUE
  )
})

mahal_resids <- compute_mahalanobis_resids(sim$W, w_preds3)
print(summary(mahal_resids))
plot(mahal_resids)
qqnorm(mahal_resids)
abline(0, 1, col = "red", lwd = 2)
beta_chain_array <- simplify2array(res$beta_chain)
plot_trace_and_summary(beta_chain_array, "beta")
sim$beta
sigma_chain_array <- simplify2array(res$sigma_chain)
plot_trace_and_summary(sigma_chain_array, "sigma")
sim$Sigma
phi_chain_array <- simplify2array(res$phi_chain)
plot_trace_and_summary(phi_chain_array, "phi")
sim$Phi
