library(mvnfast)
# working fixed effects MLN model

run_fixed_gibbs_sampler <- function(W, X, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = c("norm", "beta", "normbeta")) {
  proposal <- match.arg(proposal)
  N <- nrow(W)
  d <- ncol(W) - 1
  p <- ncol(X)
  PA <- rowSums(W)

  if (is.null(prior_settings)) {
    prior_settings <- list(
      beta_var = 10,
      nu_S = d + 1,
      Lambda_S = diag(d)
    )
  }

  beta <- matrix(0, nrow = p, ncol = d)
  Sigma <- diag(d)

  keep_iters <- seq(burn_in + 1, n_iter, by = thin)
  n_save <- length(keep_iters)

  beta_chain <- vector("list", n_save)
  sigma_chain <- vector("list", n_save)
  y_chain <- vector("list", n_save)

  warned_na_ratio <- FALSE
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  start_time <- Sys.time()
  save_idx <- 1

  Y <- alr(compress_counts(W))
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv <- chol2inv(chol(crossprod(X)))

  for (i in seq_len(n_iter)) {
    Mu <- tcrossprod(X, t(beta))
    wmu <- cbind(W, Mu)

    if (proposal == "norm") {
      Y_new <- Y + rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- rep(0, N)
      log_q_prop <- rep(0, N)
    } else if (proposal == "beta") {
      Pstar <- t(apply(Y, 1, ytopstar))
      Pstar_new <- t(apply(wmu, 1, betapropdist, Sigma = mh_scale * Sigma))
      Y_new <- t(apply(Pstar_new, 1, pstartoy))
      wmuPstar <- cbind(wmu, Pstar)
      wmuPstar_new <- cbind(wmu, Pstar_new)
      log_q_old <- apply(wmuPstar, 1, betaloglike, Sigma = mh_scale * Sigma)
      log_q_prop <- apply(wmuPstar_new, 1, betaloglike, Sigma = mh_scale * Sigma)
    } else if (proposal == "normbeta") {
      Y_new <- t(apply(wmu, 1, normbetapropdist, Sigma = mh_scale * Sigma))
      wmuy <- cbind(wmu, Y)
      wmuy_new <- cbind(wmu, Y_new)
      log_q_old <- apply(wmuy, 1, normbetaloglike, Sigma = mh_scale * Sigma)
      log_q_prop <- apply(wmuy_new, 1, normbetaloglike, Sigma = mh_scale * Sigma)
    }

    Y_diff <- Y - Mu
    Y_new_diff <- Y_new - Mu
    newnormpart <- rowSums(tcrossprod(Y_new_diff, t(Sigma_inv)) * Y_new_diff) / 2
    oldnormpart <- rowSums(tcrossprod(Y_diff, t(Sigma_inv)) * Y_diff) / 2
    sexpY <- rowSums(exp(Y))
    sexpYn <- rowSums(exp(Y_new))
    newloglike <- rowSums(W[, 1:d] * Y_new[, 1:d]) - rowSums(W * log1p(sexpYn)) - newnormpart
    oldloglike <- rowSums(W[, 1:d] * Y[, 1:d]) - rowSums(W * log1p(sexpY)) - oldnormpart

    ratio <- newloglike - oldloglike + log_q_old - log_q_prop
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accept <- log(runif(N)) < ratio
    Y[accept, ] <- Y_new[accept, ]

    R <- Y
    beta_hat <- tcrossprod(S_xx_inv, t(crossprod(X, R)))
    post_beta_vec_cov <- kronecker(S_xx_inv, Sigma)
    post_beta_vec <- rmvn(1, as.vector(beta_hat), post_beta_vec_cov)
    beta <- matrix(post_beta_vec, nrow = p)

    Sigma <- chol2inv(chol(rWishart(1, df = prior_settings$nu_S + N,
                                    Sigma = solve(prior_settings$Lambda_S + crossprod(Y - X %*% beta)))[,,1]))
    Sigma_inv <- chol2inv(chol(Sigma))

    if (i %in% keep_iters) {
      beta_chain[[save_idx]] <- beta
      sigma_chain[[save_idx]] <- Sigma
      y_chain[[save_idx]] <- Y
      save_idx <- save_idx + 1
    }

    if (verbose && (i %% max(1, floor(n_iter / 100)) == 0 || i == n_iter)) {
      setTxtProgressBar(pb, i)
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / i) * (n_iter - i)
      cat(sprintf("\r ETA: %s", format(.POSIXct(eta, tz = "GMT"), "%M:%S")))
      flush.console()
    }
  }
  if (verbose) close(pb)

  list(
    beta_chain = beta_chain,
    sigma_chain = sigma_chain,
    y_chain = y_chain
  )
}


# Test the Fixed Gibbs sampler wrapper
# Using data simulated from sam_gerber_craig_code.R
# Under the True Distribution
sampler_result_fix <- run_fixed_gibbs_sampler(as.matrix(ml2[,1:3]), as.matrix(ml2[,4:5]), n_iter = 1000, proposal = c("normbeta"))
str(sampler_result_fix, max.level = 1)
beta_samples_fix <- matrix(0, nrow = length(sampler_result_fix$beta_chain), ncol = length(B2))
for(i in 1:length(sampler_result_fix$beta_chain)){
  beta_samples_fix[i,] <- c(sampler_result_fix$beta_chain[[i]])
}
plot(beta_samples_fix[,2])
abline(h = B2[2,1])
head(beta_samples_fix)
B2

w_preds <- lapply(seq_along(sampler_result_fix$y_chain), function(i) {
  sample_posterior_predictive(X = as.matrix(ml2[,4:5]),
                              beta = sampler_result_fix$beta_chain[[i]],
                              Sigma = sampler_result_fix$sigma_chain[[i]],
                              PA = 200,
                              mixed = FALSE
  )
})

mahal_resids <- compute_mahalanobis_resids(as.matrix(ml2[,1:3]), w_preds)
summary(mahal_resids)

beta_chain_array <- simplify2array(sampler_result_fix$beta_chain)
plot_trace_and_summary(beta_chain_array, "beta")
B2
sigma_chain_array <- simplify2array(sampler_result_fix$sigma_chain)
plot_trace_and_summary(sigma_chain_array, "sigma")
diag(J-1)


# with other sim strategy
df_sim2 <- function(m, ns, PA = NULL) {
  set.seed(42)
  m <- m; n_i <- rep(ns, m); N <- sum(n_i)
  X <- cbind(1, rnorm(N))
  id <- rep(seq_len(m), each=ns);
  if(is.null(PA)){
    PA <- rep(200, N)
  }
  d <- 2
  beta_true <- matrix(c(-2, 0.5, 1, -0.3), nrow=2, ncol=2)  # p x d
  Sigma_true <- diag(c(1.2, .7))
  Y_true <- X %*% beta_true + rmvn(N, rep(0, d), Sigma_true)
  P_true <- alr_inv(Y_true)
  W <- t(apply(P_true, 1, function(p) rmultinom(1, PA[1], p)))
  list(W=W, X=X, id=id, PA=PA, beta_true=beta_true, Sigma_true=Sigma_true)
}

sim2 <- df_sim2(m = 20, ns = 6)
# Test the fixed Gibbs sampler wrapper
sampler_result_fix2 <- run_fixed_gibbs_sampler(sim2$W, sim2$X, n_iter = 1000, burn_in = 180, thin = 2, proposal = "normbeta")
str(sampler_result_fix2, max.level = 1)
beta_samples2 <- matrix(0, nrow = length(sampler_result_fix2$beta_chain), ncol = length(sim2$beta_true))
for(i in 1:length(sampler_result_fix2$beta_chain)){
  beta_samples2[i,] <- c(sampler_result_fix2$beta_chain[[i]])
}
plot(beta_samples2[,1])
abline(h = sim2$beta_true[1,1])
head(beta_samples2)
sim2$beta_true

w_preds2 <- lapply(seq_along(sampler_result_fix2$y_chain), function(i) {
  sample_posterior_predictive(X = sim2$X,
                              beta = sampler_result_fix2$beta_chain[[i]],
                              Sigma = sampler_result_fix2$sigma_chain[[i]],
                              PA = sim2$PA,
                              mixed = FALSE
  )
})

mahal_resids <- compute_mahalanobis_resids(sim2$W, w_preds2)
summary(mahal_resids)
plot(mahal_resids)

beta_chain_array <- simplify2array(sampler_result_fix2$beta_chain)
plot_trace_and_summary(beta_chain_array, "beta")
sim2$beta_true
sigma_chain_array <- simplify2array(sampler_result_fix2$sigma_chain)
plot_trace_and_summary(sigma_chain_array, "sigma")
sim2$Sigma_true


#Under More distributions, incl.
# Dirichlet
sampler_result_fix3 <- run_fixed_gibbs_sampler(as.matrix(ml3[,1:3]), as.matrix(ml3[,4:5]), n_iter = 1000, proposal = c("normbeta"))
str(sampler_result_fix3, max.level = 1)
beta_samples_fix3 <- matrix(0, nrow = length(sampler_result_fix3$beta_chain), ncol = length(B2))
for(i in 1:length(sampler_result_fix3$beta_chain)){
  beta_samples_fix3[i,] <- c(sampler_result_fix3$beta_chain[[i]])
}

w_preds3 <- lapply(seq_along(sampler_result_fix3$y_chain), function(i) {
  sample_posterior_predictive(X = as.matrix(ml3[,4:5]),
                              beta = sampler_result_fix3$beta_chain[[i]],
                              Sigma = sampler_result_fix3$sigma_chain[[i]],
                              PA = 200,
                              mixed = FALSE
  )
})

mahal_resids <- compute_mahalanobis_resids(as.matrix(ml3[,1:3]), w_preds3)
summary(mahal_resids)

beta_chain_array <- simplify2array(sampler_result_fix3$beta_chain)
plot_trace_and_summary(beta_chain_array, "beta")
sigma_chain_array <- simplify2array(sampler_result_fix3$sigma_chain)
plot_trace_and_summary(sigma_chain_array, "sigma")

# And with positive correlations
sampler_result_fix4 <- run_fixed_gibbs_sampler(as.matrix(ml4[,1:3]), as.matrix(ml4[,4:5]), n_iter = 1000, proposal = c("normbeta"))
str(sampler_result_fix, max.level = 1)
beta_samples_fix4 <- matrix(0, nrow = length(sampler_result_fix4$beta_chain), ncol = length(B2))
for(i in 1:length(sampler_result_fix4$beta_chain)){
  beta_samples_fix4[i,] <- c(sampler_result_fix4$beta_chain[[i]])
}
plot(beta_samples_fix4[,1])
abline(h = B2[1,1])
head(beta_samples_fix4)
B2

w_preds4 <- lapply(seq_along(sampler_result_fix4$y_chain), function(i) {
  sample_posterior_predictive(X = as.matrix(ml4[,4:5]),
                              beta = sampler_result_fix4$beta_chain[[i]],
                              Sigma = sampler_result_fix4$sigma_chain[[i]],
                              PA = 200,
                              mixed = FALSE
  )
})

mahal_resids <- compute_mahalanobis_resids(as.matrix(ml4[,1:3]), w_preds4)
print(summary(mahal_resids))
plot(mahal_resids)
qqnorm(mahal_resids)
abline(0, 1, col = "red", lwd = 2)
beta_chain_array <- simplify2array(sampler_result_fix4$beta_chain)
plot_trace_and_summary(beta_chain_array, "beta")
B2
sigma_chain_array <- simplify2array(sampler_result_fix4$sigma_chain)
plot_trace_and_summary(sigma_chain_array, "sigma")
matrix(c(1,-.9,-.9,4), ncol = 2)

