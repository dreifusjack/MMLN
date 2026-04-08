library(parallel)

# --- Posterior Diagnostics and DIC ---

plot_trace_and_summary <- function(chain_array, param_name = "param", max_rows = 4) {
  dims <- dim(chain_array)
  par(mfrow = c(2, 2))
  n_rows <- min(dims[1], max_rows)
  for (i in 1:n_rows) {
    for (j in 1:dims[2]) {
      trace <- chain_array[i, j, ]
      plot(trace, type = "l", main = paste0("Trace: ", param_name, "[", i, ",", j, "]"), ylab = "Value", xlab = "Iteration")
      abline(h = mean(trace), col = "red")
    }
  }

  summary_array <- array(NA, dim = c(dim(chain_array)[1], dim(chain_array)[2], 4))
  for (i in seq_len(dim(chain_array)[1])) {
    for (j in seq_len(dim(chain_array)[2])) {
      trace <- chain_array[i, j, ]
      summary_array[i, j, ] <- c(
        mean = mean(trace),
        sd = sd(trace),
        `2.5%` = quantile(trace, 0.025),
        `97.5%` = quantile(trace, 0.975)
      )
    }
  }
  dimnames(summary_array) <- list(
    paste0("row", seq_len(dim(chain_array)[1])),
    paste0("col", seq_len(dim(chain_array)[2])),
    c("mean", "sd", "2.5%", "97.5%")
  )

  par(mfrow = c(1, 1))
  return(aperm(summary_array, c(2, 3, 1)))
}

compute_dic <- function(loglik_chain, loglik_hat) {
  pD <- 2 * (loglik_hat - mean(loglik_chain))
  DIC <- -2 * loglik_hat + 2 * pD
  list(DIC = DIC, pD = pD)
}

sample_posterior_predictive <- function(X, beta, Sigma, PA,
                                        Z = NULL, psi = NULL,
                                        mixed = TRUE) {
  N <- nrow(X)
  d <- ncol(Sigma)
  # linear predictor
  mu <- X %*% beta
  if(mixed) {
    if(is.null(Z) || is.null(psi)) stop("For mixed = TRUE, must supply Z and psi")
    mu <- mu + Z %*% psi
  }
  # simulate latent alr-scale outcomes
  Y_hat <- mu + mvnfast::rmvn(N, mu = rep(0, d), sigma = Sigma)
  # transform back to simplex
  P_hat <- alr_inv(Y_hat)
  # sample counts per observation
  W <- matrix(NA, nrow = N, ncol = d + 1)
  for(i in seq_len(N)) {
    size_i <- if(length(PA) > 1) PA[i] else PA
    W[i, ] <- rmultinom(1, size = size_i, prob = P_hat[i, ])
  }
  W
}


compute_mahalanobis_resids <- function(W, W_pred_list) {

  W_obs <- compress_counts(W)
  alr_obs <- alr(W_obs)

  N <- nrow(W_obs)
  d <- ncol(W_obs) - 1
  P <- length(W_pred_list)

  pred_array <- array(NA, dim = c(N, d, P))
  for (j in seq_len(P)) {
    pred_array[,,j] <- alr(compress_counts(W_pred_list[[j]]))
  }

  mu_all <- apply(pred_array, 1:2, mean)
  Sigma_all <- lapply(seq_len(N), function(i) cov(t(pred_array[i,,])))

  mds_list <- vector("list", N)
  if (interactive()) cat("Computing Mahalanobis distances:\n")
  start_time <- Sys.time()
  for (i in seq_len(N)) {
    if (interactive() && (i %% max(1, floor(N / 100)) == 0 || i == N)) {
      pct <- floor(100 * i / N)
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / i) * (N - i)
      cat(sprintf("\r[%3d%%] ETA: %s", pct, format(.POSIXct(eta, tz="GMT"), "%M:%S")))
      flush.console()
    }
    pred_i <- t(pred_array[i,,])
    y_obs_i <- alr_obs[i, ]
    obsi <- rbind(y_obs_i, pred_i)
    mds_list[[i]] <- apply(obsi, 1, mahalanobis, center = mu_all[i,], cov = Sigma_all[[i]])
  }

  u_resids <- vapply(seq_len(N), function(i) {
    obs_val <- mds_list[[i]][1]
    post_vals <- mds_list[[i]][-1]
    ecdf_i <- ecdf(post_vals)

    if (obs_val <= min(post_vals)) {
      minpct <- 0
      maxpct <- ecdf_i(min(post_vals))
      if (maxpct == 0) maxpct <- (length(post_vals) - 1) / length(post_vals)
    } else if (obs_val > max(post_vals)) {
      minpct <- ecdf_i(max(post_vals))
      maxpct <- 1
      if (minpct == 1) minpct <- (length(post_vals) - 1) / length(post_vals)
    } else {
      sorted_vals <- sort(post_vals)
      lower <- max(which(sorted_vals < obs_val))
      minpct <- ecdf_i(sorted_vals[lower])
      maxpct <- ecdf_i(obs_val)
      if (minpct == 1) minpct <- (length(post_vals) - 1) / length(post_vals)
      if (maxpct == 0) maxpct <- (length(post_vals) - 1) / length(post_vals)
    }

    runif(1, minpct, maxpct)
  }, numeric(1))

  z_resids <- qnorm(u_resids)

  class(z_resids) <- "mahalanobis_residuals"
  return(z_resids)

}

summary.mahalanobis_residuals <- function(object, ...) {
  # `object` is just a numeric vector of residuals
  cat("Kolmogorov–Smirnov test for normality of Mahalanobis residuals:\n")
  ks <- suppressWarnings(
    ks.test(object,
            "pnorm",
            mean = 0, sd = 1,
            exact = FALSE)   # use asymptotic p‐value to avoid ties warning
  )
  print(ks)

  # QQ‐plot
  qqnorm(object, main = "QQ Plot of Mahalanobis Residuals")
  abline(0, 1, col='red', lwd = 2)

  invisible(ks)
}
