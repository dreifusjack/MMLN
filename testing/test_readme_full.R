library(MMLN)

# End-to-end smoke test mirroring README "Quick Start" workflow.
# Goal: ensure the full pipeline runs and that outputs look sane (no outliers),
# while keeping runtime reasonable for local debugging.
#
# To run:
#   devtools::load_all()
#   source("testing/test_readme_full.R")

run_readme_smoke <- function(
  label = NULL,
  seed = 42,
  m = 10,
  n_i = 10,
  p = 3,
  d = 2,
  beta = matrix(c(0.5, -1, 0.2, 0.3, 0.7, -0.4), 3, 2),
  n_iter = 250,
  burn_in = 50,
  thin = 2,
  proposal = "normbeta",
  verbose = FALSE,
  P_pred = 30
) {
  if (!is.null(seed)) set.seed(seed)

  sim <- simulate_mixed_mln_data(
    m = m,
    n_i = n_i,
    p = p,
    d = d,
    beta = beta,
    Sigma = diag(2),
    Phi = 5 * diag(2),
    n_mean = 200
  )

  cat("\n== Simulated truth beta ==\n")
  print(sim$beta)

  cat("\n== Fit fixed-effects FMLN ==\n")
  res_f <- FMLN(
    Y = sim$Y,
    X = sim$X,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    proposal = proposal,
    verbose = verbose
  )

  cat("\n== Fit mixed-effects MMLN ==\n")
  res_m <- MMLN(
    Y = sim$Y,
    X = sim$X,
    Z = sim$Z,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    proposal = proposal,
    verbose = verbose
  )

  # Posterior beta mean (p x d)
  beta_arr <- simplify2array(res_m$beta_chain)
  beta_mean_m <- apply(beta_arr, c(1, 2), mean)

  beta_arr_f <- simplify2array(res_f$beta_chain)
  beta_mean_f <- apply(beta_arr_f, c(1, 2), mean)

  cat("\n== Posterior beta means (MMLN vs FMLN) ==\n")
  print(list(beta_mean_m = beta_mean_m, beta_mean_f = beta_mean_f))

  # Compare to simulated truth (sanity; not a strict equivalence criterion).
  beta_diff_m <- beta_mean_m - sim$beta
  beta_diff_f <- beta_mean_f - sim$beta
  max_abs_diff_m <- max(abs(beta_diff_m))
  max_abs_diff_f <- max(abs(beta_diff_f))

  cat("\n== Max abs beta error vs simulated truth ==\n")
  cat(sprintf("  MMLN max|beta_mean - sim$beta| = %.4f\n", max_abs_diff_m))
  cat(sprintf("  FMLN max|beta_mean - sim$beta| = %.4f\n", max_abs_diff_f))

  # DIC computations follow README logic
  cat("\n== DIC computations ==\n")
  ll_chain_m <- sapply(res_m$w_chain, function(W) dmnl_loglik(W, sim$Y))
  W_hat <- alr(compress_counts(sim$Y) / rowSums(sim$Y))
  ll_hat <- dmnl_loglik(W_hat, sim$Y)
  dic_res_m <- compute_dic(ll_chain_m, ll_hat)

  ll_chain_f <- sapply(res_f$w_chain, function(W) dmnl_loglik(W, sim$Y))
  dic_res_f <- compute_dic(ll_chain_f, ll_hat)

  cat("\nDIC (lower is better):\n")
  cat(sprintf("  MMLN DIC = %.3f\n", dic_res_m$DIC))
  cat(sprintf("  FMLN DIC = %.3f\n", dic_res_f$DIC))

  # Posterior predictive checks + MDres (subsample posterior draws for speed)
  cat("\n== Posterior predictive + MDres (subsample) ==\n")
  n_save_m <- length(res_m$w_chain)
  n_save_f <- length(res_f$w_chain)
  P_m <- min(P_pred, n_save_m)
  P_f <- min(P_pred, n_save_f)

  idx_m <- seq_len(P_m)
  idx_f <- seq_len(P_f)

  Y_pred_m <- lapply(idx_m, function(i) {
    sample_posterior_predictive(
      X = sim$X,
      beta = res_m$beta_chain[[i]],
      Sigma = res_m$sigma_chain[[i]],
      n = sim$n,
      Z = sim$Z,
      psi = res_m$psi_chain[[i]],
      mixed = TRUE,
      verbose = FALSE
    )
  })

  resids_m <- MDres(sim$Y, Y_pred_m)
  ks_m <- suppressWarnings(ks.test(resids_m, "pnorm", mean = 0, sd = 1, exact = FALSE))

  Y_pred_f <- lapply(idx_f, function(i) {
    sample_posterior_predictive(
      X = sim$X,
      beta = res_f$beta_chain[[i]],
      Sigma = res_f$sigma_chain[[i]],
      n = sim$n,
      mixed = FALSE,
      verbose = FALSE
    )
  })

  resids_f <- MDres(sim$Y, Y_pred_f)
  ks_f <- suppressWarnings(ks.test(resids_f, "pnorm", mean = 0, sd = 1, exact = FALSE))

  cat("\nKS test for normality of Mahalanobis residuals:\n")
  cat(sprintf("  MMLN: D=%.4f, p-value=%.4f\n", ks_m$statistic, ks_m$p.value))
  cat(sprintf("  FMLN: D=%.4f, p-value=%.4f\n", ks_f$statistic, ks_f$p.value))

  # Outlier sanity checks (deterministic bounds; not a correctness proof)
  max_abs_res_m <- max(abs(resids_m), na.rm = TRUE)
  max_abs_res_f <- max(abs(resids_f), na.rm = TRUE)

  cat("\nResidual sanity:\n")
  cat(sprintf("  max|resids| MMLN = %.3f\n", max_abs_res_m))
  cat(sprintf("  max|resids| FMLN = %.3f\n", max_abs_res_f))

  if (!is.null(label)) {
    cat("\n== Compact comparison row ==\n")
    comp_row <- data.frame(
      label = label,
      seed = seed,
      n_iter = n_iter,
      burn_in = burn_in,
      thin = thin,
      proposal = proposal,
      P_pred = P_pred,
      max_abs_beta_err_mmln = max_abs_diff_m,
      max_abs_beta_err_fmln = max_abs_diff_f,
      dic_mmln = dic_res_m$DIC,
      dic_fmln = dic_res_f$DIC,
      ks_p_mmln = ks_m$p.value,
      ks_p_fmln = ks_f$p.value,
      max_abs_resids_mmln = max_abs_res_m,
      max_abs_resids_fmln = max_abs_res_f
    )
    print(comp_row, row.names = FALSE)
  }

  # Return everything for programmatic inspection
  list(
    sim = sim,
    res_f = res_f,
    res_m = res_m,
    beta_mean_m = beta_mean_m,
    beta_mean_f = beta_mean_f,
    dic_res_m = dic_res_m,
    dic_res_f = dic_res_f,
    ks_m = ks_m,
    ks_f = ks_f,
    resids_m = resids_m,
    resids_f = resids_f
  )
}

if (isTRUE(getOption("MMLN_RUN_README_SMOKE", FALSE))) {
  out <- run_readme_smoke(verbose = FALSE)

  cat("\n== Summary ==\n")
  cat(sprintf("MMLN DIC: %.3f\n", out$dic_res_m$DIC))
  cat(sprintf("FMLN DIC: %.3f\n", out$dic_res_f$DIC))
  cat(sprintf("MMLN KS p-value: %.4f\n", out$ks_m$p.value))
  cat(sprintf("FMLN KS p-value: %.4f\n", out$ks_f$p.value))

  cat("\n(If max abs beta error >> 1, the chain likely isn't mixing well or a correctness bug is still present.)\n")
}

