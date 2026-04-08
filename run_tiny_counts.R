library(MMLN)

# 1. Simulate tiny_counts dataset
set.seed(42003)

sim <- simulate_mixed_mln_data(
  m      = 40,
  n_i    = 30,
  p      = 5,
  d      = 3,
  beta   = matrix(c(
    -0.5,  0.3,  0.1,  0.2, -0.1,
     0.4, -0.2,  0.3, -0.1,  0.2,
     0.1,  0.5, -0.2,  0.0,  0.3
  ), nrow = 5, ncol = 3),
  Sigma  = diag(3),
  Phi    = 0.3 * diag(3),
  n_mean = 5
)

# 2. Fit a fixed-effects MLN model
res_f <- FMLN(
  Y        = sim$Y,
  X        = sim$X,
  n_iter   = 1000,
  burn_in  = 300,
  thin     = 2,
  proposal = "normbeta",
  verbose  = TRUE
)

# 3. Fit a mixed-effects MLN model
res_m <- MMLN(
  Y        = sim$Y,
  X        = sim$X,
  Z        = sim$Z,
  n_iter   = 1000,
  burn_in  = 300,
  thin     = 2,
  proposal = "normbeta",
  verbose  = TRUE
)

# 4. Trace plots & posterior summaries
beta_chain_array <- simplify2array(res_m$beta_chain)
trace_stats      <- plot_trace_and_summary(beta_chain_array, "beta")
trace_stats
sim$beta
par(mfrow = c(1, 1))

# 5. Compute model DICs
ll_chain <- sapply(res_m$w_chain,
                   function(W) dmnl_loglik(W, sim$Y))
W_hat   <- alr(compress_counts(sim$Y) / rowSums(sim$Y))
ll_hat  <- dmnl_loglik(W_hat, sim$Y)
dic_res <- compute_dic(ll_chain, ll_hat)

# 6. Posterior predictive simulation and Mahalanobis residuals
Y_pred_list <- lapply(seq_along(res_m$w_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_m$beta_chain[[i]],
                              Sigma = res_m$sigma_chain[[i]],
                              n = sim$n,
                              Z = sim$Z,
                              psi = res_m$psi_chain[[i]],
                              mixed = TRUE,
                              verbose = FALSE
  )
})
resids <- MDres(sim$Y, Y_pred_list)
summary(resids)

# 7. Compare to incorrect model fit (fixed model fit to mixed data)
Y_pred_list_ovd <- lapply(seq_along(res_f$w_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_f$beta_chain[[i]],
                              Sigma = res_f$sigma_chain[[i]],
                              n = sim$n,
                              mixed = FALSE,
                              verbose = FALSE
  )
})
resids_ovd <- MDres(sim$Y, Y_pred_list_ovd)
summary(resids_ovd)

# 8. Incorrect model fit should have higher DIC
ll_chain_ovd <- sapply(res_f$w_chain,
                       function(W) dmnl_loglik(W, sim$Y))
dic_res_ovd <- compute_dic(ll_chain_ovd, ll_hat)
dic_res_ovd$DIC > dic_res$DIC
