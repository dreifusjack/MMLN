#' Real‐Data Example: Pollen Data – MN, DM & MLN Models
#'
#' This function reproduces the pollen‐data example from Gerber & Craig, fitting three
#' models to the “pollen” counts (multinomial logit, Dirichlet‐multinomial, and
#' multinomial-logistic‐normal), drawing replicate from fitted sampling distributions, and computing
#' Mahalanobis residuals for each model.
#'
#' @param n_iter   Integer; total number of MCMC iterations for the MLN model (default 1000)
#' @param burn_in  Integer; number of initial MLN iterations to discard (default 400)
#' @param thin     Integer; MLN thinning interval (default 2)
#' @param proposal Character; one of `"norm"`, `"beta"`, or `"normbeta"` for the MLN sampler (default `"normbeta"`)
#' @param P        Integer; number of fitted sampling distribution replicates per model (default 1000)
#'
#' @return A list with components:
#' \describe{
#'   \item{fit_mlr}{The `MGLMreg` object for the multinomial logit fit.}
#'   \item{fit_dm}{The `MGLMreg` object for the Dirichlet‐multinomial fit.}
#'   \item{fit_mln}{The list returned by `FMLN()` for the fixed‐effects MLN fit.}
#'   \item{Y_pred_mlr}{List of length P of \(N\times J\) count matrices sampled from the MN model.}
#'   \item{Y_pred_dm}{List of length P of \(N\times J\) count matrices sampled from the DM model.}
#'   \item{Y_pred_mln}{List of length P of \(N\times J\) count matrices sampled from the MLN model (using posterior‐mean parameters).}
#'   \item{resids_mlr}{Mahalanobis residuals (`mdres`) for the multinomial logit model.}
#'   \item{resids_dm}{Mahalanobis residuals (`mdres`) for the Dirichlet‐multinomial model.}
#'   \item{resids_mln}{Mahalanobis residuals (`mdres`) for the MLN model.}
#' }
#'
#' @examples
#' \dontrun{
#' # run all three fits & diagnostics
#' pollen_res <- run_pollen_models(n_iter = 1500, burn_in = 500, thin = 5, P = 500)
#'
#' # view KS‐tests & QQ‐plots
#' summary(pollen_res$resids_mlr)
#' summary(pollen_res$resids_dm)
#' summary(pollen_res$resids_mln)
#' }
#'
#' @importFrom stats rmultinom
#' @importFrom mvnfast rmvn
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#' @importFrom MGLM MGLMreg rdirmn
#' @export
run_pollen_models <- function(n_iter   = 1000,
                              burn_in  = 400,
                              thin     = 2,
                              proposal = "normbeta",
                              P        = 1000) {
  ## dependencies
  if (!requireNamespace("MM",   quietly = TRUE)) stop("Install the 'MM' package to load pollen data")
  if (!requireNamespace("MGLM", quietly = TRUE)) stop("Install the 'MGLM' package for MGLMreg()")

  ## 1) load data
  data(pollen, package = "MM")
  Y  <- as.matrix(pollen)
  N  <- nrow(Y)
  PA <- rowSums(Y)
  J  <- ncol(Y)

  X <- matrix(1, nrow = N, ncol = 1)  # intercept‐only design

  ## 2) fit M-Logit (multinomial) and Dirichlet-multinomial via MGLM
  message("▶ Fitting multinomial-logistic model using MGLMreg()")
  fit_mlr <- suppressWarnings(MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1,
                     data = as.data.frame(pollen), dist = "MN"))
  cat("Done: \n")
  message("▶ Fitting Dirichlet-lmultinomial model using MGLMreg()")
  fit_dm  <- suppressWarnings(MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1,
                     data = as.data.frame(pollen), dist = "DM"))

  ## 3) fit fixed‐effects MLN
  cat("Done: \n")
  message("▶ Fitting multinomial-logistic-normal model using FMLN()")
  fit_mln <- FMLN(
    Y              = Y,
    X              = X,
    n_iter         = n_iter,
    burn_in        = burn_in,
    thin           = thin,
    proposal       = proposal,
    verbose        = TRUE
  )

  ## 4) prepare predictive replicates
  cat("Done: \n")
  message("▶ Drawing sampling distributions of predicted values from each model")
  # 4a) MN model: draw P full-dataset replicates
  probs_mlr <- fit_mlr@fitted
  Y_pred_mlr <- vector("list", P)
  for (p_i in seq_len(P)) {
    M <- t(sapply(seq_len(N),
                  function(i) rmultinom(1, size = PA[i], prob = probs_mlr[i, ])))
    Y_pred_mlr[[p_i]] <- M
  }

  # 4b) DM model: draw P replicates with Dirichlet‐multinomial
  alpha_hat <- exp(fit_dm@coefficients)
  Y_pred_dm <- vector("list", P)
  for (p_i in seq_len(P)) {
    M <- t(sapply(seq_len(N),
                  function(i) MGLM::rdirmn(n    = 1,
                                           size = PA[i],
                                           alpha = alpha_hat)))
    Y_pred_dm[[p_i]] <- M
  }

  # 4c) MLN model: use posterior‐mean β and Σ to draw P replicates
  # posterior means:
  beta_arr  <- simplify2array(fit_mln$beta_chain)    # p × d × n_saves
  Sigma_arr <- simplify2array(fit_mln$sigma_chain)   # d × d × n_saves
  beta_mean  <- apply(beta_arr,  c(1,2), mean)
  Sigma_mean <- apply(Sigma_arr, c(1,2), mean)

  Y_pred_mln <- vector("list", P)
  for (p_i in seq_len(P)) {
    Y_pred_mln[[p_i]] <- sample_posterior_predictive(
      X     = X,
      beta  = beta_mean,
      Sigma = Sigma_mean,
      PA    = PA,
      mixed = FALSE
    )
  }

  ## 5) compute Mahalanobis‐residuals
  cat("Done: \n")  # ensure we start on a fresh line
  message("▶ Computing MDRes for multinomial‐logit model")
  resids_mlr <- MDres(Y, Y_pred_mlr)

  cat("Done: \n")  # ensure we start on a fresh line
  message("▶ Computing MDRes for Dirichlet-multinomial model")
  resids_dm  <- MDres(Y, Y_pred_dm)

  cat("Done: \n")  # ensure we start on a fresh line
  message("▶ Computing MDRes for multinomial‐logistic-normal model")
  resids_mln <- MDres(Y, Y_pred_mln)

  ## return all results
  output = list(
    fit_mlr    = fit_mlr,
    fit_dm     = fit_dm,
    fit_mln    = fit_mln,
    Y_pred_mlr = Y_pred_mlr,
    Y_pred_dm  = Y_pred_dm,
    Y_pred_mln = Y_pred_mln,
    resids_mlr = resids_mlr,
    resids_dm  = resids_dm,
    resids_mln = resids_mln
  )

  return(output)
}
