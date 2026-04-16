#!/usr/bin/env Rscript
# run_real_data_reference.R -- Run upstream FMLN/MMLN/MDres on real pollen and MLB datasets.
#
# Runs:
#   - FMLN on pollen data        (intercept-only, from MM package)
#   - FMLN on MLB data           (fixed effects only, no random intercepts)
#   - MMLN on MLB data           (player-level random effects)
#   - MDres on MLB               (using upstream MMLN chain)
#
# Saves timing + results to fixtures/ for later comparison with the fork.
#
# Usage: Rscript run_real_data_reference.R <repo_root> [<n_iter>] [<burn_in>] [<mdres_p>]

cat("=== Real-data reference run (upstream) ===\n\n")

library(MMLN)

args      <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
n_iter    <- if (length(args) >= 2) as.integer(args[2]) else 1000L
burn_in   <- if (length(args) >= 3) as.integer(args[3]) else 300L
mdres_p   <- if (length(args) >= 4) as.integer(args[4]) else 50L
thin      <- 2L

fix_dir <- file.path(repo_root, "optimization-testing/fixtures")
dir.create(fix_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  n_iter=%d  burn_in=%d  thin=%d  mdres_p=%d\n\n",
            n_iter, burn_in, thin, mdres_p))

# ---------------------------------------------------------------------------
# Helper: save a fitted model with timing metadata
# ---------------------------------------------------------------------------
save_fit <- function(fit, elapsed, file, N, J, m = NA_integer_) {
  saveRDS(list(fit = fit, elapsed = elapsed,
               n_iter = n_iter, burn_in = burn_in,
               N = N, J = J, m = m),
          file = file)
}

# ---------------------------------------------------------------------------
# 1. Pollen -- FMLN
# ---------------------------------------------------------------------------
pollen_file <- file.path(fix_dir, "real_pollen_fmln.rds")

if (file.exists(pollen_file)) {
  cat("  Pollen FMLN: SKIP (fixture exists)\n")
} else {
  if (!requireNamespace("MM", quietly = TRUE))
    stop("Install the 'MM' package first: install.packages('MM')\n  (requires: sudo apt-get install libgmp-dev)")

  data(pollen, package = "MM", envir = environment())
  Y_pol <- as.matrix(as.data.frame(pollen))
  X_pol <- matrix(1, nrow = nrow(Y_pol), ncol = 1)

  cat(sprintf("  Pollen FMLN | N=%d  J=%d  n_iter=%d\n", nrow(Y_pol), ncol(Y_pol), n_iter))
  flush.console()

  set.seed(99001)
  t0 <- Sys.time()
  fit <- FMLN(Y = Y_pol, X = X_pol, n_iter = n_iter, burn_in = burn_in,
              thin = thin, proposal = "normbeta", verbose = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  save_fit(fit, elapsed, pollen_file, nrow(Y_pol), ncol(Y_pol))
  cat(sprintf("\n  Pollen FMLN done: %.1fs\n\n", elapsed))
}

# ---------------------------------------------------------------------------
# 2 & 3. MLB -- load data once, run FMLN then MMLN
# ---------------------------------------------------------------------------
mlb_fmln_file <- file.path(fix_dir, "real_mlb_fmln.rds")
mlb_mmln_file <- file.path(fix_dir, "real_mlb_mmln.rds")
mlb_mdres_file <- file.path(fix_dir, "real_mlb_mdres.rds")

need_mlb <- !file.exists(mlb_fmln_file) || !file.exists(mlb_mmln_file) || !file.exists(mlb_mdres_file)

if (!need_mlb) {
  cat("  MLB FMLN:  SKIP (fixture exists)\n")
  cat("  MLB MMLN:  SKIP (fixture exists)\n")
  cat("  MLB MDres: SKIP (fixture exists)\n")
} else {
  if (!requireNamespace("Lahman", quietly = TRUE))
    stop("Install the 'Lahman' package first: install.packages('Lahman')")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Install the 'dplyr' package first: install.packages('dplyr')")

  cat("  Loading MLB data via clean_Lahman_data()...\n")
  flush.console()
  mlb   <- clean_Lahman_data()
  Y_mlb <- mlb$Y
  X_mlb <- mlb$X
  Z_mlb <- mlb$Z
  m_mlb <- ncol(Z_mlb)
  n_counts <- rowSums(Y_mlb)

  # --- MLB FMLN ---
  if (file.exists(mlb_fmln_file)) {
    cat("  MLB  FMLN: SKIP (fixture exists)\n")
  } else {
    cat(sprintf("  MLB  FMLN | N=%d  J=%d  p=%d  n_iter=%d\n",
                nrow(Y_mlb), ncol(Y_mlb), ncol(X_mlb), n_iter))
    flush.console()

    set.seed(99003)
    t0 <- Sys.time()
    fit <- FMLN(Y = Y_mlb, X = X_mlb, n_iter = n_iter, burn_in = burn_in,
                thin = thin, proposal = "normbeta", verbose = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    save_fit(fit, elapsed, mlb_fmln_file, nrow(Y_mlb), ncol(Y_mlb))
    cat(sprintf("\n  MLB  FMLN done: %.1fs\n\n", elapsed))
  }

  # --- MLB MMLN ---
  if (file.exists(mlb_mmln_file)) {
    cat("  MLB  MMLN: SKIP (fixture exists)\n")
  } else {
    cat(sprintf("  MLB  MMLN | N=%d  J=%d  p=%d  m=%d  n_iter=%d\n",
                nrow(Y_mlb), ncol(Y_mlb), ncol(X_mlb), m_mlb, n_iter))
    flush.console()

    set.seed(99002)
    t0 <- Sys.time()
    fit <- MMLN(Y = Y_mlb, X = X_mlb, Z = Z_mlb, n_iter = n_iter, burn_in = burn_in,
                thin = thin, mh_scale = 1, proposal = "normbeta", verbose = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    save_fit(fit, elapsed, mlb_mmln_file, nrow(Y_mlb), ncol(Y_mlb), m_mlb)
    cat(sprintf("\n  MLB  MMLN done: %.1fs\n\n", elapsed))
  }

  # --- MLB MDres (using the upstream MMLN chain just saved) ---
  if (file.exists(mlb_mdres_file)) {
    cat("  MLB  MDres: SKIP (fixture exists)\n")
  } else {
    mmln_fix <- readRDS(mlb_mmln_file)
    n_saved  <- length(mmln_fix$fit$beta_chain)
    beta_hat  <- mmln_fix$fit$beta_chain[[n_saved]]
    sigma_hat <- mmln_fix$fit$sigma_chain[[n_saved]]
    psi_hat   <- mmln_fix$fit$psi_chain[[n_saved]]

    cat(sprintf("  MLB  MDres | P=%d\n", mdres_p))
    flush.console()

    set.seed(99004)
    t0 <- Sys.time()
    Y_pred_list <- vector("list", mdres_p)
    for (j in seq_len(mdres_p)) {
      Y_pred_list[[j]] <- sample_posterior_predictive(
        X = X_mlb, beta = beta_hat, Sigma = sigma_hat,
        n = n_counts, Z = Z_mlb, psi = psi_hat,
        mixed = TRUE, verbose = FALSE
      )
    }
    mdres_result <- MDres(Y_mlb, Y_pred_list)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    saveRDS(list(mdres = mdres_result, Y_pred_list = Y_pred_list,
                 elapsed = elapsed, mdres_p = mdres_p),
            file = mlb_mdres_file)
    cat(sprintf("  MLB  MDres done: %.1fs\n\n", elapsed))
  }
}

cat("Done. Reference fixtures saved to fixtures/\n")
