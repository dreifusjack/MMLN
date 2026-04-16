#!/usr/bin/env Rscript
# run_real_data_reference.R -- Run upstream MMLN/FMLN on real pollen and MLB datasets.
#
# Installs upstream first (via Makefile), then runs:
#   - FMLN on pollen data (intercept-only, from MM package)
#   - MMLN on MLB data   (full player-level random effects, from Lahman package)
#
# Saves timing + chain results to fixtures/ for later comparison with the fork.
#
# Usage: Rscript run_real_data_reference.R <repo_root> [<n_iter>] [<burn_in>]
#   n_iter:   MCMC iterations (default: 2000)
#   burn_in:  burn-in to discard (default: 500)

cat("=== Real-data reference run (upstream) ===\n\n")

library(MMLN)

args     <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
n_iter   <- if (length(args) >= 2) as.integer(args[2]) else 2000L
burn_in  <- if (length(args) >= 3) as.integer(args[3]) else 500L
thin     <- 2L

fix_dir <- file.path(repo_root, "optimization-testing/fixtures")
dir.create(fix_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  n_iter=%d  burn_in=%d  thin=%d\n\n", n_iter, burn_in, thin))

# ---------------------------------------------------------------------------
# Helper: format seconds nicely
# ---------------------------------------------------------------------------
fmt_time <- function(secs) {
  secs <- round(secs)
  if (secs >= 3600) {
    sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
  } else {
    sprintf("%02d:%02d", secs %/% 60, secs %% 60)
  }
}

# ---------------------------------------------------------------------------
# 1. Pollen -- FMLN (no random effects)
# ---------------------------------------------------------------------------
pollen_file <- file.path(fix_dir, "real_pollen_fmln.rds")

if (file.exists(pollen_file)) {
  cat("  Pollen FMLN: SKIP (fixture exists)\n")
} else {
  if (!requireNamespace("MM", quietly = TRUE))
    stop("Install the 'MM' package first: install.packages('MM', repos='https://cloud.r-project.org')")

  data(pollen, package = "MM", envir = environment())
  pollen_df <- as.data.frame(pollen)
  Y_pol <- as.matrix(pollen_df)
  X_pol <- matrix(1, nrow = nrow(Y_pol), ncol = 1)

  cat(sprintf("  Pollen FMLN | N=%d  J=%d  n_iter=%d\n", nrow(Y_pol), ncol(Y_pol), n_iter))
  flush.console()

  set.seed(99001)
  start <- Sys.time()
  fit_pollen <- FMLN(
    Y        = Y_pol,
    X        = X_pol,
    n_iter   = n_iter,
    burn_in  = burn_in,
    thin     = thin,
    proposal = "normbeta",
    verbose  = TRUE
  )
  elapsed_pol <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  saveRDS(
    list(
      fit     = fit_pollen,
      elapsed = elapsed_pol,
      n_iter  = n_iter,
      burn_in = burn_in,
      N       = nrow(Y_pol),
      J       = ncol(Y_pol),
      m       = NA_integer_
    ),
    file = pollen_file
  )
  cat(sprintf("\n  Pollen FMLN done: %.1fs  [saved to fixtures/real_pollen_fmln.rds]\n\n",
              elapsed_pol))
}

# ---------------------------------------------------------------------------
# 2. MLB -- MMLN (player-level random intercepts)
# ---------------------------------------------------------------------------
mlb_file <- file.path(fix_dir, "real_mlb_mmln.rds")

if (file.exists(mlb_file)) {
  cat("  MLB MMLN: SKIP (fixture exists)\n")
} else {
  if (!requireNamespace("Lahman", quietly = TRUE))
    stop("Install the 'Lahman' package first: install.packages('Lahman')")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Install the 'dplyr' package first: install.packages('dplyr')")

  cat("  Loading MLB data via clean_Lahman_data()...\n")
  flush.console()
  mlb <- clean_Lahman_data()
  Y_mlb <- mlb$Y
  X_mlb <- mlb$X
  Z_mlb <- mlb$Z
  m_mlb <- ncol(Z_mlb)

  cat(sprintf("  MLB MMLN | N=%d  J=%d  p=%d  m=%d  n_iter=%d\n",
              nrow(Y_mlb), ncol(Y_mlb), ncol(X_mlb), m_mlb, n_iter))
  flush.console()

  set.seed(99002)
  start <- Sys.time()
  fit_mlb <- MMLN(
    Y        = Y_mlb,
    X        = X_mlb,
    Z        = Z_mlb,
    n_iter   = n_iter,
    burn_in  = burn_in,
    thin     = thin,
    mh_scale = 1,
    proposal = "normbeta",
    verbose  = TRUE
  )
  elapsed_mlb <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  saveRDS(
    list(
      fit     = fit_mlb,
      elapsed = elapsed_mlb,
      n_iter  = n_iter,
      burn_in = burn_in,
      N       = nrow(Y_mlb),
      J       = ncol(Y_mlb),
      m       = m_mlb
    ),
    file = mlb_file
  )
  cat(sprintf("\n  MLB MMLN done: %.1fs  [saved to fixtures/real_mlb_mmln.rds]\n\n",
              elapsed_mlb))
}

cat("Done. Reference fixtures saved to fixtures/\n")
