#!/usr/bin/env Rscript
# run_real_data_test.R -- Run fork MMLN/FMLN on real pollen and MLB datasets,
#                         then compare timing and posterior means to upstream fixtures.
#
# Expects:
#   - The fork MMLN package is already installed (done by Makefile before this script)
#   - fixtures/real_pollen_fmln.rds and fixtures/real_mlb_mmln.rds exist
#     (run `make real-data-reference` first)
#
# Usage: Rscript run_real_data_test.R <repo_root> [<n_iter>] [<burn_in>]
#   n_iter:  MCMC iterations (default: 2000)
#   burn_in: burn-in to discard (default: 500)

cat("=== Real-data fork test ===\n\n")

library(MMLN)

args      <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
n_iter    <- if (length(args) >= 2) as.integer(args[2]) else 2000L
burn_in   <- if (length(args) >= 3) as.integer(args[3]) else 500L
thin      <- 2L

fix_dir <- file.path(repo_root, "optimization-testing/fixtures")
res_dir <- file.path(repo_root, "optimization-testing/results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  n_iter=%d  burn_in=%d  thin=%d\n\n", n_iter, burn_in, thin))

# Collected timing for summary table
timing <- list()

# ---------------------------------------------------------------------------
# 1. Pollen -- FMLN
# ---------------------------------------------------------------------------
pollen_res_file <- file.path(res_dir, "real_pollen_fmln.rds")
pollen_fix_file <- file.path(fix_dir, "real_pollen_fmln.rds")

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
  file = pollen_res_file
)
cat(sprintf("\n  Pollen FMLN done: %.1fs\n\n", elapsed_pol))
timing$pollen <- elapsed_pol

# ---------------------------------------------------------------------------
# 2. MLB -- MMLN
# ---------------------------------------------------------------------------
mlb_res_file <- file.path(res_dir, "real_mlb_mmln.rds")
mlb_fix_file <- file.path(fix_dir, "real_mlb_mmln.rds")

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
  file = mlb_res_file
)
cat(sprintf("\n  MLB MMLN done: %.1fs\n\n", elapsed_mlb))
timing$mlb <- elapsed_mlb

# ---------------------------------------------------------------------------
# Summary comparison table
# ---------------------------------------------------------------------------
cat("\n")
cat("=================================================================\n")
cat("  Real Data Benchmark\n")
cat("=================================================================\n")
cat(sprintf("  %-8s  %-5s  %6s  %10s  %10s  %8s\n",
            "Dataset", "Model", "n_iter", "Reference", "Fork", "Speedup"))
cat("  -------- -----  ------  ----------  ----------  --------\n")

# Pollen comparison
if (file.exists(pollen_fix_file)) {
  fix_pol <- readRDS(pollen_fix_file)
  speedup_pol <- fix_pol$elapsed / elapsed_pol
  cat(sprintf("  %-8s  %-5s  %6d  %9.1fs  %9.1fs  %6.2fx\n",
              "Pollen", "FMLN", n_iter, fix_pol$elapsed, elapsed_pol, speedup_pol))
} else {
  cat(sprintf("  %-8s  %-5s  %6d  %10s  %9.1fs  %8s\n",
              "Pollen", "FMLN", n_iter, "NO REF", elapsed_pol, "N/A"))
}

# MLB comparison
if (file.exists(mlb_fix_file)) {
  fix_mlb <- readRDS(mlb_fix_file)
  speedup_mlb <- fix_mlb$elapsed / elapsed_mlb
  cat(sprintf("  %-8s  %-5s  %6d  %9.1fs  %9.1fs  %6.2fx\n",
              "MLB", "MMLN", n_iter, fix_mlb$elapsed, elapsed_mlb, speedup_mlb))
} else {
  cat(sprintf("  %-8s  %-5s  %6d  %10s  %9.1fs  %8s\n",
              "MLB", "MMLN", n_iter, "NO REF", elapsed_mlb, "N/A"))
}

cat("\n")

# Posterior mean comparison (diagnostic -- RNG streams differ so nonzero diff is expected)
cat("  Posterior mean comparison (max |upstream - fork|):\n")

if (file.exists(pollen_fix_file)) {
  fix_pol <- readRDS(pollen_fix_file)
  n_saved_ref <- length(fix_pol$fit$beta_chain)
  n_saved_frk <- length(fit_pollen$beta_chain)
  n_compare   <- min(n_saved_ref, n_saved_frk)
  beta_ref_mean <- Reduce("+", fix_pol$fit$beta_chain[seq_len(n_compare)]) / n_compare
  beta_frk_mean <- Reduce("+", fit_pollen$beta_chain[seq_len(n_compare)]) / n_compare
  diff_pol <- max(abs(beta_ref_mean - beta_frk_mean))
  cat(sprintf("    Pollen FMLN beta: %.6f\n", diff_pol))
} else {
  cat("    Pollen FMLN beta: NO REF (run real-data-reference first)\n")
}

if (file.exists(mlb_fix_file)) {
  fix_mlb_loaded <- readRDS(mlb_fix_file)
  n_saved_ref <- length(fix_mlb_loaded$fit$beta_chain)
  n_saved_frk <- length(fit_mlb$beta_chain)
  n_compare   <- min(n_saved_ref, n_saved_frk)
  beta_ref_mean <- Reduce("+", fix_mlb_loaded$fit$beta_chain[seq_len(n_compare)]) / n_compare
  beta_frk_mean <- Reduce("+", fit_mlb$beta_chain[seq_len(n_compare)]) / n_compare
  diff_mlb <- max(abs(beta_ref_mean - beta_frk_mean))
  cat(sprintf("    MLB    MMLN beta: %.6f\n", diff_mlb))
} else {
  cat("    MLB    MMLN beta: NO REF (run real-data-reference first)\n")
}

cat("\n")
cat("  Fork results saved to results/\n")
cat("=================================================================\n")
