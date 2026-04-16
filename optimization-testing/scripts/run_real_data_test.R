#!/usr/bin/env Rscript
# run_real_data_test.R -- Run fork FMLN/MMLN/MDres on real pollen and MLB datasets,
#                         then compare timing to upstream fixtures.
#
# MDres uses the upstream MMLN chain (from fixtures/) to isolate the MDres
# implementation from any chain differences — same pattern as the synthetic suite.
#
# Usage: Rscript run_real_data_test.R <repo_root> [<n_iter>] [<burn_in>] [<mdres_p>]

cat("=== Real-data fork test ===\n\n")

library(MMLN)

args      <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
n_iter    <- if (length(args) >= 2) as.integer(args[2]) else 1000L
burn_in   <- if (length(args) >= 3) as.integer(args[3]) else 300L
mdres_p   <- if (length(args) >= 4) as.integer(args[4]) else 50L
thin      <- 2L

fix_dir <- file.path(repo_root, "optimization-testing/fixtures")
res_dir <- file.path(repo_root, "optimization-testing/results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  n_iter=%d  burn_in=%d  thin=%d  mdres_p=%d\n\n",
            n_iter, burn_in, thin, mdres_p))

# ---------------------------------------------------------------------------
# 1. Pollen -- FMLN
# ---------------------------------------------------------------------------
pollen_res_file <- file.path(res_dir, "real_pollen_fmln.rds")
pollen_fix_file <- file.path(fix_dir, "real_pollen_fmln.rds")

if (!requireNamespace("MM", quietly = TRUE))
  stop("Install the 'MM' package first: install.packages('MM')\n  (requires: sudo apt-get install libgmp-dev)")

data(pollen, package = "MM", envir = environment())
Y_pol <- as.matrix(as.data.frame(pollen))
X_pol <- matrix(1, nrow = nrow(Y_pol), ncol = 1)

cat(sprintf("  Pollen FMLN | N=%d  J=%d  n_iter=%d\n", nrow(Y_pol), ncol(Y_pol), n_iter))
flush.console()

set.seed(99001)
t0 <- Sys.time()
fit_pollen <- FMLN(Y = Y_pol, X = X_pol, n_iter = n_iter, burn_in = burn_in,
                   thin = thin, proposal = "normbeta", verbose = TRUE)
elapsed_pol_fmln <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

saveRDS(list(fit = fit_pollen, elapsed = elapsed_pol_fmln,
             n_iter = n_iter, burn_in = burn_in,
             N = nrow(Y_pol), J = ncol(Y_pol), m = NA_integer_),
        file = pollen_res_file)
cat(sprintf("\n  Pollen FMLN done: %.1fs\n\n", elapsed_pol_fmln))

# ---------------------------------------------------------------------------
# 2, 3 & 4. MLB -- load once, run FMLN, MMLN, then MDres
# ---------------------------------------------------------------------------
mlb_fmln_res_file  <- file.path(res_dir, "real_mlb_fmln.rds")
mlb_mmln_res_file  <- file.path(res_dir, "real_mlb_mmln.rds")
mlb_mdres_res_file <- file.path(res_dir, "real_mlb_mdres.rds")
mlb_fmln_fix_file  <- file.path(fix_dir, "real_mlb_fmln.rds")
mlb_mmln_fix_file  <- file.path(fix_dir, "real_mlb_mmln.rds")
mlb_mdres_fix_file <- file.path(fix_dir, "real_mlb_mdres.rds")

if (!requireNamespace("Lahman", quietly = TRUE))
  stop("Install the 'Lahman' package first: install.packages('Lahman')")
if (!requireNamespace("dplyr", quietly = TRUE))
  stop("Install the 'dplyr' package first: install.packages('dplyr')")

cat("  Loading MLB data via clean_Lahman_data()...\n")
flush.console()
mlb      <- clean_Lahman_data()
Y_mlb    <- mlb$Y
X_mlb    <- mlb$X
Z_mlb    <- mlb$Z
m_mlb    <- ncol(Z_mlb)
n_counts <- rowSums(Y_mlb)

# --- MLB FMLN ---
cat(sprintf("  MLB  FMLN | N=%d  J=%d  p=%d  n_iter=%d\n",
            nrow(Y_mlb), ncol(Y_mlb), ncol(X_mlb), n_iter))
flush.console()

set.seed(99003)
t0 <- Sys.time()
fit_mlb_fmln <- FMLN(Y = Y_mlb, X = X_mlb, n_iter = n_iter, burn_in = burn_in,
                     thin = thin, proposal = "normbeta", verbose = TRUE)
elapsed_mlb_fmln <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

saveRDS(list(fit = fit_mlb_fmln, elapsed = elapsed_mlb_fmln,
             n_iter = n_iter, burn_in = burn_in,
             N = nrow(Y_mlb), J = ncol(Y_mlb), m = NA_integer_),
        file = mlb_fmln_res_file)
cat(sprintf("\n  MLB  FMLN done: %.1fs\n\n", elapsed_mlb_fmln))

# --- MLB MMLN ---
cat(sprintf("  MLB  MMLN | N=%d  J=%d  p=%d  m=%d  n_iter=%d\n",
            nrow(Y_mlb), ncol(Y_mlb), ncol(X_mlb), m_mlb, n_iter))
flush.console()

set.seed(99002)
t0 <- Sys.time()
fit_mlb_mmln <- MMLN(Y = Y_mlb, X = X_mlb, Z = Z_mlb, n_iter = n_iter, burn_in = burn_in,
                     thin = thin, mh_scale = 1, proposal = "normbeta", verbose = TRUE)
elapsed_mlb_mmln <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

saveRDS(list(fit = fit_mlb_mmln, elapsed = elapsed_mlb_mmln,
             n_iter = n_iter, burn_in = burn_in,
             N = nrow(Y_mlb), J = ncol(Y_mlb), m = m_mlb),
        file = mlb_mmln_res_file)
cat(sprintf("\n  MLB  MMLN done: %.1fs\n\n", elapsed_mlb_mmln))

# --- MLB MDres (using upstream MMLN chain from fixtures for isolation) ---
if (!file.exists(mlb_mmln_fix_file)) {
  cat("  MLB  MDres: SKIP (upstream MMLN fixture missing -- run real-data-reference first)\n\n")
  elapsed_mlb_mdres <- NA_real_
} else {
  mmln_fix  <- readRDS(mlb_mmln_fix_file)
  n_saved   <- length(mmln_fix$fit$beta_chain)
  beta_hat  <- mmln_fix$fit$beta_chain[[n_saved]]
  sigma_hat <- mmln_fix$fit$sigma_chain[[n_saved]]
  psi_hat   <- mmln_fix$fit$psi_chain[[n_saved]]

  cat(sprintf("  MLB  MDres | P=%d  (using upstream MMLN chain)\n", mdres_p))
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
  elapsed_mlb_mdres <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  saveRDS(list(mdres = mdres_result, Y_pred_list = Y_pred_list,
               elapsed = elapsed_mlb_mdres, mdres_p = mdres_p),
          file = mlb_mdres_res_file)
  cat(sprintf("  MLB  MDres done: %.1fs\n\n", elapsed_mlb_mdres))
}

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

print_row <- function(label, model, fix_file, elapsed_fork, iter = n_iter) {
  if (is.na(elapsed_fork)) {
    cat(sprintf("  %-8s  %-5s  %6s  %10s  %10s  %8s\n",
                label, model, iter, "—", "SKIP", "N/A"))
  } else if (file.exists(fix_file)) {
    fix     <- readRDS(fix_file)
    speedup <- fix$elapsed / elapsed_fork
    cat(sprintf("  %-8s  %-5s  %6s  %9.1fs  %9.1fs  %6.2fx\n",
                label, model, iter, fix$elapsed, elapsed_fork, speedup))
  } else {
    cat(sprintf("  %-8s  %-5s  %6s  %10s  %9.1fs  %8s\n",
                label, model, iter, "NO REF", elapsed_fork, "N/A"))
  }
}

print_row("Pollen", "FMLN",  pollen_fix_file,   elapsed_pol_fmln)
print_row("MLB",    "FMLN",  mlb_fmln_fix_file,  elapsed_mlb_fmln)
print_row("MLB",    "MMLN",  mlb_mmln_fix_file,  elapsed_mlb_mmln)
print_row("MLB",    "MDres", mlb_mdres_fix_file, elapsed_mlb_mdres, mdres_p)

cat("\n")

# Posterior mean comparison (diagnostic -- chains differ by RNG stream, nonzero diff expected)
cat("  Posterior mean comparison (max |upstream - fork|):\n")

beta_diff <- function(fix_file, fit_fork) {
  if (!file.exists(fix_file)) return("NO REF")
  fix       <- readRDS(fix_file)
  n_compare <- min(length(fix$fit$beta_chain), length(fit_fork$beta_chain))
  ref_mean  <- Reduce("+", fix$fit$beta_chain[seq_len(n_compare)]) / n_compare
  frk_mean  <- Reduce("+", fit_fork$beta_chain[seq_len(n_compare)]) / n_compare
  sprintf("%.6f", max(abs(ref_mean - frk_mean)))
}

cat(sprintf("    Pollen FMLN beta: %s\n", beta_diff(pollen_fix_file,  fit_pollen)))
cat(sprintf("    MLB    FMLN beta: %s\n", beta_diff(mlb_fmln_fix_file, fit_mlb_fmln)))
cat(sprintf("    MLB    MMLN beta: %s\n", beta_diff(mlb_mmln_fix_file, fit_mlb_mmln)))

cat("\n")
cat("  Fork results saved to results/\n")
cat("=================================================================\n")
