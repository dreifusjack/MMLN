#!/usr/bin/env Rscript
# generate_data.R -- Create all synthetic test datasets with fixed seeds.
#
# Sources mln_helpers.R directly (only needs mvnfast) so it does not depend
# on any particular MMLN package installation.
#
# Usage: Rscript generate_data.R <repo_root>

cat("=== Generating synthetic test datasets ===\n")

if (!requireNamespace("mvnfast", quietly = TRUE)) {
  stop("Package 'mvnfast' is required. Install with: install.packages('mvnfast')")
}

# Resolve paths
args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")

source(file.path(repo_root, "optimization-testing/config.R"))
source(file.path(repo_root, "R/mln_helpers.R"))

data_dir <- file.path(repo_root, "optimization-testing/data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

for (name in names(SCENARIOS)) {
  sc <- SCENARIOS[[name]]
  n_i_str <- if (length(sc$n_i) > 1) "varied" else as.character(sc$n_i)
  cat(sprintf("  %-15s  m=%-3d  n_i=%-6s  d=%d  n_mean=%d\n",
              name, sc$m, n_i_str, sc$d, sc$n_mean))

  set.seed(sc$seed)
  dat <- simulate_mixed_mln_data(
    m      = sc$m,
    n_i    = sc$n_i,
    p      = sc$p,
    d      = sc$d,
    beta   = sc$beta,
    Sigma  = sc$Sigma,
    Phi    = sc$Phi,
    n_mean = sc$n_mean
  )

  saveRDS(dat, file = file.path(data_dir, paste0(name, ".rds")))
}

cat(sprintf("\nDone. %d datasets saved to %s\n", length(SCENARIOS), data_dir))
