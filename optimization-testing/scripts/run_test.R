#!/usr/bin/env Rscript
# run_test.R -- Run fork tests for selected functions and scenarios.
#
# Usage: Rscript run_test.R <repo_root> <func> <scenario>
#   func:     "fmln", "mmln", "mdres", or "all"
#   scenario: scenario name, space-separated names, or "all"
#
# Also runs the Tier 4 unit-level mh_update_cpp isolation test when func
# includes fmln or is "all".

cat("=== Running fork tests ===\n\n")

library(MMLN)

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
func_arg  <- if (length(args) >= 2) tolower(args[2]) else "all"
scen_arg  <- if (length(args) >= 3) args[3] else "all"

source(file.path(repo_root, "optimization-testing/config.R"))

data_dir <- file.path(repo_root, "optimization-testing/data")
fix_dir  <- file.path(repo_root, "optimization-testing/fixtures")
res_dir  <- file.path(repo_root, "optimization-testing/results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Parse which functions to run
run_fmln  <- func_arg %in% c("all", "fmln")
run_mmln  <- func_arg %in% c("all", "mmln")
run_mdres <- func_arg %in% c("all", "mdres")

# Parse which scenarios to run
if (scen_arg == "all") {
  scenario_names <- names(SCENARIOS)
} else {
  scenario_names <- strsplit(scen_arg, "\\s+")[[1]]
  invalid <- setdiff(scenario_names, names(SCENARIOS))
  if (length(invalid) > 0) {
    stop("Unknown scenario(s): ", paste(invalid, collapse = ", "),
         "\nValid: ", paste(names(SCENARIOS), collapse = ", "))
  }
}

cat(sprintf("  Functions: %s\n", paste(c("fmln", "mmln", "mdres")[c(run_fmln, run_mmln, run_mdres)], collapse = ", ")))
cat(sprintf("  Scenarios: %s\n\n", paste(scenario_names, collapse = ", ")))

# Helper: format seconds as MM:SS or HH:MM:SS
fmt_time <- function(secs) {
  secs <- round(secs)
  if (secs >= 3600) {
    sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
  } else {
    sprintf("%02d:%02d", secs %/% 60, secs %% 60)
  }
}

# Helper: erase N lines above cursor (ANSI escape codes)
erase_up <- function(n) {
  for (i in seq_len(n)) cat("\033[A\033[2K")
}

# ---------------------------------------------------------------------------
# FMLN and MMLN
# ---------------------------------------------------------------------------
if (run_fmln || run_mmln) {
  cat("--- FMLN/MMLN model runs ---\n\n")

  # Count total jobs for progress
  n_funcs <- sum(run_fmln, run_mmln)
  total_jobs <- length(scenario_names) * length(PROPOSALS) * n_funcs
  job_num <- 0
  phase_start <- Sys.time()

  for (name in scenario_names) {
    sc  <- SCENARIOS[[name]]
    dat <- readRDS(file.path(data_dir, paste0(name, ".rds")))

    for (prop in PROPOSALS) {
      fmln_elapsed <- NA
      mmln_elapsed <- NA
      lines_above  <- 0

      # FMLN
      if (run_fmln) {
        job_num <- job_num + 1
        cat(sprintf("  FMLN | %-20s | %-10s\n", name, prop))
        flush.console()
        lines_above <- lines_above + 1

        start <- Sys.time()
        set.seed(sc$seed + 1000)
        fmln_res <- FMLN(
          Y = dat$Y, X = dat$X,
          n_iter = MCMC$n_iter, burn_in = MCMC$burn_in,
          thin = MCMC$thin, mh_scale = MCMC$mh_scale,
          verbose = TRUE, proposal = prop
        )
        saveRDS(fmln_res, file = file.path(res_dir, paste0("fmln_", name, "_", prop, ".rds")))
        fmln_elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
        lines_above <- lines_above + 1  # progress bar

        # Erase progress bar + header, replace with timed header
        erase_up(lines_above)
        cat(sprintf("\r  FMLN | %-20s | %-10s  %5.1fs\n", name, prop, fmln_elapsed))
        flush.console()
        lines_above <- 1
      }

      # MMLN
      if (run_mmln) {
        job_num <- job_num + 1
        cat(sprintf("  MMLN | %-20s | %-10s\n", name, prop))
        flush.console()
        lines_above <- lines_above + 1

        start <- Sys.time()
        set.seed(sc$seed + 2000)
        mmln_res <- MMLN(
          Y = dat$Y, X = dat$X, Z = dat$Z,
          n_iter = MCMC$n_iter, burn_in = MCMC$burn_in,
          thin = MCMC$thin, mh_scale = MCMC$mh_scale,
          verbose = TRUE, proposal = prop
        )
        saveRDS(mmln_res, file = file.path(res_dir, paste0("mmln_", name, "_", prop, ".rds")))
        mmln_elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
        lines_above <- lines_above + 1  # progress bar
      }

      # Summary line: erase all temp lines, print final
      erase_up(lines_above)
      parts <- c()
      if (!is.na(fmln_elapsed)) parts <- c(parts, sprintf("FMLN: %5.1fs", fmln_elapsed))
      if (!is.na(mmln_elapsed)) parts <- c(parts, sprintf("MMLN: %5.1fs", mmln_elapsed))
      cat(sprintf("\r  [%2d/%d] %-20s %-10s  %s\n",
                  job_num, total_jobs, name, prop, paste(parts, collapse = "  ")))
    }
  }
}

# ---------------------------------------------------------------------------
# MDres (uses precomputed chain fixtures for isolation)
# ---------------------------------------------------------------------------
if (run_mdres) {
  cat("\n--- MDres runs (isolated, using upstream chain fixtures) ---\n\n")

  mdres_total <- length(scenario_names) * length(PROPOSALS)
  mdres_num <- 0
  mdres_start <- Sys.time()

  for (name in scenario_names) {
    sc  <- SCENARIOS[[name]]
    dat <- readRDS(file.path(data_dir, paste0(name, ".rds")))
    n_counts <- rowSums(dat$Y)

    for (prop in PROPOSALS) {
      mdres_num <- mdres_num + 1

      # Load the upstream MMLN fixture to get chain inputs
      mmln_fix_file <- file.path(fix_dir, paste0("mmln_", name, "_", prop, ".rds"))
      if (!file.exists(mmln_fix_file)) {
        cat(sprintf("  [%2d/%d] %-20s %-10s  SKIP (fixture missing)\n", mdres_num, mdres_total, name, prop))
        next
      }

      cat(sprintf("  MDres | %-20s | %-10s\n", name, prop))
      flush.console()

      start <- Sys.time()

      mmln_fix <- readRDS(mmln_fix_file)
      n_saved  <- length(mmln_fix$beta_chain)
      beta_hat  <- mmln_fix$beta_chain[[n_saved]]
      sigma_hat <- mmln_fix$sigma_chain[[n_saved]]
      psi_hat   <- mmln_fix$psi_chain[[n_saved]]

      # Generate P posterior predictive datasets with progress bar
      n_steps <- MDRES_P + 1  # P replicates + 1 MDres call
      pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
      set.seed(sc$seed + 3000 + match(prop, PROPOSALS))
      Y_pred_list <- vector("list", MDRES_P)
      for (j in seq_len(MDRES_P)) {
        Y_pred_list[[j]] <- sample_posterior_predictive(
          X = dat$X, beta = beta_hat, Sigma = sigma_hat,
          n = n_counts, Z = dat$Z, psi = psi_hat,
          mixed = TRUE, verbose = FALSE
        )
        setTxtProgressBar(pb, j)
        elapsed_sec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
        eta_sec <- elapsed_sec / j * (n_steps - j)
        cat(sprintf("\r ETA: %s", fmt_time(eta_sec)))
        flush.console()
      }

      # Run MDres with same seed as reference
      set.seed(sc$seed + 4000 + match(prop, PROPOSALS))
      mdres_out <- MDres(dat$Y, Y_pred_list)
      setTxtProgressBar(pb, n_steps)
      cat(sprintf("\r ETA: %s", fmt_time(0)))
      flush.console()
      close(pb)

      saveRDS(
        list(Y_pred_list = Y_pred_list, mdres = mdres_out),
        file = file.path(res_dir, paste0("mdres_", name, "_", prop, ".rds"))
      )

      elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)

      # Erase progress bar + header, print summary
      erase_up(2)
      cat(sprintf("\r  [%2d/%d] %-20s %-10s  %5.1fs\n", mdres_num, mdres_total, name, prop, elapsed))
    }
  }
}

cat("\nDone. Fork results saved to ", res_dir, "\n", sep = "")
