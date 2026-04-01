#!/usr/bin/env Rscript
# generate_reference.R -- Run upstream MMLN on all scenarios, save fixtures.
#
# Expects:
#   - The upstream MMLN package (eaegerber/MMLN) is already installed
#   - Datasets exist in optimization-testing/data/ (run generate_data.R first)
#
# Saves:
#   - fixtures/fmln_{scenario}_{proposal}.rds   (FMLN chain outputs)
#   - fixtures/mmln_{scenario}_{proposal}.rds   (MMLN chain outputs)
#   - fixtures/mdres_{scenario}.rds             (Y_pred_list + MDres output)
#
# Usage: Rscript generate_reference.R <repo_root> [<func>] [<scenario>]
#   func:     "fmln", "mmln", "mdres", or "all" (default: "all")
#   scenario: scenario name, space-separated names, or "all" (default: "all")
#
# Existing fixture files are skipped (resume support).

cat("=== Generating reference fixtures (upstream) ===\n\n")

library(MMLN)

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
func_arg  <- if (length(args) >= 2) tolower(args[2]) else "all"
scen_arg  <- if (length(args) >= 3) args[3] else "all"

source(file.path(repo_root, "optimization-testing/config.R"))

data_dir <- file.path(repo_root, "optimization-testing/data")
fix_dir  <- file.path(repo_root, "optimization-testing/fixtures")
dir.create(fix_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Parse which functions and scenarios to run
# ---------------------------------------------------------------------------
run_fmln  <- func_arg %in% c("all", "fmln")
run_mmln  <- func_arg %in% c("all", "mmln")
run_mdres <- func_arg %in% c("all", "mdres")

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

# ---------------------------------------------------------------------------
# Check required data files exist
# ---------------------------------------------------------------------------
required_data <- file.path(data_dir, paste0(scenario_names, ".rds"))
missing <- required_data[!file.exists(required_data)]
if (length(missing) > 0) {
  cat("ERROR: Missing data files:\n")
  for (f in missing) cat("  ", basename(f), "\n")
  cat("\nRun 'make generate-data' first.\n")
  quit(status = 1)
}

total_start <- Sys.time()

# ---------------------------------------------------------------------------
# Helper: format seconds as MM:SS or HH:MM:SS
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
# Helper: erase N lines above cursor (ANSI escape codes)
# ---------------------------------------------------------------------------
erase_up <- function(n) {
  for (i in seq_len(n)) cat("\033[A\033[2K")
}

# ---------------------------------------------------------------------------
# Phase 1: FMLN and MMLN for selected scenarios x proposals
# ---------------------------------------------------------------------------
if (run_fmln || run_mmln) {
  cat("--- Phase 1: FMLN & MMLN chains ---\n\n")

  n_funcs <- sum(run_fmln, run_mmln)
  total_jobs <- length(scenario_names) * length(PROPOSALS) * n_funcs
  job_num <- 0

  for (name in scenario_names) {
    sc  <- SCENARIOS[[name]]
    dat <- readRDS(file.path(data_dir, paste0(name, ".rds")))

    for (prop in PROPOSALS) {
      fmln_elapsed <- NA
      mmln_elapsed <- NA
      fmln_skipped <- FALSE
      mmln_skipped <- FALSE
      lines_above  <- 0   # track lines printed that need cleanup

      # --- FMLN ---
      if (run_fmln) {
        job_num <- job_num + 1
        fmln_file <- file.path(fix_dir, paste0("fmln_", name, "_", prop, ".rds"))

        if (file.exists(fmln_file)) {
          fmln_skipped <- TRUE
        } else {
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
          saveRDS(fmln_res, file = fmln_file)
          fmln_elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
          lines_above <- lines_above + 1  # progress bar added a line

          # Erase progress bar + header, replace with timed header
          erase_up(lines_above)
          cat(sprintf("\r  FMLN | %-20s | %-10s  %5.1fs\n", name, prop, fmln_elapsed))
          flush.console()
          lines_above <- 1  # now just 1 line (the timed header)
        }
      }

      # --- MMLN ---
      if (run_mmln) {
        job_num <- job_num + 1
        mmln_file <- file.path(fix_dir, paste0("mmln_", name, "_", prop, ".rds"))

        if (file.exists(mmln_file)) {
          mmln_skipped <- TRUE
        } else {
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
          saveRDS(mmln_res, file = mmln_file)
          mmln_elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
          lines_above <- lines_above + 1  # progress bar added a line
        }
      }

      # --- Summary line ---
      if (fmln_skipped && mmln_skipped) {
        cat(sprintf("  [%2d/%d] %-20s %-10s  SKIP (exists)\n",
                    job_num, total_jobs, name, prop))
      } else {
        erase_up(lines_above)
        parts <- c()
        if (run_fmln) {
          if (fmln_skipped) {
            parts <- c(parts, "FMLN: SKIP")
          } else {
            parts <- c(parts, sprintf("FMLN: %5.1fs", fmln_elapsed))
          }
        }
        if (run_mmln) {
          if (mmln_skipped) {
            parts <- c(parts, "MMLN: SKIP")
          } else {
            parts <- c(parts, sprintf("MMLN: %5.1fs", mmln_elapsed))
          }
        }
        cat(sprintf("\r  [%2d/%d] %-20s %-10s  %s\n",
                    job_num, total_jobs, name, prop, paste(parts, collapse = "  ")))
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Phase 2: MDres for selected scenarios x all MMLN proposals
# ---------------------------------------------------------------------------
if (run_mdres) {
  cat("\n--- Phase 2: MDres (P=", MDRES_P, " replicates) ---\n\n", sep = "")

  mdres_total <- length(scenario_names) * length(PROPOSALS)
  mdres_num <- 0
  mdres_start <- Sys.time()

  for (name in scenario_names) {
    sc  <- SCENARIOS[[name]]
    dat <- readRDS(file.path(data_dir, paste0(name, ".rds")))
    n_counts <- rowSums(dat$Y)

    for (prop in PROPOSALS) {
      mdres_num <- mdres_num + 1
      tag <- paste0("mdres_", name, "_", prop)
      mdres_file <- file.path(fix_dir, paste0(tag, ".rds"))

      # Skip if fixture already exists (resume support)
      if (file.exists(mdres_file)) {
        cat(sprintf("  [%2d/%d] %-20s %-10s  SKIP (exists)\n", mdres_num, mdres_total, name, prop))
        next
      }

      # Load the upstream MMLN fixture to get chain inputs
      mmln_fix_file <- file.path(fix_dir, paste0("mmln_", name, "_", prop, ".rds"))
      if (!file.exists(mmln_fix_file)) {
        cat(sprintf("  [%2d/%d] %-20s %-10s  SKIP (mmln fixture missing)\n", mdres_num, mdres_total, name, prop))
        next
      }

      cat(sprintf("  MDres | %-20s | %-10s\n", name, prop))
      flush.console()

      start <- Sys.time()

      mmln_res <- readRDS(mmln_fix_file)
      n_saved <- length(mmln_res$beta_chain)
      beta_hat  <- mmln_res$beta_chain[[n_saved]]
      sigma_hat <- mmln_res$sigma_chain[[n_saved]]
      psi_hat   <- mmln_res$psi_chain[[n_saved]]

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

      # Run MDres
      set.seed(sc$seed + 4000 + match(prop, PROPOSALS))
      mdres_out <- MDres(dat$Y, Y_pred_list)
      setTxtProgressBar(pb, n_steps)
      cat(sprintf("\r ETA: %s", fmt_time(0)))
      flush.console()
      close(pb)

      saveRDS(
        list(Y_pred_list = Y_pred_list, mdres = mdres_out),
        file = mdres_file
      )

      elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)

      # Erase progress bar + header, print summary
      erase_up(2)
      cat(sprintf("\r  [%2d/%d] %-20s %-10s  %5.1fs\n", mdres_num, mdres_total, name, prop, elapsed))
    }
  }
}

total_elapsed <- round(as.numeric(difftime(Sys.time(), total_start, units = "mins")), 1)
cat(sprintf("\nDone. Fixtures saved to %s (%.1f min total)\n", fix_dir, total_elapsed))
