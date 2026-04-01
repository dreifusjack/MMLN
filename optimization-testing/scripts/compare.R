#!/usr/bin/env Rscript
# compare.R -- Load fixtures (upstream) and fork results, run assertions, produce report.
#
# No MMLN dependency needed -- only reads .rds files and does arithmetic.
# Exit code 1 on any failure.
#
# Usage: Rscript compare.R <repo_root> <func> <scenario>

cat("=== Comparing upstream fixtures vs fork results ===\n\n")

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath("../..")
func_arg  <- if (length(args) >= 2) tolower(args[2]) else "all"
scen_arg  <- if (length(args) >= 3) args[3] else "all"

source(file.path(repo_root, "optimization-testing/config.R"))

fix_dir <- file.path(repo_root, "optimization-testing/fixtures")
res_dir <- file.path(repo_root, "optimization-testing/results")
rep_dir <- file.path(repo_root, "optimization-testing/report")
dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)

# Parse args
run_fmln  <- func_arg %in% c("all", "fmln")
run_mmln  <- func_arg %in% c("all", "mmln")
run_mdres <- func_arg %in% c("all", "mdres")

if (scen_arg == "all") {
  scenario_names <- names(SCENARIOS)
} else {
  scenario_names <- strsplit(scen_arg, "\\s+")[[1]]
}

all_pass <- TRUE
report_rows <- list()

# ---------------------------------------------------------------------------
# Helper: compare two chain lists element-by-element
# ---------------------------------------------------------------------------
compare_chains <- function(chain_a, chain_b, tol, label) {
  n <- length(chain_a)
  if (n != length(chain_b)) {
    return(list(label = label, max_abs = Inf, max_rel = Inf,
                pass = FALSE, trajectory = NULL, note = "length mismatch"))
  }
  if (n == 0) {
    return(list(label = label, max_abs = 0, max_rel = 0,
                pass = TRUE, trajectory = NULL, note = "empty"))
  }

  max_abs <- 0
  max_rel <- 0
  trajectory <- numeric(n)

  for (i in seq_len(n)) {
    a <- as.numeric(chain_a[[i]])
    b <- as.numeric(chain_b[[i]])
    if (length(a) != length(b)) {
      return(list(label = label, max_abs = Inf, max_rel = Inf,
                  pass = FALSE, trajectory = NULL,
                  note = sprintf("dim mismatch at iter %d", i)))
    }
    abs_diff <- abs(a - b)
    iter_max <- max(abs_diff)
    trajectory[i] <- iter_max
    max_abs <- max(max_abs, iter_max)
    denom <- pmax(abs(a), abs(b), 1e-15)
    max_rel <- max(max_rel, max(abs_diff / denom))
  }

  list(label = label, max_abs = max_abs, max_rel = max_rel,
       pass = max_abs <= tol, trajectory = trajectory, note = "")
}

# ---------------------------------------------------------------------------
# Helper: check for exponential divergence
# ---------------------------------------------------------------------------
check_divergence_growth <- function(trajectory) {
  n <- length(trajectory)
  if (n < 10) return(TRUE)
  mid <- floor(n / 2)
  first_half  <- max(trajectory[1:mid])
  second_half <- max(trajectory[(mid + 1):n])
  if (first_half > 0 && second_half / first_half > 100) return(FALSE)
  TRUE
}

# ---------------------------------------------------------------------------
# Helper: add a row to the report
# ---------------------------------------------------------------------------
add_row <- function(scenario, proposal, model, chain, max_abs, max_rel, tol, pass) {
  report_rows[[length(report_rows) + 1]] <<- data.frame(
    scenario = scenario, proposal = proposal, model = model,
    chain = chain, max_abs_diff = max_abs, max_rel_diff = max_rel,
    tol = tol, pass = pass, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# FMLN / MMLN chain comparison
# ---------------------------------------------------------------------------
FMLN_CHAINS <- c("beta_chain", "sigma_chain", "w_chain", "mhaccept_chain")
MMLN_CHAINS <- c("beta_chain", "sigma_chain", "w_chain", "mhaccept_chain",
                 "phi_chain", "psi_chain")

models_to_check <- c()
if (run_fmln) models_to_check <- c(models_to_check, "fmln")
if (run_mmln) models_to_check <- c(models_to_check, "mmln")

if (length(models_to_check) > 0) {
  cat("--- Tiers 1-3: Chain comparison ---\n\n")

  for (name in scenario_names) {
    for (model_name in models_to_check) {
      chains <- if (model_name == "fmln") FMLN_CHAINS else MMLN_CHAINS

      for (prop in PROPOSALS) {
        tag <- paste0(model_name, "_", name, "_", prop)
        fix_file <- file.path(fix_dir, paste0(tag, ".rds"))
        res_file <- file.path(res_dir, paste0(tag, ".rds"))

        if (!file.exists(fix_file) || !file.exists(res_file)) {
          cat(sprintf("  SKIP: %s (files missing)\n", tag))
          next
        }

        fix_data <- readRDS(fix_file)
        res_data <- readRDS(res_file)

        is_cpp <- (prop == "normbeta")

        for (ch in chains) {
          tol <- if (!is_cpp) TOL_EXACT
                 else if (ch == "mhaccept_chain") TOL_ACCEPT
                 else TOL_CHAIN

          result <- compare_chains(fix_data[[ch]], res_data[[ch]], tol,
                                   paste(name, prop, toupper(model_name), ch))

          status <- if (result$pass) "PASS" else "FAIL"
          cat(sprintf("  %-55s max|diff|=%.2e  tol=%.0e  %s\n",
                      result$label, result$max_abs, tol, status))

          if (!result$pass) all_pass <- FALSE

          # Tier 3: divergence trajectory check
          if (is_cpp && !is.null(result$trajectory) && length(result$trajectory) > 10) {
            growth_ok <- check_divergence_growth(result$trajectory)
            if (!growth_ok) {
              cat(sprintf("    *** EXPONENTIAL DIVERGENCE in %s ***\n", result$label))
              all_pass <- FALSE
            }
            tryCatch({
              png_file <- file.path(rep_dir, paste0("divergence_", tag, "_", ch, ".png"))
              png(png_file, width = 600, height = 300)
              plot(result$trajectory, type = "l", main = result$label,
                   xlab = "Saved iteration", ylab = "max|diff|",
                   col = ifelse(growth_ok, "darkgreen", "red"))
              abline(h = tol, lty = 2, col = "gray50")
              dev.off()
            }, error = function(e) NULL)
          }

          add_row(name, prop, toupper(model_name), ch,
                  result$max_abs, result$max_rel, tol, result$pass)
        }
      }
    }
  }
}

# ---------------------------------------------------------------------------
# MDres comparison
# ---------------------------------------------------------------------------
if (run_mdres) {
  cat("\n--- MDres comparison ---\n\n")

  for (name in scenario_names) {
    for (prop in PROPOSALS) {
      fix_file <- file.path(fix_dir, paste0("mdres_", name, "_", prop, ".rds"))
      res_file <- file.path(res_dir, paste0("mdres_", name, "_", prop, ".rds"))

      if (!file.exists(fix_file) || !file.exists(res_file)) {
        cat(sprintf("  SKIP: mdres_%s_%s (files missing)\n", name, prop))
        next
      }

      fix_data <- readRDS(fix_file)
      res_data <- readRDS(res_file)

      label <- paste("mdres", name, prop)

      # Compare Y_pred_list element-by-element
      ypred_max_diff <- 0
      for (j in seq_along(fix_data$Y_pred_list)) {
        d <- max(abs(fix_data$Y_pred_list[[j]] - res_data$Y_pred_list[[j]]))
        ypred_max_diff <- max(ypred_max_diff, d)
      }
      ypred_pass <- ypred_max_diff <= TOL_CHAIN
      cat(sprintf("  %-45s max|diff|=%.2e  tol=%.0e  %s\n",
                  paste(label, "Y_pred_list"), ypred_max_diff, TOL_CHAIN,
                  ifelse(ypred_pass, "PASS", "FAIL")))
      if (!ypred_pass) all_pass <- FALSE
      add_row(name, prop, "MDRES", "Y_pred_list",
              ypred_max_diff, NA, TOL_CHAIN, ypred_pass)

      # Compare MDres residuals (ignoring NAs in both)
      fix_r <- fix_data$mdres
      res_r <- res_data$mdres
      both_valid <- !is.na(fix_r) & !is.na(res_r)
      if (sum(both_valid) > 0) {
        mdres_diff <- max(abs(fix_r[both_valid] - res_r[both_valid]))
      } else {
        mdres_diff <- 0
      }
      # Check NA pattern matches
      na_match <- identical(is.na(fix_r), is.na(res_r))
      mdres_pass <- mdres_diff <= TOL_MDRES && na_match

      cat(sprintf("  %-45s max|diff|=%.2e  tol=%.0e  %s%s\n",
                  paste(label, "residuals"), mdres_diff, TOL_MDRES,
                  ifelse(mdres_pass, "PASS", "FAIL"),
                  ifelse(na_match, "", " (NA pattern mismatch)")))
      if (!mdres_pass) all_pass <- FALSE
      add_row(name, prop, "MDRES", "residuals",
              mdres_diff, NA, TOL_MDRES, mdres_pass)
    }
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
report_df <- do.call(rbind, report_rows)

cat("\n=== SUMMARY ===\n")
n_pass <- sum(report_df$pass)
n_fail <- sum(!report_df$pass)
cat(sprintf("  Total checks: %d  |  PASS: %d  |  FAIL: %d\n",
            nrow(report_df), n_pass, n_fail))

if (n_fail > 0) {
  cat("\nFailed checks:\n")
  fails <- report_df[!report_df$pass, ]
  for (i in seq_len(nrow(fails))) {
    cat(sprintf("  %s / %s / %s / %s  max|diff|=%.2e  tol=%.0e\n",
                fails$scenario[i], fails$proposal[i], fails$model[i],
                fails$chain[i], fails$max_abs_diff[i], fails$tol[i]))
  }
}

# Save report
summary_file <- file.path(rep_dir, "summary.txt")
sink(summary_file)
cat("MMLN Optimization Test Report\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("FUNC=%s  SCENARIO=%s\n\n", func_arg, scen_arg))
print(report_df, row.names = FALSE)
cat(sprintf("\nOverall: %s\n", ifelse(all_pass, "ALL PASS", "FAILURES DETECTED")))
sink()

saveRDS(report_df, file = file.path(rep_dir, "summary.rds"))
cat(sprintf("\nReport saved to %s\n", summary_file))

if (!all_pass) {
  cat("\n*** COMPARISON FAILED ***\n")
  quit(status = 1)
} else {
  cat("\nAll comparisons passed.\n")
}
