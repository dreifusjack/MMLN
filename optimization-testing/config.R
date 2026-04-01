# config.R -- Scenario definitions and tolerance constants for upstream vs fork comparison
#
# Each scenario is a named list consumed by generate_data.R, generate_reference.R,
# run_test.R, and compare.R. Tolerances are used by compare.R.
#
# Scenarios are scaled to match the MLB baseball example from the README:
#   N ~ 2000 obs, m ~ 500 groups, p = 10, d = 3, n_iter = 1000

# ---------------------------------------------------------------------------
# Tolerance constants
# ---------------------------------------------------------------------------
TOL_EXACT   <- 0       # pure-R proposals ("norm", "beta") must be bitwise identical
TOL_FLOAT   <- 1e-8    # single-operation C++ vs R precision (digamma, exp, ...)
TOL_CHAIN   <- 1e-6    # cumulative drift after 1000 MCMC iterations
TOL_ACCEPT  <- 1e-10   # acceptance ratios (boolean-derived, nearly exact)
TOL_MDRES   <- 1e-6    # MDres residuals (downstream from chain differences)

# ---------------------------------------------------------------------------
# MCMC settings shared across all scenarios
# ---------------------------------------------------------------------------
MCMC <- list(
  n_iter   = 1000,
  burn_in  = 300,
  thin     = 2,
  mh_scale = 1
)

# Proposal types to test for each scenario
PROPOSALS <- c("normbeta", "norm", "beta")

# MDres settings
MDRES_P <- 50  # number of posterior predictive replicates

# ---------------------------------------------------------------------------
# Helper: near-singular positive-definite matrix
# ---------------------------------------------------------------------------
make_near_singular <- function(d, rho = 0.92) {
  S <- matrix(rho, d, d)
  diag(S) <- 1.0
  S
}

# ---------------------------------------------------------------------------
# Helper: resolve repo root from command-line or source context
# ---------------------------------------------------------------------------
get_repo_root <- function() {
  # Check if passed as first command-line arg
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && dir.exists(args[1])) {
    return(normalizePath(args[1]))
  }
  # Try Rscript --file= invocation
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    # scripts/ is one level below optimization-testing/, which is one level below repo
    return(normalizePath(file.path(script_dir, "../..")))
  }
  # Try source() invocation
  for (i in rev(seq_len(sys.nframe()))) {
    if (!is.null(sys.frame(i)$ofile)) {
      return(normalizePath(file.path(dirname(sys.frame(i)$ofile), "../..")))
    }
  }
  # Fallback
  normalizePath(".")
}

# ---------------------------------------------------------------------------
# Scenario definitions
#
# Baseball reference: N~2000, m~500, p=10, d=3(J=4), n_mean~200
# All scenarios target N~1200-2000 to match that computational load.
# ---------------------------------------------------------------------------
SCENARIOS <- list(

  # 1. Baseball-scale sanity check
  # N=2000, m=50, p=10, d=3 -- mirrors the MLB example dimensions
  baseline = list(
    seed   = 42001,
    m      = 50,
    n_i    = 40,
    p      = 10,
    d      = 3,
    n_mean = 200,
    beta   = matrix(c(
    #  int    wt    ht    age   age2  batL  batR  lgAL  posIF posOF
       0.5, -0.3,  0.2,  0.4, -0.1,  0.3, -0.2,  0.1,  0.15, -0.25,
      -0.5,  0.8, -0.1, -0.3,  0.2, -0.4,  0.5, -0.1,  0.20,  0.10,
       0.3, -0.4,  0.6,  0.1, -0.3,  0.2, -0.1,  0.3, -0.10,  0.30
    ), nrow = 10, ncol = 3),
    Sigma  = diag(3),
    Phi    = 0.5 * diag(3)
  ),

  # 2. Many zero cells -- 8 categories with low counts guarantees zeros everywhere
  # N=1500, m=30, p=5, d=7
  zero_counts = list(
    seed   = 42002,
    m      = 30,
    n_i    = 50,
    p      = 5,
    d      = 7,
    n_mean = 10,
    beta   = matrix(c(
       0.5,  0.0,  0.0,  0.0,  0.0,
       0.0,  0.5,  0.0,  0.0,  0.0,
       0.0,  0.0,  0.3,  0.0,  0.0,
       0.0,  0.0,  0.0,  0.2,  0.0,
       0.0,  0.0,  0.0,  0.0,  0.1,
      -0.1,  0.0,  0.0,  0.0,  0.0,
       0.0, -0.1,  0.0,  0.0,  0.0
    ), nrow = 5, ncol = 7),
    Sigma  = 0.3 * diag(7),
    Phi    = 0.2 * diag(7)
  ),

  # 3. Very small total counts -- extreme ALR values, wide W variance
  # N=1200, m=40, p=5, d=3
  tiny_counts = list(
    seed   = 42003,
    m      = 40,
    n_i    = 30,
    p      = 5,
    d      = 3,
    n_mean = 5,
    beta   = matrix(c(
      -0.5,  0.3,  0.1,  0.2, -0.1,
       0.4, -0.2,  0.3, -0.1,  0.2,
       0.1,  0.5, -0.2,  0.0,  0.3
    ), nrow = 5, ncol = 3),
    Sigma  = diag(3),
    Phi    = 0.3 * diag(3)
  ),

  # 4. Large total counts -- precision in exp()/log1p() accumulates
  # N=1200, m=30, p=5, d=3
  large_counts = list(
    seed   = 42004,
    m      = 30,
    n_i    = 40,
    p      = 5,
    d      = 3,
    n_mean = 5000,
    beta   = matrix(c(
       0.3, -0.2,  0.1,  0.2, -0.1,
      -0.1,  0.4, -0.2,  0.1,  0.3,
       0.2, -0.3,  0.1, -0.2,  0.0
    ), nrow = 5, ncol = 3),
    Sigma  = diag(3),
    Phi    = 0.5 * diag(3)
  ),

  # 5. High-dimensional -- 9 outcome categories (d=8), 8x8 covariance
  # N=1200, m=30, p=8, d=8
  high_dim = list(
    seed   = 42005,
    m      = 30,
    n_i    = 40,
    p      = 8,
    d      = 8,
    n_mean = 100,
    beta   = matrix(c(
       0.10, -0.20,  0.15, -0.10,  0.05,  0.20, -0.15,  0.10,
      -0.25,  0.05,  0.30, -0.20,  0.10, -0.05,  0.25, -0.10,
      -0.10,  0.20, -0.15,  0.25, -0.05, -0.30,  0.10,  0.15,
       0.25, -0.05, -0.30,  0.15,  0.20,  0.10, -0.05,  0.20,
       0.20,  0.10, -0.20,  0.05, -0.15,  0.25,  0.30, -0.25,
      -0.15,  0.25,  0.05,  0.30, -0.10,  0.05, -0.20,  0.10,
       0.30, -0.10,  0.10, -0.25,  0.15, -0.20,  0.05,  0.30,
      -0.05,  0.15, -0.25,  0.10,  0.25,  0.15, -0.10, -0.05
    ), nrow = 8, ncol = 8),
    Sigma  = diag(8),
    Phi    = diag(8)
  ),

  # 6. Near-singular covariance -- rho=0.92, condition number ~100
  # N=1400, m=40, p=5, d=4
  near_singular = list(
    seed   = 42006,
    m      = 40,
    n_i    = 35,
    p      = 5,
    d      = 4,
    n_mean = 200,
    beta   = matrix(c(
       0.3, -0.1,  0.2,  0.1, -0.2,
      -0.2,  0.4, -0.1,  0.3,  0.0,
       0.1, -0.3,  0.2, -0.1,  0.3,
       0.2,  0.1, -0.2,  0.0,  0.1
    ), nrow = 5, ncol = 4),
    Sigma  = make_near_singular(4, 0.92),
    Phi    = make_near_singular(4, 0.85)
  ),

  # 7. Large coefficients -- overflow risk in exp(W)
  # N=1200, m=30, p=8, d=3
  large_coefs = list(
    seed   = 42007,
    m      = 30,
    n_i    = 40,
    p      = 8,
    d      = 3,
    n_mean = 100,
    beta   = matrix(c(
       3.0, -4.0,  2.0, -1.5,  3.5, -2.5,  1.0, -3.0,
      -3.0,  4.0, -2.0,  2.5, -3.5,  1.5, -1.0,  2.0,
       1.0, -3.0,  2.0, -2.0,  1.5, -1.0,  3.0, -4.0
    ), nrow = 8, ncol = 3),
    Sigma  = 2 * diag(3),
    Phi    = diag(3)
  ),

  # 8. Unbalanced groups -- n_i from 1 to 80, total N~2000
  # N~2025, m=50, p=5, d=3
  unbalanced = list(
    seed   = 42008,
    m      = 50,
    n_i    = c(1, 1, 1, 2, 2, 3, 3, 4, 5, 5,
               6, 7, 8, 10, 10, 12, 14, 15, 18, 20,
               22, 25, 28, 30, 32, 35, 38, 40, 42, 45,
               48, 50, 52, 55, 58, 60, 60, 62, 65, 65,
               68, 68, 70, 70, 72, 72, 75, 75, 78, 80),
    p      = 5,
    d      = 3,
    n_mean = 150,
    beta   = matrix(c(
       0.4, -0.3,  0.2,  0.1, -0.2,
      -0.2,  0.5, -0.1,  0.3,  0.0,
       0.3, -0.1,  0.2, -0.2,  0.1
    ), nrow = 5, ncol = 3),
    Sigma  = diag(3),
    Phi    = 0.5 * diag(3)
  ),

  # 9. Many groups, few observations -- mirrors baseball's ~500 players
  # N=2000, m=500, p=5, d=2
  many_groups = list(
    seed   = 42009,
    m      = 500,
    n_i    = 4,
    p      = 5,
    d      = 2,
    n_mean = 100,
    beta   = matrix(c(
       0.5, -0.3,  0.2,  0.1, -0.2,
      -0.4,  0.6, -0.1,  0.3,  0.0
    ), nrow = 5, ncol = 2),
    Sigma  = diag(2),
    Phi    = diag(2)
  ),

  # 10. Dominant category -- reference category gets ~90%+ probability
  # N=1200, m=30, p=5, d=2
  dominant_cat = list(
    seed   = 42010,
    m      = 30,
    n_i    = 40,
    p      = 5,
    d      = 2,
    n_mean = 200,
    beta   = matrix(c(
      -3.0,  0.0,  0.1,  0.0, -0.1,
      -3.0,  0.0, -0.1,  0.0,  0.1
    ), nrow = 5, ncol = 2),
    Sigma  = 0.3 * diag(2),
    Phi    = 0.2 * diag(2)
  )
)
