library(profvis)
library(MMLN)

set.seed(42)

# Match similar dataset dimensions to your MMLN profiling
# (20 groups × 20 obs = 400 total, 3 covariates, 3 categories)
sim <- simulate_mixed_mln_data(
  m = 20, n_i = 20, p = 3, d = 2,
  beta = matrix(c(0.5, -1, 0.2, 0.3, 0.7, -0.4), 3, 2),
  Sigma = diag(2),
  Phi = 5 * diag(2),
  n_mean = 200
)

# Profile FMLN with same parameters as your MMLN test
res_profile <- profvis({
  res_f <- FMLN(
    Y        = sim$Y,
    X        = sim$X,
    n_iter   = 500,
    burn_in  = 100,
    thin     = 2,
    proposal = "normbeta",
    verbose  = FALSE
  )
})

# View interactively
res_profile

# Save to HTML for sharing
script_dir <- dirname(sys.frame(1)$ofile)
htmlwidgets::saveWidget(res_profile, file.path(script_dir, "fmln_profile.html"), selfcontained = TRUE)