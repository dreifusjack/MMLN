library(profvis)
library(MMLN)

# 0. Load the cleaned MLB data (post-1960, min. PA of 200, four categories: HR, BB, SO, Other)
baseball_example <- clean_Lahman_data()

# 1. The fixed model
fmln_res_profile <- profvis({
  mlb_f <- FMLN(
    Y = baseball_example$Y,
    X = cbind(1, baseball_example$X),
    n_iter = 100,
    burn_in = 30,
    proposal = "normbeta",
    mh_scale = .2,
    verbose = TRUE
  )
})

# 2. Adding in random effect (will take longer)
mmln_res_profile <- profvis({
  mlb_m <- MMLN(
    Y = baseball_example$Y,
    X = cbind(1, baseball_example$X),
    Z = baseball_example$Z,
    n_iter = 100,
    burn_in = 30,
    proposal = "normbeta",
    mh_scale = .2,
    verbose = TRUE
  )
})

script_dir <- dirname(sys.frame(1)$ofile)
htmlwidgets::saveWidget(fmln_res_profile, file.path(script_dir, "fmln_baseball_profile.html"), selfcontained = TRUE)
htmlwidgets::saveWidget(mmln_res_profile, file.path(script_dir, "mmln_baseball_profile.html"), selfcontained = TRUE)

# # 3. Posterior predictive simulation and Mahalanobis residuals
# Y_pred_list_f <- lapply(seq_along(mlb_f$w_chain), function(i) {
#   sample_posterior_predictive(X = cbind(1, baseball_example$X),
#                               beta = mlb_f$beta_chain[[i]],
#                               Sigma = mlb_f$sigma_chain[[i]],
#                               n = baseball_example$PA,
#                               mixed = FALSE,
#                               verbose = FALSE
#   )
# })
# resids_f <- MDres(baseball_example$Y, Y_pred_list_f)
# summary(resids_f)


# Y_pred_list_m <- lapply(seq_along(mlb_m$w_chain), function(i) {
#   sample_posterior_predictive(X = cbind(1, baseball_example$X),
#                               beta = mlb_m$beta_chain[[i]],
#                               Sigma = mlb_m$sigma_chain[[i]],
#                               n = baseball_example$PA,
#                               Z = baseball_example$Z,
#                               psi = mlb_m$psi_chain[[i]],
#                               mixed = TRUE,
#                               verbose = FALSE
#   )
# })
# resids_m <- MDres(baseball_example$Y, Y_pred_list_m)
# summary(resids_m)

# # 4. Evaluate posterior parameter means of the fixed effects model
# post_beta <- apply(simplify2array(mlb_f$beta_chain), c(1,2), mean)
# row.names(post_beta) <- c("int", colnames(baseball_example$X))
# colnames(post_beta) <- colnames(baseball_example$Y)[1:3]
# post_beta

# # Note: first column represents covariate effect on HR relative to Other
# # Example: Right handed hitters tend to hit more HR than Left handed hitters
# #          NL batters tend to walk more than AL batters
# #          Taller batters tend SO more than shorter batters
# #          The Average Batter (avg. weight, height, age, Lefty, AL, C)

# avg_batter_pi <- alr_inv(post_beta[1,])
# avg_batter_pi