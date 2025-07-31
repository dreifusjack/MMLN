#' Real-Data Example: Pollen Data - MN, DM & MLN Models
#'
#' This function reproduces the pollen-data example from Gerber & Craig (2024), fitting three
#' models to the "pollen" counts (multinomial logit, Dirichlet-multinomial, and
#' multinomial-logistic-normal), drawing replicate from fitted sampling distributions, and computing
#' Mahalanobis residuals for each model.
#'
#' @param n_iter   Integer; total number of MCMC iterations for the MLN model (default 1000)
#' @param burn_in  Integer; number of initial MLN iterations to discard (default 400)
#' @param thin     Integer; MLN thinning interval (default 2)
#' @param proposal Character; one of `"norm"`, `"beta"`, or `"normbeta"` for the MLN sampler (default `"normbeta"`)
#' @param P        Integer; number of fitted sampling distribution replicates per model (default 1000)
#'
#' @return A list with components:
#' \describe{
#'   \item{fit_mlr}{The `MGLMreg` object for the multinomial logit fit.}
#'   \item{fit_dm}{The `MGLMreg` object for the Dirichlet-multinomial fit.}
#'   \item{fit_mln}{The list returned by `FMLN()` for the fixed-effects MLN fit.}
#'   \item{Y_pred_mlr}{List of length P of \eqn{N \times J} count matrices sampled from the MN model.}
#'   \item{Y_pred_dm}{List of length P of \eqn{N \times J} count matrices sampled from the DM model.}
#'   \item{Y_pred_mln}{List of length P of \eqn{N \times J} count matrices sampled from the MLN model (using posterior-mean parameters).}
#'   \item{resids_mlr}{Mahalanobis residuals (`mdres`) for the multinomial logit model.}
#'   \item{resids_dm}{Mahalanobis residuals (`mdres`) for the Dirichlet-multinomial model.}
#'   \item{resids_mln}{Mahalanobis residuals (`mdres`) for the MLN model.}
#' }
#'
#' @examples
#' \dontrun{
#' # run all three fits & diagnostics
#' pollen_res <- run_pollen_models(n_iter = 1500, burn_in = 500, thin = 5, P = 500)
#'
#' # view KS-tests & QQ-plots
#' summary(pollen_res$resids_mlr)
#' summary(pollen_res$resids_dm)
#' summary(pollen_res$resids_mln)
#' }
#'
#' @importFrom stats rmultinom
#' @importFrom mvnfast rmvn
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console data
#' @importFrom MGLM MGLMreg rdirmn
#' @export
run_pollen_models <- function(n_iter   = 1000,
                              burn_in  = 400,
                              thin     = 2,
                              proposal = "normbeta",
                              P        = 1000) {
  ## dependencies
  if (!requireNamespace("MM",   quietly = TRUE)) stop("Install the 'MM' package to load pollen data")
  if (!requireNamespace("MGLM", quietly = TRUE)) stop("Install the 'MGLM' package for MGLMreg()")

  ## 1) load data
  data(pollen, package = "MM", envir = environment())
  pollen <- as.data.frame(pollen)
  Y  <- as.matrix(pollen)
  N  <- nrow(Y)
  n <- rowSums(Y)
  J  <- ncol(Y)

  X <- matrix(1, nrow = N, ncol = 1)  # intercept-only design

  ## 2) fit M-Logit (multinomial) and Dirichlet-multinomial via MGLM
  message("? Fitting multinomial-logistic model using MGLMreg()")
  fit_mlr <- suppressWarnings(MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1,
                     data = as.data.frame(pollen), dist = "MN"))
  cat("Done: \n")
  message("? Fitting Dirichlet-lmultinomial model using MGLMreg()")
  fit_dm  <- suppressWarnings(MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1,
                     data = as.data.frame(pollen), dist = "DM"))

  ## 3) fit fixed-effects MLN
  cat("Done: \n")
  message("? Fitting multinomial-logistic-normal model using FMLN()")
  fit_mln <- FMLN(
    Y              = Y,
    X              = X,
    n_iter         = n_iter,
    burn_in        = burn_in,
    thin           = thin,
    proposal       = proposal,
    verbose        = TRUE
  )

  ## 4) prepare predictive replicates
  cat("Done: \n")
  message("? Drawing sampling distributions of predicted values from each model")
  # 4a) MN model: draw P full-dataset replicates
  probs_mlr <- fit_mlr@fitted
  Y_pred_mlr <- vector("list", P)
  for (p_i in seq_len(P)) {
    M <- t(sapply(seq_len(N),
                  function(i) rmultinom(1, size = n[i], prob = probs_mlr[i, ])))
    Y_pred_mlr[[p_i]] <- M
  }

  # 4b) DM model: draw P replicates with Dirichlet-multinomial
  alpha_hat <- exp(fit_dm@coefficients)
  Y_pred_dm <- vector("list", P)
  for (p_i in seq_len(P)) {
    M <- t(sapply(seq_len(N),
                  function(i) MGLM::rdirmn(n    = 1,
                                           size = n[i],
                                           alpha = alpha_hat)))
    Y_pred_dm[[p_i]] <- M
  }

  # 4c) MLN model: use posterior-mean ? and ? to draw P replicates
  # posterior means:
  beta_arr  <- simplify2array(fit_mln$beta_chain)    # p ? d ? n_saves
  Sigma_arr <- simplify2array(fit_mln$sigma_chain)   # d ? d ? n_saves
  beta_mean  <- apply(beta_arr,  c(1,2), mean)
  Sigma_mean <- apply(Sigma_arr, c(1,2), mean)

  Y_pred_mln <- vector("list", P)
  for (p_i in seq_len(P)) {
    Y_pred_mln[[p_i]] <- sample_posterior_predictive(
      X     = X,
      beta  = beta_mean,
      Sigma = Sigma_mean,
      n    = n,
      mixed = FALSE
    )
  }

  ## 5) compute Mahalanobis-residuals
  cat("Done: \n")  # ensure we start on a fresh line
  message("? Computing MDRes for multinomial-logit model")
  resids_mlr <- MDres(Y, Y_pred_mlr)

  cat("\nDone: \n")  # ensure we start on a fresh line
  message("? Computing MDRes for Dirichlet-multinomial model")
  resids_dm  <- MDres(Y, Y_pred_dm)

  cat("\nDone: \n")  # ensure we start on a fresh line
  message("? Computing MDRes for multinomial-logistic-normal model")
  resids_mln <- MDres(Y, Y_pred_mln)

  ## return all results
  output = list(
    fit_mlr    = fit_mlr,
    fit_dm     = fit_dm,
    fit_mln    = fit_mln,
    Y_pred_mlr = Y_pred_mlr,
    Y_pred_dm  = Y_pred_dm,
    Y_pred_mln = Y_pred_mln,
    resids_mlr = resids_mlr,
    resids_dm  = resids_dm,
    resids_mln = resids_mln
  )

  return(output)
}

#' Prepare Cleaned Lahman Baseball Data for Multinomial Modeling
#'
#' This function cleans and processes data from the \pkg{Lahman} database to produce a usable dataset
#' for modeling MLB player batting outcomes. It returns multinomial response data (`Y`), standardized predictors (`X`),
#' and a random intercept design matrix (`Z`) for player-level effects.
#'
#' The output includes only non-pitcher players with PA > 200 from the American and National Leagues
#' (post-1960; expansion era), and aggregates data by player-season (grouping across stints and teams).
#' It collapses batting outcomes into home runs, walks, strikeouts, and other outcomes.
#' Predictors include physical measurements, age (and age^2), batting side, league, and a grouped fielding position indicator.
#' Investigating the code, it should not be difficult to adjust the outcomes investigated and/or the predictors.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y}{An \eqn{N \times J} matrix of batting outcome counts (default: HR, BB, SO, and Other).}
#'   \item{X}{An \eqn{N \times p} numeric matrix of standardized predictors (default: p = 10).}
#'   \item{Z}{An \eqn{N \times m} indicator matrix for random intercepts by playerID.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- clean_Lahman_data()
#' str(data$X)
#' colnames(data$Z)[1:5]
#' }
#'
#' @importFrom Lahman Batting People Fielding
#' @importFrom stats model.matrix
#' @importFrom dplyr mutate group_by summarise filter arrange recode select
#' @export

clean_Lahman_data <- function(){
  ## dependency
  if (!requireNamespace("Lahman", quietly = TRUE)) stop("Install the 'Lahman' package to get MLB data")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install the 'dplyr' package for easy data manipulation")

  ## 1) load data
  data(Batting, package = "Lahman", envir = environment())
  data(People, package = "Lahman", envir = environment())
  data(Fielding, package = "Lahman", envir = environment())

  # Let's only care about AL and NL, and post 1900 (Modern Era) or 1960 (Expansion Era)
  # Fiddle with this for data size
  year_cutoff <- 1960  # or try 1900 for a larger data set

  # Compute debut year in AL/NL for each player
  player_debut <- Batting %>%
    filter(lgID %in% c("AL","NL")) %>%
    group_by(playerID) %>%
    summarise(debut = min(yearID), .groups="drop")

  # Filter to only those who debuted AFTER the cutoff
  valid_players <- player_debut %>%
    filter(debut > year_cutoff) %>%
    pull(playerID)

  # Continue with building data set
  MLB_Batting <- Batting %>%
    filter(playerID %in% valid_players,
           lgID     %in% c("AL","NL"),
           yearID   >  year_cutoff)

  MLB_People <- People  %>% filter(playerID %in% valid_players)

  # Now, merge all People information into MLB_Batting
  MLB_Merged <- merge(MLB_Batting, MLB_People, by = "playerID")

  # Sort by name, year, stint
  MLB_Merged <- MLB_Merged %>%
    arrange(playerID, yearID, stint)

  # Separate into Info, X, and y data frames
  # Info (keep some information which may be useful to refer to later; i.e. team or birthCountry)
  MLB_Info <- MLB_Merged[,c("playerID", "nameFirst", "nameLast", "yearID", "teamID", "lgID", "birthDate", "birthCountry")]

  # X (predictors; will need to create some and may not use all)
  MLB_Merged$age <- MLB_Merged$yearID - MLB_Merged$birthYear
  MLB_X <- MLB_Merged[,c("weight", "height", "age", "bats", "throws", "lgID")]

  # Want to add Position information (need to use Fielding)
  # Filter Fielding the same way as above
  MLB_Fielding <- Fielding %>%
    filter(playerID %in% valid_players,
           lgID %in% c("AL","NL"),
           yearID > year_cutoff)
  # Replace missing InnOuts with 0
  MLB_Fielding <- MLB_Fielding %>%
    mutate(InnOuts = ifelse(is.na(InnOuts), 0, InnOuts))
  # Get the most played position per player-year-stint (will assume any NA are DH/PH in a sec)
  Primary_Pos <- MLB_Fielding %>%
    mutate(InnOuts = coalesce(InnOuts,0)) %>%
    group_by(playerID, yearID, stint) %>%
    slice_max(InnOuts, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    select(playerID, yearID, stint, POS)
  # Add to MLB_Merged
  MLB_Merged_Pos <- MLB_Merged %>%
    left_join(Primary_Pos, by = c("playerID", "yearID", "stint"))
  # Finally put it into MLB_X and then update DH/PH for the NA accordingly
  MLB_X$pos <- MLB_Merged_Pos$POS
  MLB_X$pos[is.na(MLB_X$pos) & (
    MLB_Merged_Pos$yearID == 2020 | # COVID
      (MLB_Merged_Pos$lgID == "AL" & MLB_Merged_Pos$yearID >= 1973) | # AL DH
      (MLB_Merged_Pos$lgID == "NL" & MLB_Merged_Pos$yearID >= 2023) # NL DH
  )] <- "DH"
  MLB_X$pos[is.na(MLB_X$pos)] <- "PH"

  # y (output; may combine to reduce number of categories, but must sum up to PA)
  MLB_y <- cbind("X1B" = MLB_Merged$H - rowSums(MLB_Merged[,c("X2B", "X3B", "HR")]),
                 MLB_Merged[,c("X2B", "X3B", "HR", "BB", "SO", "HBP", "SH", "SF")])
  MLB_y <- cbind(MLB_y, "OIP" = rowSums(MLB_Merged[,c("AB", "BB", "HBP", "SH", "SF")]) - rowSums(MLB_y))

  # save the PA (plate appearances; exposure for the Multinomial)
  PA <- rowSums(MLB_y)

  # Turn the data frames into useable matrices
  # Let's focus on the three true outcomes (HR, BB, SO, Other)
  Y <- as.matrix(cbind(MLB_y[,c("HR", "BB", "SO")], "Other" = PA - rowSums(MLB_y[,c("HR", "BB", "SO")])))

  # For batting, throws is probably not important, so remove it
  # Add indicator variables for the categorical columns
  # Add a quadratic age term
  # drop unused levels first
  MLB_X$lgID <- droplevels(factor(MLB_X$lgID))
  # some unknown batting hands
  MLB_X$bats <- as.character(MLB_X$bats)
  MLB_X$bats[is.na(MLB_X$bats)] <- "U" # only necessary for post-1900, not post-1960
  MLB_X$bats[MLB_X$bats == "B"] <- "S"
  MLB_X$bats <- factor(MLB_X$bats)
  MLB_X$bats <- droplevels(factor(MLB_X$bats))
  # for ease, combine Positions into C, INF, OF, P, DH
  MLB_X$pos_group <- recode(MLB_X$pos,
                            "C" = "C",
                            "1B" = "IF", "2B" = "IF", "SS" = "IF", "3B" = "IF",
                            "OF" = "OF",
                            "P" = "P",
                            "DH" = "DPH", "PH" = "DPH",
                            .default = "C")
  MLB_X$pos_group <- factor(MLB_X$pos_group)
  # scale numeric features
  MLB_numeric <- cbind(MLB_X[,c("weight", "height", "age")],
                       "age2" = MLB_X$age^2)

  X <- as.matrix(cbind(scale(MLB_numeric),
                       model.matrix(~ bats + lgID + pos_group, data = MLB_X)[,-1]))

  # Now we need Z for the random effects; these are the player intercepts
  # based on playerID
  Z <- model.matrix(~ playerID - 1, data = MLB_Merged)

  # To make it take less time, and actually fit, Update Y, X, and Z to:
  # If a player played in the same league in the same year, combine those rows
  # Ignore pitchers
  # Keep only rows where PA > 30
  MLB_Combo <- cbind(MLB_Info[,c("playerID", "yearID", "lgID")],
                     PA = PA,
                     Y,
                     X)
  MLB_Combo$pos <- MLB_X$pos
  MLB_Combo <- MLB_Combo %>%
    filter(pos != "P", PA >= 200) # adjust PA for bigger/smaller data sets

  MLB_Combo <- MLB_Combo %>%
    group_by(playerID, yearID, lgID) %>%
    summarise(
      across(c(HR, BB, SO, Other, PA), sum),
      across(where(is.numeric) & !c(HR, BB, SO, Other, PA), mean),
      .groups = "drop"
    )

  MLB_Combo2 <- MLB_Combo %>%
    group_by(playerID) %>%      # group by player
    filter(n() > 1) %>%         # keep only players with more than one row
    ungroup()

  # New matrices (ignore batsU and pos_groupP, since they are empty now; though don't need to do batsU if post-1960)
  Y <- as.matrix(MLB_Combo2[, c("HR", "BB", "SO", "Other")])
  PA <- MLB_Combo2$PA
  if(year_cutoff == 1900){
    X <- as.matrix(MLB_Combo2 %>% select(-playerID, -yearID, -lgID, -HR, -BB, -SO, -Other, -PA, -batsU, -pos_groupP))
  } else{
    X <- as.matrix(MLB_Combo2 %>% select(-playerID, -yearID, -lgID, -HR, -BB, -SO, -Other, -PA, -pos_groupP))
  }

  Z <- model.matrix(~ playerID - 1, data = MLB_Combo2)

  output = list(
    Y = Y,
    X = X,
    Z = Z,
    PA = PA
  )

  return(output)

}

