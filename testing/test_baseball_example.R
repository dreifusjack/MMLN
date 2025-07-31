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



# Now, see how long it takes to fit...
library(MMLN)

baseball_example <- clean_Lahman_data()

# let's try it, both fixed and mixed
mlb_f <- FMLN(
  Y = baseball_example$Y,
  X = cbind(1, baseball_example$X),
  n_iter = 500,
  burn_in = 100,
  proposal = "normbeta",
  verbose = TRUE
)

# adding in random effect (may take much longer)
mlb_m <- MMLN(
  Y = baseball_example$Y,
  X = cbind(1, baseball_example$X),
  Z = baseball_example$Z,
  n_iter = 500,
  burn_in = 100,
  proposal = "normbeta",
  verbose = TRUE
)
# 3. Posterior predictive simulation and Mahalanobis residuals
Y_pred_list_f <- lapply(seq_along(mlb_f$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_f$beta_chain[[i]],
                              Sigma = mlb_f$sigma_chain[[i]],
                              n = baseball_example$PA,
                              mixed = FALSE
  )
})
resids_f <- MDres(baseball_example$Y, Y_pred_list_f)
summary(resids_f)


Y_pred_list_m <- lapply(seq_along(mlb_m$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_m$beta_chain[[i]],
                              Sigma = mlb_m$sigma_chain[[i]],
                              n = baseball_example$PA,
                              Z = baseball_example$Z,
                              psi = mlb_m$psi_chain[[i]],
                              mixed = TRUE
  )
})
resids_m <- MDres(baseball_example$Y, Y_pred_list_m)
summary(resids_m)
