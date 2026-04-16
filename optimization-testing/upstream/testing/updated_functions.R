# ------------------------------------------------------------------------------
# Test Script for Comparing Old vs New Versions of Utility Functions
# ------------------------------------------------------------------------------

####Packages################################
##  Packages  ##
library(car)          #for confidence envelope in qqPlot
library(doParallel)   #for parallel computing
library(dirmult)      #for estimating dirichlet parameters
library(EnvStats)     #for uniform qqPlot
library(iterators)    #for some parallel options
library(MCMCpack)     #for inverse wishart sampling
library(MGLM)         #for multinomial logistic regression with counts
library(mvnfast)      #for multivariate normal sampling

# 1. Deterministic Inverse via Cholesky
#    We can verify cholinv() just by comparing with solve() or cross-checking
#Inverse of matrix using Cholesky Decomposition (faster than solve())
cholinv <- function(x){
  chol2inv(chol(x))
}

set.seed(123)
test_mat <- matrix(rnorm(16), 4, 4)
test_mat <- crossprod(test_mat)  # make it positive definite
myInv1 <- cholinv(test_mat)
myInv2 <- solve(test_mat)

# Should be numerically identical:
all.equal(myInv1, myInv2, tolerance = 1e-12)


# 2. Compare old_pstartoy() and new pstartoy()
old_pstartoy <- function(pstarvec){
  y <- numeric(length=length(pstarvec))
  for(i in 1:length(pstarvec)){
    y[i] <- log(pstarvec[i]/(1-pstarvec[i]))
  }
  return(y)
}

pstartoy <- function(pstarvec){
  log(pstarvec / (1 - pstarvec))
}

pstar_test <- c(.2, .3, .5)

old_ps <- old_pstartoy(pstar_test)
new_ps <- pstartoy(pstar_test)
cat("old_pstartoy vs. pstartoy:\n")
print(old_ps)
print(new_ps)
cat("Are they identical? =>", all.equal(old_ps, new_ps), "\n\n")


# 3. Compare old_ytopstar() and new ytopstar()
old_ytopstar <- function(yvec){
  pstar <- numeric(length=length(yvec))
  for(i in 1:length(yvec)){
    pstar[i] <- exp(yvec[i])/(1 + exp(yvec[i]))
  }
  return(pstar)
}

ytopstar <- function(yvec){
  expy <- exp(yvec)
  expy / (1 + expy)
}

y_test <- old_pstartoy(pstar_test)  # some log-odds
old_pstar <- old_ytopstar(y_test)
new_pstar <- ytopstar(y_test)

cat("old_ytopstar vs. ytopstar:\n")
print(old_pstar)
print(new_pstar)
cat("Are they identical? =>", all.equal(old_pstar, new_pstar), "\n\n")


# 4. Compare old_betapropdist() and new betapropdist() with same seed
old_betapropdist <- function(WMu_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMu_vec[1:(k+1)]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  result <- numeric(length=k)
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) -
      (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    if(alpha < 0){
      alpha <- alpha + 1
    }
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star  <- w_vec[k+1] + beta
    result[i]  <- rbeta(1, alpha_star, beta_star)
  }
  return(result)
}

betapropdist <- function(WMu_vec, Sigma){
  k <- ncol(Sigma)
  w_vec  <- WMu_vec[1:k]
  w_base <- WMu_vec[k+1]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  
  exp_mu     <- exp(mu_vec)
  diagSigma  <- diag(Sigma)
  alpha_vec  <- (1 + exp_mu)/diagSigma - (exp_mu/(1 + exp_mu))
  
  alpha_vec[alpha_vec < 0] <- alpha_vec[alpha_vec < 0] + 1
  
  beta_vec   <- alpha_vec * exp(-mu_vec)
  alpha_star <- w_vec + alpha_vec
  beta_star  <- w_base + beta_vec
  
  rbeta(k, alpha_star, beta_star)
}

# Example inputs
set.seed(999)
Sigma_test <- crossprod(matrix(rnorm(4), 2, 2))  # 2x2 PSD
test_vec <- c(0.5, 0.7, 0.1, 0.9, -1.0)     # length = k+1 + k = 2+1 + 2 = 5

# The old approach calls rbeta(1, ...) k times, the new approach calls rbeta(k, ...).
# We expect the same seed => same outcome if the parameters line up in the same order.
set.seed(1234)
old_draw <- old_betapropdist(test_vec, Sigma_test)

set.seed(1234)
new_draw <- betapropdist(test_vec, Sigma_test)

cat("old_betapropdist vs. betapropdist:\n")
cat("Old draw: ", paste0(round(old_draw, 5), collapse = ", "), "\n")
cat("New draw: ", paste0(round(new_draw, 5), collapse = ", "), "\n")
cat("All.equal? =>", all.equal(old_draw, new_draw, tol = 1e-15), "\n\n")

# --------------------------------------------------------------------
# Test Script: Old vs. New betaloglike()
# --------------------------------------------------------------------

# Old version
old_betaloglike <- function(WMuPstar_vec, Sigma){
  k <- ncol(Sigma)
  w_vec <- WMuPstar_vec[1:(k+1)]
  mu_vec <- WMuPstar_vec[(k+2):(2*k+1)]
  pstar_vec <- WMuPstar_vec[(2*k+2):length(WMuPstar_vec)]
  loglike <- 0
  logjac <- 0
  for(i in 1:k){
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) -
      (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    if(alpha < 0){
      alpha <- alpha + 1
    }
    beta <- alpha*exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star  <- w_vec[k+1] + beta
    
    loglike <- loglike + dbeta(pstar_vec[i], alpha_star, beta_star, log = TRUE)
    logjac  <- logjac  + log(pstar_vec[i]) + log(1 - pstar_vec[i])
  }
  return(loglike + logjac)
}

# New vectorized version
betaloglike_vec <- function(WMuPstar_vec, Sigma){
  k <- ncol(Sigma)
  
  w_vec     <- WMuPstar_vec[1:k]
  w_base    <- WMuPstar_vec[k+1]
  mu_vec    <- WMuPstar_vec[(k+2):(2*k+1)]
  pstar_vec <- WMuPstar_vec[(2*k+2):length(WMuPstar_vec)]
  
  diagSigma <- diag(Sigma)
  
  exp_mu    <- exp(mu_vec)
  alpha_vec <- (1 + exp_mu)/diagSigma - (exp_mu / (1 + exp_mu))
  
  alpha_vec[alpha_vec < 0] <- alpha_vec[alpha_vec < 0] + 1
  
  beta_vec       <- alpha_vec * exp(-mu_vec)
  alpha_star_vec <- w_vec + alpha_vec
  beta_star_vec  <- w_base + beta_vec
  
  loglike <- sum(dbeta(pstar_vec, alpha_star_vec, beta_star_vec, log = TRUE))
  logjac  <- sum(log(pstar_vec) + log(1 - pstar_vec))
  
  loglike + logjac
}

# --------------------------------------------------------------------
# Compare outputs on random inputs
# --------------------------------------------------------------------
set.seed(999)

num_tests <- 2  # run multiple tests
for(tt in seq_len(num_tests)){
  # Random dimension
  k <- sample(2:5, 1)  # pick a dimension between 2 and 5
  # Build a random Sigma
  A <- matrix(rnorm(k^2), nrow = k)
  Sigma_test <- crossprod(A)  # ensures positive-definite
  
  # Construct a random WMuPstar_vec with length = 3*k + 1
  # old code indexing: 
  #  w_vec => k+1,
  #  mu_vec => k,
  #  pstar => k
  # total = k+1 + k + k = 3k + 1
  WMuPstar_test <- rnorm(3*k + 1)
  # But we want pstar in (0,1). Let's forcibly do that:
  # last k entries are pstar:
  pstar_indices <- (2*k+2):(3*k+1)
  WMuPstar_test[pstar_indices] <- runif(k, 0.01, 0.99)
  
  # Evaluate old vs. new
  old_val <- old_betaloglike(WMuPstar_test, Sigma_test)
  new_val <- betaloglike_vec(WMuPstar_test, Sigma_test)
  
  are_equal <- all.equal(old_val, new_val, tol = 1e-12)
  cat(sprintf("Test %d: k=%d => old=%.6f, new=%.6f, match=%s\n", 
              tt, k, old_val, new_val, are_equal))
}




















# End of test script
cat("Done. All tests completed.\n")
