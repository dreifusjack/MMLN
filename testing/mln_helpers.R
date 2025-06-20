# working helper functions
# Smithson-Verkuilen Correction
compress_counts <- function(w) {
  if (is.matrix(w)) {
    zero_rows <- apply(w == 0, 1, any)
    rs <- rowSums(w[zero_rows, , drop=FALSE])
    K <- ncol(w)
    w_new <- w
    w_new[zero_rows, ] <- ((w[zero_rows, , drop=FALSE] * (rs - 1)) + 1/K) / rs
    return(w_new)
  } else {
    if (any(w == 0)) {
      N <- sum(w); K <- length(w)
      return((w * (N - 1) + 1/K) / N)
    } else {
      return(w)
    }
  }
}

# Additive Log Ratio
alr <- function(P) {
  if(!is.matrix(P)) P <- matrix(P, nrow=1)
  log(P[,-ncol(P),drop=FALSE] / P[,ncol(P)])
}

# Inverse Additive Log Ratio
alr_inv <- function(Y) {
  if(!is.matrix(Y)) Y <- matrix(Y, nrow=1)
  expY <- exp(Y)
  den <- rowSums(expY) + 1
  cbind(expY/den, 1/den)
}

# MLN Log Likelihood
dmnl_loglik <- function(Y, W) {
  if(!is.matrix(Y)) Y <- matrix(Y, nrow=1)
  if(!is.matrix(W)) W <- matrix(W, nrow=nrow(Y), byrow=TRUE)
  P <- alr_inv(Y)
  sum(W * log(P))
}

# Metropolis-Hastings proposal helpers
# For the Beta proposal
# Logit Function
pstartoy <- function(pstarvec) {
  log(pstarvec / (1 - pstarvec))
}

# Inverse Logit
ytopstar <- function(yvec) {
  exp(yvec) / (1 + exp(yvec))
}

# Beta proposal distribution
betapropdist <- function(WMu_vec, Sigma) {
  k <- ncol(Sigma)
  w_vec <- WMu_vec[1:(k+1)]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  exp_mu <- exp(mu_vec)
  denom <- diag(Sigma)
  alpha <- pmax((1 + exp_mu) / denom - exp_mu / (1 + exp_mu), 1e-3)
  beta <- alpha * exp(-mu_vec)
  alpha_star <- w_vec[1:k] + alpha
  beta_star <- w_vec[k+1] + beta
  rbeta(k, alpha_star, beta_star)
}

# Beta log-likelihood
betaloglike <- function(WMuPstar_vec, Sigma) {
  k <- ncol(Sigma)
  w_vec <- WMuPstar_vec[1:(k+1)]
  mu_vec <- WMuPstar_vec[(k+2):(2*k+1)]
  pstar_vec <- WMuPstar_vec[(2*k+2):length(WMuPstar_vec)]
  exp_mu <- exp(mu_vec)
  denom <- diag(Sigma)
  alpha <- pmax((1 + exp_mu) / denom - exp_mu / (1 + exp_mu), 1e-3)
  beta <- alpha * exp(-mu_vec)
  alpha_star <- w_vec[1:k] + alpha
  beta_star <- w_vec[k+1] + beta
  loglike <- dbeta(pstar_vec, alpha_star, beta_star, log = TRUE)
  logjac <- log(pstar_vec) + log(1 - pstar_vec)
  sum(loglike + logjac)
}

# Normal Approximation to the Beta Proposal Distribution (preferred)
normbetapropdist <- function(WMu_vec, Sigma) {
  k <- ncol(Sigma)
  w_vec <- WMu_vec[1:(k+1)]
  mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
  result <- numeric(length = k)
  for (i in 1:k) {
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i, i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    beta <- alpha * exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k + 1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result[i] <- rnorm(1, muprop, sigprop)
  }
  result
}

# Normal Approximation to the Beta Log-Likelihood
normbetaloglike <- function(WMuY_vec, Sigma) {
  k <- ncol(Sigma)
  w_vec <- WMuY_vec[1:(k+1)]
  mu_vec <- WMuY_vec[(k+2):(2*k+1)]
  y_vec <- WMuY_vec[(2*k+2):length(WMuY_vec)]
  result <- 0
  for (i in 1:k) {
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i, i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
    beta <- alpha * exp(-mu_vec[i])
    alpha_star <- w_vec[i] + alpha
    beta_star <- w_vec[k + 1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result <- result + dnorm(y_vec[i], muprop, sigprop, log = TRUE)
  }
  result
}
