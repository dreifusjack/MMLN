library(MMLN)
library(mvnfast)

# Minimal Gibbs sampler that runs both psi implementations in lockstep.
# Same seed, same W updates, same beta/Sigma/Phi updates.
# Only psi differs: one uses update_psi_cpp, the other uses the R loop.

set.seed(42)
sim <- simulate_mixed_mln_data(
  m=10, n_i=10, p=3, d=2,
  beta=matrix(c(0.5,-1,0.2,0.3,0.7,-0.4),3,2),
  Sigma=diag(2), Phi=5*diag(2), n_mean=200)

Y <- sim$Y; X <- sim$X; Z <- sim$Z
N <- nrow(Y); d <- ncol(Y)-1; p <- ncol(X); m <- ncol(Z)
n_iter <- 3000; burn_in <- 1000; thin <- 2

psi_r_loop <- function(group_idx, R_tot, Phi, Sigma_inv) {
  psi <- matrix(0, length(group_idx), ncol(R_tot))
  Phi_inv <- chol2inv(chol(Phi))
  for (j in seq_along(group_idx)) {
    idx <- group_idx[[j]] + 1L  # 0-based -> 1-based
    R_j <- R_tot[idx, , drop=FALSE]
    V_j <- solve(Phi_inv + length(idx) * Sigma_inv)
    M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
    psi[j, ] <- mvnfast::rmvn(1, mu=as.vector(M_j), sigma=V_j)
  }
  psi
}

run_gibbs <- function(use_cpp) {
  set.seed(99)
  W         <- alr(compress_counts(Y))
  beta      <- matrix(0, p, d)
  Sigma     <- diag(d); Phi <- diag(d)
  psi       <- matrix(0, m, d)
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv  <- chol2inv(chol(crossprod(X)))
  prior     <- list(beta_var=10, nu_S=d+1, Lambda_S=diag(d),
                    nu_P=d+1, Lambda_P=diag(d))
  group_idx <- build_group_idx(Z)

  keep <- seq(burn_in+1, n_iter, by=thin)
  beta_chain <- vector("list", length(keep)); save_i <- 1

  for (it in seq_len(n_iter)) {
    Mu  <- X %*% beta + Z %*% psi
    ymu <- cbind(Y, Mu)
    W_prop    <- t(apply(ymu, 1, normbetapropdist, Sigma=Sigma))
    ymuw      <- cbind(ymu, W)
    ymuw_new  <- cbind(ymu, W_prop)
    log_q_old <- apply(ymuw,     1, normbetaloglike, Sigma=Sigma)
    log_q_new <- apply(ymuw_new, 1, normbetaloglike, Sigma=Sigma)
    expW  <- rowSums(exp(W));  expWp <- rowSums(exp(W_prop))
    ll_old <- rowSums(Y[,1:d]*W[,1:d])     - rowSums(Y*log1p(expW))  - 0.5*rowSums((W-Mu)%*%Sigma_inv*(W-Mu))
    ll_new <- rowSums(Y[,1:d]*W_prop[,1:d])- rowSums(Y*log1p(expWp)) - 0.5*rowSums((W_prop-Mu)%*%Sigma_inv*(W_prop-Mu))
    ratio  <- ll_new - ll_old + log_q_new - log_q_old
    ratio[is.na(ratio)] <- -Inf
    acc <- log(runif(N)) < ratio
    W[acc,] <- W_prop[acc,]; W[is.na(W)] <- 0

    R_tot <- W - X %*% beta
    psi <- if (use_cpp) update_psi_cpp(group_idx, R_tot, Phi, Sigma_inv)
           else         psi_r_loop(group_idx, R_tot, Phi, Sigma_inv)

    S1 <- prior$Lambda_P + t(psi)%*%psi
    S1 <- (S1+t(S1))/2; diag(S1) <- diag(S1)+1e-8
    W1 <- chol2inv(chol(S1)); W1 <- (W1+t(W1))/2; diag(W1) <- diag(W1)+1e-8
    Phi <- chol2inv(chol(rWishart(1, df=prior$nu_P+m, Sigma=W1)[,,1]))

    R     <- W - Z %*% psi
    beta0 <- S_xx_inv %*% (t(X) %*% R)
    beta  <- matrix(mvnfast::rmvn(1, as.vector(beta0), kronecker(S_xx_inv,Sigma)), nrow=p)

    Eps <- R - X %*% beta
    S2  <- prior$Lambda_S + t(Eps)%*%Eps
    S2  <- (S2+t(S2))/2; diag(S2) <- diag(S2)+1e-8
    W2  <- chol2inv(chol(S2)); W2 <- (W2+t(W2))/2; diag(W2) <- diag(W2)+1e-8
    Sigma <- chol2inv(chol(rWishart(1, df=prior$nu_S+N, Sigma=W2)[,,1]))
    Sigma_inv <- chol2inv(chol(Sigma))

    if (it %in% keep) { beta_chain[[save_i]] <- beta; save_i <- save_i+1 }
  }
  Reduce("+", beta_chain) / length(beta_chain)
}

cat("Running C++ psi...\n")
mean_cpp <- run_gibbs(use_cpp=TRUE)
cat("Running R psi loop...\n")
mean_r <- run_gibbs(use_cpp=FALSE)

cat("\nC++ psi posterior mean:\n"); print(mean_cpp)
cat("\nR psi loop posterior mean:\n"); print(mean_r)
cat("\nTrue beta:\n"); print(sim$beta)
cat("\nMax abs diff between C++ and R:\n"); print(max(abs(mean_cpp - mean_r)))
