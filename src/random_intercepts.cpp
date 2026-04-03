// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// Replaces the per-group random intercept update loop in MMLN.
//
// The R loop this replaces (runs every MCMC iteration):
//   for(j in seq_len(m)) {
//     idx <- which(Z[, j] == 1)
//     R_j <- R_tot[idx, , drop = FALSE]
//     V_j <- solve(chol2inv(chol(Phi)) + length(idx) * Sigma_inv)
//     M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
//     psi[j, ] <- rmvn(1, mu = as.vector(M_j), sigma = V_j)
//   }
//
// Key changes vs the R version:
//   1. Group indices passed in pre-computed. Z never changes across iterations.
//   2. Eigen LLT for Phi_inv. Armadillo inv_sympd for V_j (avoids Arma/Eigen round-trip in loop).
//   3. Phi_inv computed once outside the group loop.
//   4. Armadillo Cholesky MVN sampler for psi draws.

//' compute_psi_moments_cpp: Return posterior M_j and V_j for each group.
//'
//' Deterministic test helper. No random draws. Use this to verify the C++
//' posterior mean and covariance match the R reference to 1e-8.
//'
//' @param group_idx List of 0-based integer vectors from build_group_idx
//' @param R_tot     N x d residual matrix: W - X %*% beta
//' @param Phi       d x d random intercept covariance
//' @param Sigma_inv d x d inverse within-group covariance
//' @return List with M (list of d-vectors) and V (list of d x d matrices)
//' @export
// [[Rcpp::export]]
Rcpp::List compute_psi_moments_cpp(
    const Rcpp::List& group_idx,
    const arma::mat& R_tot,
    const arma::mat& Phi,
    const arma::mat& Sigma_inv
) {
  int m = group_idx.size();
  int d = R_tot.n_cols;

  Eigen::MatrixXd Phi_e = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Phi));
  Eigen::LLT<Eigen::MatrixXd> llt_phi(Phi_e);
  Eigen::MatrixXd Phi_inv_e = llt_phi.solve(Eigen::MatrixXd::Identity(d, d));
  arma::mat Phi_inv = Rcpp::as<arma::mat>(Rcpp::wrap(Phi_inv_e));

  Rcpp::List M_list(m), V_list(m);
  for (int j = 0; j < m; j++) {
    arma::uvec idx = Rcpp::as<arma::uvec>(group_idx[j]);
    int n_j = idx.n_elem;
    arma::mat A   = Phi_inv + n_j * Sigma_inv;
    arma::mat V_j = arma::inv_sympd(A);
    arma::mat R_j = R_tot.rows(idx);
    arma::vec s_j = arma::trans(arma::sum(R_j, 0));
    arma::vec M_j = V_j * (Sigma_inv * s_j);
    M_list[j] = Rcpp::wrap(M_j);
    V_list[j] = Rcpp::wrap(V_j);
  }
  return Rcpp::List::create(Rcpp::Named("M") = M_list,
                            Rcpp::Named("V") = V_list);
}

// Draw one sample from MVN(mu, Sigma) using Cholesky decomposition.
// Internal helper. Not exported to R.
static arma::rowvec rmvn_chol(const arma::vec& mu, const arma::mat& Sigma) {
  int d = Sigma.n_rows;
  arma::mat L = arma::chol(Sigma, "lower");
  arma::vec z = arma::randn<arma::vec>(d);
  return (mu + L * z).t();
}

//' Pre-compute group membership indices from Z.
//'
//' Call this once before the MCMC loop. Pass the result to update_psi_cpp
//' each iteration. Z never changes, so recomputing inside update_psi_cpp
//' every iteration wastes cycles.
//'
//' @param Z  N x m group indicator matrix
//' @return List of integer vectors, one per group, with 0-based row indices
//' @export
// [[Rcpp::export]]
Rcpp::List build_group_idx(const arma::mat& Z) {
  int m = Z.n_cols;
  Rcpp::List idx_list(m);
  for (int j = 0; j < m; j++) {
    // Store as Rcpp::IntegerVector for easy R interop
    arma::uvec idx = arma::find(Z.col(j) > 0.5);
    idx_list[j] = Rcpp::wrap(idx);
  }
  return idx_list;
}

//' update_psi_cpp: Update random intercepts for all m groups.
//'
//' Call build_group_idx once before the MCMC loop and pass group_idx here.
//' Uses Eigen LLT for Phi_inv. Uses Armadillo inv_sympd for V_j inside
//' the loop to avoid Arma/Eigen round-trip copies per group.
//'
//' @param group_idx List of 0-based integer vectors from build_group_idx
//' @param R_tot     N x d residual matrix: W - X %*% beta
//' @param Phi       d x d random intercept covariance
//' @param Sigma_inv d x d inverse within-group covariance
//' @return m x d matrix of updated random intercepts (psi)
//' @export
// [[Rcpp::export]]
arma::mat update_psi_cpp(
    const Rcpp::List& group_idx,
    const arma::mat& R_tot,
    const arma::mat& Phi,
    const arma::mat& Sigma_inv
) {
  int m = group_idx.size();
  int d = R_tot.n_cols;

  // Compute Phi_inv once. Use Eigen LLT for the small d x d solve.
  Eigen::MatrixXd Phi_e = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Phi));
  Eigen::LLT<Eigen::MatrixXd> llt_phi(Phi_e);
  Eigen::MatrixXd Phi_inv_e = llt_phi.solve(Eigen::MatrixXd::Identity(d, d));
  arma::mat Phi_inv = Rcpp::as<arma::mat>(Rcpp::wrap(Phi_inv_e));

  arma::mat psi(m, d);

  for (int j = 0; j < m; j++) {
    arma::uvec idx = Rcpp::as<arma::uvec>(group_idx[j]);
    int n_j = idx.n_elem;

    // Posterior covariance: V_j = inv(Phi_inv + n_j * Sigma_inv).
    // Stay in Armadillo -- avoids Arma/Eigen copy per group iteration.
    arma::mat A = Phi_inv + n_j * Sigma_inv;
    arma::mat V_j = arma::inv_sympd(A);

    // Posterior mean: M_j = V_j * (Sigma_inv * colSums(R_j)).
    arma::mat R_j = R_tot.rows(idx);
    arma::vec s_j = arma::trans(arma::sum(R_j, 0));
    arma::vec M_j = V_j * (Sigma_inv * s_j);

    // Draw psi[j,] ~ MVN(M_j, V_j).
    psi.row(j) = rmvn_chol(M_j, V_j);
  }

  return psi;
}
