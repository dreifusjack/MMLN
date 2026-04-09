// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppEigen.h>
#pragma GCC diagnostic pop

// Pre-compute 0-based row indices for each group column of Z.
// Called once before the MCMC loop — Z does not change across iterations.
//
// @param Z  N x m group indicator matrix (0/1 entries)
// @return   Rcpp::List of length m, each element an integer vector of 0-based row indices
//
// [[Rcpp::export]]
Rcpp::List build_group_idx(const arma::mat &Z)
{
    const int N = static_cast<int>(Z.n_rows);
    const int m = static_cast<int>(Z.n_cols);

    Rcpp::List out(m);
    for (int j = 0; j < m; ++j) {
        std::vector<int> idx;
        idx.reserve(N);
        for (int i = 0; i < N; ++i) {
            if (Z(i, j) == 1.0)
                idx.push_back(i);
        }
        out[j] = Rcpp::wrap(idx);
    }
    return out;
}

// Compute posterior mean M_j and covariance V_j for each group without drawing.
// Used by acceptance tests to verify math to 1e-8 without RNG noise.
//
// Math (per group j with n_j observations):
//   Phi_inv = inv(Phi)                        [computed once via Eigen LLT]
//   A_j     = Phi_inv + n_j * Sigma_inv
//   V_j     = inv(A_j)                        [via Eigen LLT]
//   s_j     = colSums(R_tot[idx_j, ])
//   M_j     = V_j %*% (Sigma_inv %*% s_j)
//
// @param group_idx  Rcpp::List of 0-based index vectors (from build_group_idx)
// @param R_tot      N x d residual matrix (W - X %*% beta)
// @param Phi        d x d random-intercept covariance
// @param Sigma_inv  d x d inverse within-group covariance
// @return Rcpp::List with "M" (list of d-vectors) and "V" (list of d x d matrices)
//
// [[Rcpp::export]]
Rcpp::List compute_psi_moments_cpp(
    const Rcpp::List   &group_idx,
    const arma::mat    &R_tot,
    const arma::mat    &Phi,
    const arma::mat    &Sigma_inv)
{
    const int m = static_cast<int>(group_idx.size());
    const int d = static_cast<int>(R_tot.n_cols);

    // Compute Phi_inv once via Eigen LLT (matches chol2inv(chol(Phi)) in R)
    Eigen::MatrixXd Phi_eig = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Phi));
    Eigen::LLT<Eigen::MatrixXd> llt_phi(Phi_eig);
    Eigen::MatrixXd Phi_inv_eig = llt_phi.solve(Eigen::MatrixXd::Identity(d, d));

    // Convert Sigma_inv to Eigen for A_j construction
    Eigen::MatrixXd Sig_inv_eig = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Sigma_inv));

    Rcpp::List M_list(m);
    Rcpp::List V_list(m);

    for (int j = 0; j < m; ++j) {
        Rcpp::IntegerVector idx_r = group_idx[j];
        const int n_j = idx_r.size();

        // Build arma uvec of 0-based indices for .rows()
        arma::uvec idx_arma(n_j);
        for (int k = 0; k < n_j; ++k)
            idx_arma(k) = static_cast<arma::uword>(idx_r[k]);

        // s_j = colSums of R_tot rows in this group (Armadillo)
        arma::vec s_j = arma::sum(R_tot.rows(idx_arma), 0).t();

        // A_j = Phi_inv + n_j * Sigma_inv  (Eigen)
        Eigen::MatrixXd A_j = Phi_inv_eig + static_cast<double>(n_j) * Sig_inv_eig;

        // V_j = inv(A_j) via LLT  (Eigen)
        Eigen::LLT<Eigen::MatrixXd> llt_A(A_j);
        Eigen::MatrixXd V_j_eig = llt_A.solve(Eigen::MatrixXd::Identity(d, d));

        // Back to Armadillo for M_j = V_j %*% (Sigma_inv %*% s_j)
        arma::mat V_j(V_j_eig.data(), d, d);          // copy from Eigen buffer
        arma::vec M_j = V_j * (Sigma_inv * s_j);

        M_list[j] = Rcpp::wrap(M_j);
        V_list[j] = Rcpp::wrap(V_j);
    }

    return Rcpp::List::create(
        Rcpp::Named("M") = M_list,
        Rcpp::Named("V") = V_list);
}

// Sample updated random intercepts for all groups.
// Implements the conjugate normal-normal posterior draw:
//   psi[j, ] ~ MVN(M_j, V_j)
// using Armadillo's Cholesky-based MVN sampler.
//
// @param group_idx  Rcpp::List of 0-based index vectors (from build_group_idx)
// @param R_tot      N x d residual matrix (W - X %*% beta)
// @param Phi        d x d random-intercept covariance
// @param Sigma_inv  d x d inverse within-group covariance
// @return m x d arma::mat of updated random intercepts
//
// [[Rcpp::export]]
arma::mat update_psi_cpp(
    const Rcpp::List   &group_idx,
    const arma::mat    &R_tot,
    const arma::mat    &Phi,
    const arma::mat    &Sigma_inv)
{
    const int m = static_cast<int>(group_idx.size());
    const int d = static_cast<int>(R_tot.n_cols);

    // Compute Phi_inv once via Eigen LLT
    Eigen::MatrixXd Phi_eig = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Phi));
    Eigen::LLT<Eigen::MatrixXd> llt_phi(Phi_eig);
    Eigen::MatrixXd Phi_inv_eig = llt_phi.solve(Eigen::MatrixXd::Identity(d, d));

    Eigen::MatrixXd Sig_inv_eig = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(Sigma_inv));

    arma::mat psi(m, d);

    for (int j = 0; j < m; ++j) {
        Rcpp::IntegerVector idx_r = group_idx[j];
        const int n_j = idx_r.size();

        arma::uvec idx_arma(n_j);
        for (int k = 0; k < n_j; ++k)
            idx_arma(k) = static_cast<arma::uword>(idx_r[k]);

        // s_j = colSums of residuals for this group
        arma::vec s_j = arma::sum(R_tot.rows(idx_arma), 0).t();

        // A_j = Phi_inv + n_j * Sigma_inv  (Eigen)
        Eigen::MatrixXd A_j = Phi_inv_eig + static_cast<double>(n_j) * Sig_inv_eig;

        // V_j = inv(A_j) via LLT  (Eigen)
        Eigen::LLT<Eigen::MatrixXd> llt_A(A_j);
        Eigen::MatrixXd V_j_eig = llt_A.solve(Eigen::MatrixXd::Identity(d, d));

        // Back to Armadillo for M_j and the MVN draw
        arma::mat V_j(V_j_eig.data(), d, d);
        arma::vec M_j = V_j * (Sigma_inv * s_j);

        // Draw psi_j ~ MVN(M_j, V_j) via Armadillo Cholesky sampler
        // arma::mvnrnd(mu, Sigma) draws one sample from MVN(mu, Sigma)
        psi.row(j) = arma::mvnrnd(M_j, V_j).t();
    }

    return psi;
}
