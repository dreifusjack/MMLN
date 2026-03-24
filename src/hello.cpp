// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#pragma GCC diagnostic pop

using namespace Rcpp;

// DEV NOTE:
// {[[Rcpp::export]]} Makes this function callable from R.
// Remove this annotation for any internal C++ helpers that R doesn't need to call directly.

// [[Rcpp::export]]
void hello_world()
{
  Rcpp::Rcout << "C++ is working w/ Rcpp" << std::endl;
}