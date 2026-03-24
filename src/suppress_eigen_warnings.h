// Suppress [-Wignored-attributes] emitted by Eigen 3.3.x SSE/SIMD headers under GCC 13+.
// Injected via -include in Makevars.win so it applies to all translation units,
// including the auto-generated RcppExports.cpp which cannot be edited directly.
// Remove once RcppEigen ships Eigen >= 3.4 on CRAN.
#pragma GCC diagnostic ignored "-Wignored-attributes"
