#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix qqconf(const int n, int conf_points, double conf_alpha) {

  NumericMatrix mpts(conf_points * 2, 2);

  for (int i = 0; i < conf_points; ++i) {
    double expected_log10p = -std::log10(((double)i + 0.5) / n);
    mpts(i, 0) = expected_log10p;
    mpts(i, 1) = -std::log10(R::qbeta(1 - conf_alpha / 2.0, i + 1, n - (i + 1), 1.0, 0.0));
    mpts(conf_points * 2 - 1 - i, 0) = expected_log10p;
    mpts(conf_points * 2 - 1 - i, 1) = -std::log10(R::qbeta(conf_alpha / 2.0, i + 1, n - (i + 1), 1.0, 0.0));
  }

  return mpts;
}
