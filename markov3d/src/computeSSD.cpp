#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat computeSSD(arma::mat m, double epsilon, int max_iter) {
  //' Compute steady-state probability from transition matrix by iteration
  //'
  //' @param m Arma::mat, the transition matrix
  //' @param epsilon Double, the tolerence of the computation
  //' @param max_iter Int, the max number of iteration
  //' @return Arma::mat, steady-state probability
  int vec_length = m.n_cols;
  double min_diff = 100.0;
  double this_diff;
  arma::mat miu0(1, vec_length);
  arma::mat miu1(1, vec_length);
  arma::mat miu_optim(1, vec_length);

  miu0.fill(1.0 / vec_length);
  for (int i = 0; i < max_iter; i++) {
    miu1 = miu0 * m;
    this_diff = abs((miu1 - miu0)).max();
    if (this_diff < epsilon) {
      return miu1;
    } else if (this_diff < min_diff) {
      min_diff = this_diff;
      miu_optim = miu1;
    }
    miu0 = miu1;
  }
  return miu_optim;
}
