#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace Eigen;

// relu function
//
// Returns g if g > 0, else 0.
//
// param g Double.
inline double relu(double g) {
  return(g>0.0 ? g : 0.0);
}

// sign function
//
// Returns 1 if x is positive, -1 if negative, 0 if 0.
//
// param x Double.
inline double sign(double x) {
  return (x > 0.) ? 1. : ((x < 0.) ? -1. : 0.);
}

//' Fit lasso model (one iteration)
//'
//' Minimize .5 || Y - XB ||^2_2 + gamma |B|_1
//'
//' @param X NxP matrix of covariates.
//' @param Y N-vector of response.
//' @param B P-vector of coefficents (to be updated).
//' @param lambda P-vector L1 regulation coeffs.
//' @param niter Integer, number of lasso iterations to perform.
//' @export
// [[Rcpp::export]]
double lasso(const Map<MatrixXd> X,
                           const Map<VectorXd> Y,
                           Map<VectorXd> B,
                           const Map<VectorXd> lambda,
                           int niter) {
  int P = X.cols();
  VectorXd E = Y - X * B;
  double lambda_B_L1 = 0.0;
  for(int i=0; i<niter; i++){
    for(int p=0; p<P; p++){
      VectorXd x = X.col(p);
      if (B[p] != 0)
        E += B[p] * x;
      double mp = x.dot(E);
      B[p] = sign(mp) * relu(std::abs(mp) - lambda[p]) / x.squaredNorm();
      E -= B[p] * x;
      lambda_B_L1 += lambda[p] * std::abs(B[p]);
    }
  }
  return(.5*E.squaredNorm() + lambda_B_L1);
}
