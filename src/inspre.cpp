#include <RcppEigen.h>
#include <iostream>

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
//' @param fixd Integer, passed from R to indicate the column of the data
//'   during the inspre loop. Set to current column index to fix the diagonal
//'   of the result to 1.
//' @export
// [[Rcpp::export]]
double lasso(const Map<MatrixXd> X, const Map<VectorXd> Y, Map<VectorXd> B,
             const Map<VectorXd> lambda, int niter, int fixd = -1) {
  int P = X.cols();
  VectorXd E = Y - X * B;
  double lambda_B_L1 = 0.0;
  for(int i=0; i<niter; i++){
    for(int p=0; p<P; p++){
      VectorXd x = X.col(p);
      if(fixd == p+1){
        B[p] = 1;
        E += x;
      } else {
        if (B[p] != 0){
          E += B[p] * x;
        }
        double mp = x.dot(E);
        B[p] = sign(mp) * relu(std::abs(mp) - lambda[p]) / x.squaredNorm();
      }
      E -= B[p] * x;
      lambda_B_L1 += lambda[p] * std::abs(B[p]);
    }
  }
  return(.5*E.squaredNorm() + lambda_B_L1);
}


//' Fit lasso model for matrix of responses (one iteration)
//'
//' Minimize .5 || Y - XB ||^2_2 + gamma |B|_1
//'
//' This is ~4x faster than the above without multithreading but
//' the above is still faster if you use 8+ threads in my testing.
//' Not being used for now but leaving it in for reference.
//'
//' @param X NxP matrix of covariates.
//' @param Y NxD matrix of response.
//' @param B PxD matrix of coefficents (to be updated).
//' @param lambda PxD matrix. L1 regulation coeffs.
//' @param niter Integer, number of lasso iterations to perform.
//' @export
// [[Rcpp::export]]
MatrixXd matrix_lasso(const Map<MatrixXd> X, const Map<MatrixXd> Y, const Map<MatrixXd> B,
             const Map<ArrayXXd> lambda, int niter) {
  int P = X.cols();
  MatrixXd E = Y - X * B;
  ArrayXXd B_next = B.array();
  double lambda_B_L1 = 0.0;
  VectorXd X_col_norms = X.colwise().squaredNorm();
  for(int i=0; i<niter; i++){
    for(int p=0; p<P; p++){
      VectorXd x = X.col(p);
      E += x * B_next.matrix().row(p);
      ArrayXd mp = x.transpose() * E;
      ArrayXd B_p = mp.abs();
      B_p -= lambda.row(p);
      B_p = B_p.unaryExpr([](double x){return x>0.0 ? x : 0.0;});
      B_p = mp.sign() * B_p / X_col_norms[p];
      B_next.row(p) = B_p;
      E -= x * B_next.matrix().row(p);
    }
    // TODO(brielin): This doesn't work and should be removed eventually but
    //   I'm leaving it in for now in case I need to reference it later.
    // E += X * B_next.matrix();
    // ArrayXXd mp = X.transpose() * E;
    // B_next = mp.abs();
    // B_next -= lambda;
    // B_next = B_next.unaryExpr([](double x){return x>0.0 ? x : 0.0;});
    // B_next *= mp.sign();
    // B_next.colwise() /= X_col_norms.array();
    // E -= X*B_next.matrix();
    // lambda_B_L1 = lambda.cwiseProduct(B_next.cwiseAbs()).sum();
  }
  return(B_next);
}
