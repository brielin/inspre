#include <RcppEigen.h>
#include <iostream>
#include <thread>

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
  return(g>1e-8 ? g : 0.0);
}

// sign function
//
// Returns 1 if x is positive, -1 if negative, 0 if 0.
//
// param x Double.
inline double sign(double x) {
  return (x > 1e-8) ? 1. : ((x < -1e-8) ? -1. : 0.);
}


//' Fit lasso model (one iteration)
//'
//' Minimize .5|| Y - XB ||^2_2 + gamma |B|_1
//'
//' @param X NxP matrix of covariates.
//' @param Y N-vector of response.
//' @param B P-vector of coefficents.
//' @param lambda P-vector L1 regulation coeffs.
//' @param niter Integer, number of lasso iterations to perform.
//' @param fixd Integer, passed from R to indicate the column of the data
//'   during the inspre loop. Set to current column index to fix the diagonal
//'   of the result to 1.
//' @export
// [[Rcpp::export]]
VectorXd lasso(const MatrixXd X, const VectorXd Y, const VectorXd B,
               const double lambda, int niter, int fixd = -1) {
  int P = X.cols();
  VectorXd E = Y - X * B;
  double lambda_B_L1 = 0.0;
  VectorXd B_next = B;
  for(int i=0; i<niter; i++){
    for(int p=0; p<P; p++){
      VectorXd x = X.col(p);
      if(fixd == p+1){
        B_next[p] = 1;
        E += x;
      } else {
        if (std::abs(B_next[p]) > 1e-8){
          E += B_next[p] * x;
        }
        double mp = x.dot(E);
        B_next[p] = sign(mp) * relu(std::abs(mp) - lambda) / x.squaredNorm();
      }
      E -= B_next[p] * x;
      lambda_B_L1 += lambda * std::abs(B_next[p]);
    }
  }
  return(B_next);
}


//' Fit U satisfying VU = I
//'
//' @param WX
//' @param W
//' @param V
//' @param U
//' @param theta
//' @param rho
//' @param gamma
//' @param solve_its
//' @param fixd
//' @param nthreads
//' @export
// [[Rcpp::export]]
MatrixXd fit_U_VU_const(const Map<MatrixXd> WX, const Map<MatrixXd> W, const Map<MatrixXd> V,
                        const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double gamma = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1){
  int D = U.cols();
  MatrixXd rVtV = rho * V.transpose() * V;
  MatrixXd rhs = WX + rho * V.transpose() - V.transpose() * theta;
  MatrixXd U_next = U;

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&](int threadId) {
      for (int d = threadId; d < D; d += nthreads) {
        MatrixXd A = (W.col(d).array() + gamma).matrix().asDiagonal();
        A += rVtV; //diag(W[, d] + 0.2) + rvtv
        U_next.col(d) = lasso(A, rhs.col(d), U.col(d), 0, solve_its, fixd==true ? d + 1: -1);
      }
    }, t);
  }
  for (std::thread& thread : threads) {
    thread.join(); // Wait for all threads to finish
  }
  return(U_next);
}


//' Fit U satisfying UV = I
//'
//' @param WX
//' @param W
//' @param V
//' @param U
//' @param theta
//' @param rho
//' @param gamma
//' @param solve_its
//' @param fixd
//' @param nthreads
//' @export
// [[Rcpp::export]]
MatrixXd fit_U_UV_const(const Map<MatrixXd> WX, const Map<MatrixXd> W, const Map<MatrixXd> V,
                        const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double gamma = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1){
  int D = U.cols();
  MatrixXd rVVt = rho * V * V.transpose();
  MatrixXd rhs = WX.transpose() + rho * V - V * theta.transpose();
  MatrixXd U_next = U;
  // std::cout << rVtV(0,0) << " " << rVtV(0,1) << "\n" << rVtV(1,0) << " " << rVtV(1,1) << " " << std::endl;
  // std::cout << rhs(0,0) << " " << rhs(0,1) << "\n" << rhs(1,0) << " " << rhs(1,1) << " " << std::endl;

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&](int threadId) {
      for (int d = threadId; d < D; d += nthreads) {
        MatrixXd A = (W.row(d).array() + gamma).matrix().asDiagonal();
        A += rVVt; //diag(W[, d] + 0.2) + rvtv
        U_next.row(d) = lasso(A, rhs.col(d), U.row(d), 0, solve_its, fixd==true ? d + 1: -1);
      }
    }, t);
  }
  for (std::thread& thread : threads) {
    thread.join(); // Wait for all threads to finish
  }
  return(U_next);
}


//' Fit V satisfying VU = I
//'
//' @param V
//' @param U
//' @param theta
//' @param rho
//' @param lambda
//' @param solve_its
//' @param fixd
//' @param nthreads
//' @export
// [[Rcpp::export]]
MatrixXd fit_V_VU_const(const Map<MatrixXd> V, const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double lambda = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1){
  int D = V.cols();
  MatrixXd V_next = V;
  MatrixXd glm_Y;
  glm_Y.setIdentity(D, D);
  glm_Y -= (theta.transpose().array()/rho).matrix();

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&](int threadId) {
      for (int d = threadId; d < D; d += nthreads) {
        V_next.row(d) = lasso(U.transpose(), glm_Y.col(d), V.row(d), lambda/rho, solve_its, fixd==true ? d + 1: -1);
      }
    }, t);
  }
  for (std::thread& thread : threads) {
    thread.join(); // Wait for all threads to finish
  }
  return(V_next);
}


//' Fit V satisfying VU = I
//'
//' @param V
//' @param U
//' @param theta
//' @param rho
//' @param lambda
//' @param solve_its
//' @param fixd
//' @param nthreads
//' @export
// [[Rcpp::export]]
MatrixXd fit_V_UV_const(const Map<MatrixXd> V, const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double lambda = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1){
  int D = V.cols();
  MatrixXd V_next = V;
  MatrixXd glm_Y;
  glm_Y.setIdentity(D, D);
  glm_Y -= (theta.array()/rho).matrix();

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&](int threadId) {
      for (int d = threadId; d < D; d += nthreads) {
        V_next.col(d) = lasso(U, glm_Y.col(d), V.col(d), lambda/rho, solve_its, fixd==true ? d + 1: -1);
      }
    }, t);
  }
  for (std::thread& thread : threads) {
    thread.join(); // Wait for all threads to finish
  }
  return(V_next);
}
