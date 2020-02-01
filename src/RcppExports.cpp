// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "inspre_types.hpp"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// lasso_one_iteration
double lasso_one_iteration(const Map<MatrixXd> X, const Map<VectorXd> Y, Map<VectorXd> B, const Map<VectorXd> lambda);
RcppExport SEXP _inspre_lasso_one_iteration(SEXP XSEXP, SEXP YSEXP, SEXP BSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Map<MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Map<VectorXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< const Map<VectorXd> >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_one_iteration(X, Y, B, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_inspre_lasso_one_iteration", (DL_FUNC) &_inspre_lasso_one_iteration, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_inspre(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
