// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rgammatr
double rgammatr(double a, double b, arma::vec range, double n);
RcppExport SEXP _fBN_rgammatr(SEXP aSEXP, SEXP bSEXP, SEXP rangeSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rgammatr(a, b, range, n));
    return rcpp_result_gen;
END_RCPP
}
// checkDAG
bool checkDAG(arma::mat G);
RcppExport SEXP _fBN_checkDAG(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(checkDAG(G));
    return rcpp_result_gen;
END_RCPP
}
// get_zcoef
Rcpp::List get_zcoef(arma::mat xobs, arma::mat zobs, arma::mat scoef, arma::mat sbasis, arma::mat fscore, arma::mat eta, arma::vec svar, arma::uvec tindex, arma::uvec cindex, arma::uvec findex);
RcppExport SEXP _fBN_get_zcoef(SEXP xobsSEXP, SEXP zobsSEXP, SEXP scoefSEXP, SEXP sbasisSEXP, SEXP fscoreSEXP, SEXP etaSEXP, SEXP svarSEXP, SEXP tindexSEXP, SEXP cindexSEXP, SEXP findexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zobs(zobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scoef(scoefSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sbasis(sbasisSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fscore(fscoreSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tindex(tindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type findex(findexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_zcoef(xobs, zobs, scoef, sbasis, fscore, eta, svar, tindex, cindex, findex));
    return rcpp_result_gen;
END_RCPP
}
// get_scoef
Rcpp::List get_scoef(arma::mat xobs, arma::mat zobs, arma::mat zcoef, arma::mat scoef, arma::mat sbasis, arma::mat penat, arma::mat fscore, arma::vec svar, arma::vec lambda, arma::uvec tindex, arma::uvec cindex, arma::uvec findex, bool orderLambdas);
RcppExport SEXP _fBN_get_scoef(SEXP xobsSEXP, SEXP zobsSEXP, SEXP zcoefSEXP, SEXP scoefSEXP, SEXP sbasisSEXP, SEXP penatSEXP, SEXP fscoreSEXP, SEXP svarSEXP, SEXP lambdaSEXP, SEXP tindexSEXP, SEXP cindexSEXP, SEXP findexSEXP, SEXP orderLambdasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zobs(zobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zcoef(zcoefSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scoef(scoefSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sbasis(sbasisSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type penat(penatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fscore(fscoreSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tindex(tindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type findex(findexSEXP);
    Rcpp::traits::input_parameter< bool >::type orderLambdas(orderLambdasSEXP);
    rcpp_result_gen = Rcpp::wrap(get_scoef(xobs, zobs, zcoef, scoef, sbasis, penat, fscore, svar, lambda, tindex, cindex, findex, orderLambdas));
    return rcpp_result_gen;
END_RCPP
}
// get_fscore
arma::mat get_fscore(arma::mat xobs, arma::mat zobs, arma::mat zcoef, arma::mat scoef, arma::mat sbasis, arma::mat fvar, arma::mat fcoef, arma::vec svar, arma::uvec tindex, arma::uvec cindex, arma::uvec findex);
RcppExport SEXP _fBN_get_fscore(SEXP xobsSEXP, SEXP zobsSEXP, SEXP zcoefSEXP, SEXP scoefSEXP, SEXP sbasisSEXP, SEXP fvarSEXP, SEXP fcoefSEXP, SEXP svarSEXP, SEXP tindexSEXP, SEXP cindexSEXP, SEXP findexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zobs(zobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zcoef(zcoefSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scoef(scoefSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sbasis(sbasisSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fvar(fvarSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fcoef(fcoefSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type svar(svarSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tindex(tindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type findex(findexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_fscore(xobs, zobs, zcoef, scoef, sbasis, fvar, fcoef, svar, tindex, cindex, findex));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _fBN_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _fBN_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _fBN_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _fBN_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fBN_rgammatr", (DL_FUNC) &_fBN_rgammatr, 4},
    {"_fBN_checkDAG", (DL_FUNC) &_fBN_checkDAG, 1},
    {"_fBN_get_zcoef", (DL_FUNC) &_fBN_get_zcoef, 10},
    {"_fBN_get_scoef", (DL_FUNC) &_fBN_get_scoef, 13},
    {"_fBN_get_fscore", (DL_FUNC) &_fBN_get_fscore, 11},
    {"_fBN_rcpparma_hello_world", (DL_FUNC) &_fBN_rcpparma_hello_world, 0},
    {"_fBN_rcpparma_outerproduct", (DL_FUNC) &_fBN_rcpparma_outerproduct, 1},
    {"_fBN_rcpparma_innerproduct", (DL_FUNC) &_fBN_rcpparma_innerproduct, 1},
    {"_fBN_rcpparma_bothproducts", (DL_FUNC) &_fBN_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fBN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
