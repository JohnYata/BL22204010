// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// simulate_bandit
NumericMatrix simulate_bandit(int N, int T, double sigma, double delta, NumericVector mu);
RcppExport SEXP _BL22204010_simulate_bandit(SEXP NSEXP, SEXP TSEXP, SEXP sigmaSEXP, SEXP deltaSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_bandit(N, T, sigma, delta, mu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BL22204010_simulate_bandit", (DL_FUNC) &_BL22204010_simulate_bandit, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_BL22204010(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
