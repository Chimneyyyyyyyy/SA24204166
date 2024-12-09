// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ARI
double ARI(int k, int n, IntegerVector true_label, IntegerVector label);
RcppExport SEXP _SA24204166_ARI(SEXP kSEXP, SEXP nSEXP, SEXP true_labelSEXP, SEXP labelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type true_label(true_labelSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type label(labelSEXP);
    rcpp_result_gen = Rcpp::wrap(ARI(k, n, true_label, label));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204166_ARI", (DL_FUNC) &_SA24204166_ARI, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204166(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
