// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikelihood_gaussian_subgroup_fiml_cpp
double logLikelihood_gaussian_subgroup_fiml_cpp(arma::mat sigma, arma::mat kappa, arma::vec mu, Rcpp::List fimldata, double epsilon);
RcppExport SEXP _psychonetrics_logLikelihood_gaussian_subgroup_fiml_cpp(SEXP sigmaSEXP, SEXP kappaSEXP, SEXP muSEXP, SEXP fimldataSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fimldata(fimldataSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihood_gaussian_subgroup_fiml_cpp(sigma, kappa, mu, fimldata, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// DWLS_wmat
arma::sp_mat DWLS_wmat(arma::mat data, arma::vec means, const int ncase, const int nvar);
RcppExport SEXP _psychonetrics_DWLS_wmat(SEXP dataSEXP, SEXP meansSEXP, SEXP ncaseSEXP, SEXP nvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type means(meansSEXP);
    Rcpp::traits::input_parameter< const int >::type ncase(ncaseSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    rcpp_result_gen = Rcpp::wrap(DWLS_wmat(data, means, ncase, nvar));
    return rcpp_result_gen;
END_RCPP
}
// WLS_wmat
arma::mat WLS_wmat(arma::mat data, arma::vec means, const int ncase, const int nvar);
RcppExport SEXP _psychonetrics_WLS_wmat(SEXP dataSEXP, SEXP meansSEXP, SEXP ncaseSEXP, SEXP nvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type means(meansSEXP);
    Rcpp::traits::input_parameter< const int >::type ncase(ncaseSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    rcpp_result_gen = Rcpp::wrap(WLS_wmat(data, means, ncase, nvar));
    return rcpp_result_gen;
END_RCPP
}
// expected_hessian_fiml_Gaussian_group_cpp
arma::mat expected_hessian_fiml_Gaussian_group_cpp(arma::mat sigma, arma::mat kappa, arma::vec mu, Rcpp::List fimldata, double epsilon);
RcppExport SEXP _psychonetrics_expected_hessian_fiml_Gaussian_group_cpp(SEXP sigmaSEXP, SEXP kappaSEXP, SEXP muSEXP, SEXP fimldataSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fimldata(fimldataSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_hessian_fiml_Gaussian_group_cpp(sigma, kappa, mu, fimldata, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// fimlEstimator_Gauss_group_cpp
double fimlEstimator_Gauss_group_cpp(arma::mat sigma, arma::mat kappa, arma::vec mu, Rcpp::List fimldata, double epsilon, double n);
RcppExport SEXP _psychonetrics_fimlEstimator_Gauss_group_cpp(SEXP sigmaSEXP, SEXP kappaSEXP, SEXP muSEXP, SEXP fimldataSEXP, SEXP epsilonSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fimldata(fimldataSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(fimlEstimator_Gauss_group_cpp(sigma, kappa, mu, fimldata, epsilon, n));
    return rcpp_result_gen;
END_RCPP
}
// jacobian_fiml_gaussian_subgroup_sigma_cpp
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp(arma::mat sigma, arma::mat kappa, arma::vec mu, Rcpp::List fimldata, double epsilon);
RcppExport SEXP _psychonetrics_jacobian_fiml_gaussian_subgroup_sigma_cpp(SEXP sigmaSEXP, SEXP kappaSEXP, SEXP muSEXP, SEXP fimldataSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fimldata(fimldataSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(jacobian_fiml_gaussian_subgroup_sigma_cpp(sigma, kappa, mu, fimldata, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// covPrepare_cpp
List covPrepare_cpp(NumericMatrix data, LogicalVector isOrdered);
RcppExport SEXP _psychonetrics_covPrepare_cpp(SEXP dataSEXP, SEXP isOrderedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type isOrdered(isOrderedSEXP);
    rcpp_result_gen = Rcpp::wrap(covPrepare_cpp(data, isOrdered));
    return rcpp_result_gen;
END_RCPP
}
// computeMean
double computeMean(NumericVector y);
RcppExport SEXP _psychonetrics_computeMean(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(computeMean(y));
    return rcpp_result_gen;
END_RCPP
}
// computeThresholds
NumericVector computeThresholds(IntegerVector y);
RcppExport SEXP _psychonetrics_computeThresholds(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(computeThresholds(y));
    return rcpp_result_gen;
END_RCPP
}
// pearsonCov
double pearsonCov(Rcpp::NumericVector y1, Rcpp::NumericVector y2, double mean1, double mean2, bool unbiased);
RcppExport SEXP _psychonetrics_pearsonCov(SEXP y1SEXP, SEXP y2SEXP, SEXP mean1SEXP, SEXP mean2SEXP, SEXP unbiasedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< double >::type mean1(mean1SEXP);
    Rcpp::traits::input_parameter< double >::type mean2(mean2SEXP);
    Rcpp::traits::input_parameter< bool >::type unbiased(unbiasedSEXP);
    rcpp_result_gen = Rcpp::wrap(pearsonCov(y1, y2, mean1, mean2, unbiased));
    return rcpp_result_gen;
END_RCPP
}
// toOrdinal
IntegerVector toOrdinal(NumericVector var);
RcppExport SEXP _psychonetrics_toOrdinal(SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(toOrdinal(var));
    return rcpp_result_gen;
END_RCPP
}
// cpp_table
IntegerMatrix cpp_table(IntegerVector y1, IntegerVector y2);
RcppExport SEXP _psychonetrics_cpp_table(SEXP y1SEXP, SEXP y2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y2(y2SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_table(y1, y2));
    return rcpp_result_gen;
END_RCPP
}
// polychoric_fit_summary
double polychoric_fit_summary(double rho, IntegerMatrix tab, NumericVector t1, NumericVector t2);
RcppExport SEXP _psychonetrics_polychoric_fit_summary(SEXP rhoSEXP, SEXP tabSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(polychoric_fit_summary(rho, tab, t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// binormal_density
double binormal_density(double x1, double x2, double rho, double sigma1, double sigma2, double mu1, double mu2);
RcppExport SEXP _psychonetrics_binormal_density(SEXP x1SEXP, SEXP x2SEXP, SEXP rhoSEXP, SEXP sigma1SEXP, SEXP sigma2SEXP, SEXP mu1SEXP, SEXP mu2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< double >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< double >::type mu2(mu2SEXP);
    rcpp_result_gen = Rcpp::wrap(binormal_density(x1, x2, rho, sigma1, sigma2, mu1, mu2));
    return rcpp_result_gen;
END_RCPP
}
// polychoric_grad_summary
double polychoric_grad_summary(double rho, IntegerMatrix tab, NumericVector t1, NumericVector t2);
RcppExport SEXP _psychonetrics_polychoric_grad_summary(SEXP rhoSEXP, SEXP tabSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(polychoric_grad_summary(rho, tab, t1, t2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_psychonetrics_logLikelihood_gaussian_subgroup_fiml_cpp", (DL_FUNC) &_psychonetrics_logLikelihood_gaussian_subgroup_fiml_cpp, 5},
    {"_psychonetrics_DWLS_wmat", (DL_FUNC) &_psychonetrics_DWLS_wmat, 4},
    {"_psychonetrics_WLS_wmat", (DL_FUNC) &_psychonetrics_WLS_wmat, 4},
    {"_psychonetrics_expected_hessian_fiml_Gaussian_group_cpp", (DL_FUNC) &_psychonetrics_expected_hessian_fiml_Gaussian_group_cpp, 5},
    {"_psychonetrics_fimlEstimator_Gauss_group_cpp", (DL_FUNC) &_psychonetrics_fimlEstimator_Gauss_group_cpp, 6},
    {"_psychonetrics_jacobian_fiml_gaussian_subgroup_sigma_cpp", (DL_FUNC) &_psychonetrics_jacobian_fiml_gaussian_subgroup_sigma_cpp, 5},
    {"_psychonetrics_covPrepare_cpp", (DL_FUNC) &_psychonetrics_covPrepare_cpp, 2},
    {"_psychonetrics_computeMean", (DL_FUNC) &_psychonetrics_computeMean, 1},
    {"_psychonetrics_computeThresholds", (DL_FUNC) &_psychonetrics_computeThresholds, 1},
    {"_psychonetrics_pearsonCov", (DL_FUNC) &_psychonetrics_pearsonCov, 5},
    {"_psychonetrics_toOrdinal", (DL_FUNC) &_psychonetrics_toOrdinal, 1},
    {"_psychonetrics_cpp_table", (DL_FUNC) &_psychonetrics_cpp_table, 2},
    {"_psychonetrics_polychoric_fit_summary", (DL_FUNC) &_psychonetrics_polychoric_fit_summary, 4},
    {"_psychonetrics_binormal_density", (DL_FUNC) &_psychonetrics_binormal_density, 7},
    {"_psychonetrics_polychoric_grad_summary", (DL_FUNC) &_psychonetrics_polychoric_grad_summary, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_psychonetrics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
