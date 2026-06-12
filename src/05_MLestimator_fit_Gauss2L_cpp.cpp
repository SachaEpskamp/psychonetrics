// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// C++ twin of the two-level (random intercept) sufficient-statistics ML fit
// function for multivariate normal data with cluster structure (ml_lvm,
// estimator = "ML", distribution "TwoLevelGaussian"). Mirrors
// R/05_MLestimator_fit_Gauss2L.R exactly; see that file for the published
// likelihood decomposition (McDonald & Goldstein, 1989; Muthen, 1990).

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "05_MLestimator_fit_Gauss2L_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// -2 l* for one group, WITHOUT the N*p*log(2*pi) constant (twin of
// minustwo_logl_Gauss2L_noconstant):
double minustwo_logl_Gauss2L_noconstant_cpp(
        const arma::vec& mu,
        const arma::mat& SW,
        const arma::mat& SB,
        const Rcpp::List& twolevel
){
    double N = Rcpp::as<double>(twolevel["N"]);
    double J = Rcpp::as<double>(twolevel["J"]);
    arma::mat S_PW = twolevel["S_PW"];
    Rcpp::List sizes = twolevel["sizes"];
    arma::vec nj_vec = Rcpp::as<arma::vec>(sizes["nj"]);
    arma::vec m_vec = Rcpp::as<arma::vec>(sizes["m"]);
    Rcpp::List mean_d = twolevel["mean_d"];
    Rcpp::List cov_d = twolevel["cov_d"];

    // Non positive-definite within covariance: return a large penalty,
    // mirroring the R twin:
    if (sympd_cpp(SW) == false){
        return(1e20);
    }

    // Pooled-within part (zero weight when all clusters have size 1):
    double res = 0;
    double logdet_val, logdet_sign;
    if (N > J){
        arma::mat iSW = solve_symmetric_cpp_matrixonly(SW);
        log_det(logdet_val, logdet_sign, SW);
        res = (N - J) * (logdet_val + accu(iSW % S_PW));
    }

    // Between part, per distinct cluster size:
    int nSizes = nj_vec.n_elem;
    for (int s = 0; s < nSizes; s++){
        double nj = nj_vec(s);
        double m = m_vec(s);
        arma::mat Sj = SW + nj * SB;
        if (sympd_cpp(Sj) == false){
            return(1e20);
        }
        arma::mat iSj = solve_symmetric_cpp_matrixonly(Sj);
        arma::vec yc = Rcpp::as<arma::vec>(mean_d[s]) - mu;
        arma::mat A = Rcpp::as<arma::mat>(cov_d[s]) + yc * yc.t();
        log_det(logdet_val, logdet_sign, Sj);
        res += m * (logdet_val + nj * accu(iSj % A));
    }

    return(res);
}


// GROUP FIT FUNCTION: (-2 l*) / J (no 2*pi constant) //
// [[Rcpp::export]]
double maxLikEstimator_Gauss2L_group_cpp(
        const Rcpp::List& grouplist
){
    arma::vec mu = grouplist["mu"];
    arma::mat SW = grouplist["sigma_within"];
    arma::mat SB = grouplist["sigma_between"];
    Rcpp::List twolevel = grouplist["twolevel"];

    double J = Rcpp::as<double>(twolevel["J"]);

    return(minustwo_logl_Gauss2L_noconstant_cpp(mu, SW, SB, twolevel) / J);
}


// full fit function
// [[Rcpp::export]]
double maxLikEstimator_Gauss2L_cpp(
        const Rcpp::List& prep
){
    Rcpp::List groupmodels = prep["groupModels"];
    int nGroup = groupmodels.length();
    arma::vec nPerGroup = prep["nPerGroup"];
    double nTotal = prep["nTotal"];

    // Result:
    double fit = 0;

    for (int i = 0; i < nGroup; i++){
        fit += (nPerGroup(i) / nTotal) * maxLikEstimator_Gauss2L_group_cpp(groupmodels[i]);
    }

    return(fit);
}
