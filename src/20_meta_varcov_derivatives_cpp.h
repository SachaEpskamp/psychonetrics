#ifndef DMETAVARCOV_H
#define DMETAVARCOV_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_phi_theta_meta_varcov_cpp(
        const Rcpp::List& prep
);


#endif
