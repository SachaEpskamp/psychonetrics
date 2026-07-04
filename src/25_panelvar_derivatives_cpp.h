#ifndef DPANELVAR_H
#define DPANELVAR_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_phi_theta_panelvar_cpp(
        const Rcpp::List& prep
);

#endif
