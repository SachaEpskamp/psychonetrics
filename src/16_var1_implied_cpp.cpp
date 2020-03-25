// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List implied_var1_cpp(
    const S4& model,
    bool all = false
){
  int i,j;
  
  S4 sample = model.slot("sample");
  Rcpp::List means = sample.slot("means");
  
  // Extra matrices:
  Rcpp::List extramats = model.slot("extramatrices");
  arma::sp_mat L = extramats["L"];
  
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);
  
  
  // Types:
  Rcpp::List types = model.slot("types");
  std::string zeta = types["zeta"];
  
  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "zeta", zeta, all);
  
  int nGroup = x.length();
  int g;
  
  for (g=0; g<nGroup; g++){
    Rcpp::List grouplist = x[g];
    bool proper = true;
    
    // Model matrices:
    arma::mat beta = grouplist["beta"];
    arma::mat sigma_zeta = grouplist["sigma_zeta"];
    arma::mat exo_cholesky = grouplist["exo_cholesky"];
    // arma::mat mu = grouplist["mu"];
    
    // Matrices I need in every model framework when estimating:
    int n = beta.n_rows;
    arma::mat I2 = eye(n*n, n*n);
    
    // Some stuff needed now:
    arma::mat BetaStar = inv(I2 - kron(beta, beta)); 
    arma::vec sigmaZetaVec = vectorise(sigma_zeta); // FIXME: not needed in fit function
    
    // Implied exogenous covariances:
    arma::mat exoCov = exo_cholesky * exo_cholesky.t();
    
    // Implied stationary distribution (vectorized)
    arma::vec vecSigma = BetaStar * sigmaZetaVec;
    arma::mat Sigma0(n,n); // Not sure if this works...
    int curind = 0;
    for (j=0;j<n;j++){
      for (i=0;i<n;i++){
        Sigma0(i,j) = vecSigma(curind);
        curind++;
      }
    }
    
    // Implied lag-1:
    arma::mat Sigma1 = beta * Sigma0;
    
    // Full implied sigma:
    arma::mat sigma = join_cols(
      join_rows(exoCov, Sigma1.t()),
      join_rows(Sigma1, Sigma0)
    );
    
    
    // Precision:
    arma::mat kappa = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
    
    // Store matrices:
    grouplist["sigma"] = sigma;
    grouplist["kappa"] = kappa;
    
    // Store more:
    if (!all){
      arma::mat L_betaStar = L * BetaStar;
      grouplist["L_betaStar"] = L_betaStar;
      grouplist["IkronBeta"] = kronecker_I_X(beta, n);
      grouplist["sigmaZetaVec"] = sigmaZetaVec;
      grouplist["BetaStar"] = BetaStar;
    } else {
      
      grouplist["PDC"] = computePDC_cpp(beta,grouplist["kappa_zeta"],grouplist["sigma_zeta"]);
      
    }
    
    // Return properness:
    grouplist["proper"] = proper;
    
    
    x[g] = grouplist;
  }
  
  return(x);
}

