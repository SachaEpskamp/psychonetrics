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
Rcpp::List implied_tsdlvm1_cpp(
    const S4& model,
    bool all = false
){
  
  S4 sample = model.slot("sample");
  Rcpp::List means = sample.slot("means");
  
  
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);
  
  
  // Types:
  Rcpp::List types = model.slot("types");
  
  std::string zeta = types["zeta"];
  std::string epsilon = types["epsilon"];
  
  
  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "zeta", zeta, all);
  x = impliedcovstructures_cpp(x, "epsilon", epsilon, all);
  
  int nGroup = x.length();
  int g;
  
  // General stuff:
  Rcpp::List extramats = model.slot("extramatrices");
  
  
  
  
  
  for (g=0; g<nGroup; g++){
    bool proper = true;
    
    Rcpp::List grouplist = x[g];
    
    // Model matrices:
    arma::mat exo_means = grouplist["exo_means"];
    arma::mat exo_cholesky = grouplist["exo_cholesky"];
    arma::mat nu = grouplist["nu"];
    arma::mat mu_eta = grouplist["mu_eta"];
    arma::mat lambda = grouplist["lambda"];
    arma::mat beta = grouplist["beta"];
    arma::mat sigma_zeta = grouplist["sigma_zeta"];
    arma::mat sigma_epsilon = grouplist["sigma_epsilon"];
    
    
    // Matrices I need in every model framework when estimating:
    // Beta star:
    int n_lat = beta.n_rows;
    arma::mat I2 = eye(n_lat*n_lat, n_lat*n_lat);
    
    // Some stuff needed now:
    arma::mat BetaStar = inv(I2 - kron(beta, beta));
    
    // Implied mean vector:
    arma::vec impMu = nu + lambda * mu_eta;
      
     arma::vec fullMu = join_cols(exo_means,impMu);
      
      // Exogenous cov part:
      arma::mat exoCov = exo_cholesky * exo_cholesky.t();
      
      // Latent lag-0:
      // int nLatent = lambda.n_cols;
      arma::mat Sigma_eta_0 = matrixform(BetaStar * vectorise(sigma_zeta));
      
      // Observed stationary:
      arma::mat Sigma_y_0 =  lambda *  Sigma_eta_0 * lambda.t() + sigma_epsilon;
      
      // Lag 1 part:
      arma::mat Sigma_eta_1 = beta * Sigma_eta_0;
      
      // Lag 1 observed:
      arma::mat Sigma_y_1 =  lambda *  Sigma_eta_1 * lambda.t();
      
      
      // Subset and add to the list:
      grouplist["mu"] = fullMu;
      
      // Full implied sigma:
      arma::mat sigma = join_cols(
          join_rows(exoCov,Sigma_y_1.t()),
          join_rows(Sigma_y_1,Sigma_y_0)
      );
       grouplist["sigma"] = sigma;
      
      // Precision:
      grouplist["kappa"] = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);

      
      // Extra matrices needed in optimization:
      if (!all){
        grouplist["BetaStar"] = BetaStar;
        grouplist["Sigma_eta_0"] = Sigma_eta_0;
        grouplist["Sigma_eta_1"] = Sigma_eta_1;
        grouplist["IkronBeta"] = kronecker_I_X(beta, n_lat);
        grouplist["lamWkronlamW"] = kron(lambda, lambda);
      } else {
        // if (!is.null(grouplist["kappa_zeta_within)){
        arma::mat kappa_zeta;
        
        if (grouplist.containsElementNamed("kappa_zeta")){
          kappa_zeta = as<arma::mat>(grouplist["kappa_zeta"]);
        } else {
          kappa_zeta = solve_symmetric_cpp_matrixonly_withcheck(sigma_zeta, proper);
        }
        
        grouplist["PDC"] = computePDC_cpp(grouplist["beta"], kappa_zeta, sigma_zeta);
      }
      
      // Return properness:
      grouplist["proper"] = proper;
      
      x[g] = grouplist;
  }
  
  return(x);
}

