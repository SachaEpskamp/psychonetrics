// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "02_algebrahelpers_RcppHelpers.h"

// mat test(
//   mat X,
//   mat Y,
//   uvec ind
// ){
//   X(ind,ind) += Y;
//   return X;
// }



// Inner function
// [[Rcpp::export]]
double fimlEstimator_Gauss_group_cpp_inner(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List dat,
    double epsilon,
    double n) {
  
  // double logdet = 0;
  double n_part;
  
  // Observed indices:
  vec obs = dat["obs"];
  uvec inds = find(obs == true);
  
  // sample size:
  n_part = dat["n"];
  
  // Subset matrices:
  arma::mat sigma_p = sigma(inds,inds);
  // arma::mat kappa_p(sigma_p);
  arma::vec mu_p = mu(inds);
  
  // Observed values:
  arma::mat S = dat["S"];
  arma::vec means = dat["means"];
  
  // log det:
  // arma::vec ev = arma::eig_sym(sigma_p);
  // bool ispos = ev[0] > epsilon;

  // bool ispos = true;
  // for (int j = 0; j < ev.size(); j++){
  //   if (ev[j] < sqrt(epsilon)){
  //     ispos = false;
  //     break;
  //   }
  // }
  // double logepsilon = log(epsilon);
  // 
  // if (ispos){
  //   kappa_p = inv(sigma_p);
  //   // logdet = log(det(kappa_p));
  //   // logdet =  real(log_det(kappa_p));
  //   logdet =  log(det(kappa_p));
  //   if (logdet < logepsilon){
  //     logdet = logepsilon;
  //   }
  // } else {
  //   kappa_p = pinv(sigma_p);
  //   logdet = logepsilon;
  // }
  List invres = solve_symmetric_cpp(sigma_p, true, epsilon);
  arma::mat kappa_p = invres["inv"];
  double logdet = invres["logdet"];
  
  // Rf_PrintValue(wrap(log_det(kappa_p)));
  
  // Likelihood:
  double result = n_part * (trace(S * kappa_p) + 
                            dot((means - mu_p).t(), kappa_p * (means - mu_p)) - 
                            logdet);

  
  // Return
  return result;
}


// Outer normal function:
// [[Rcpp::export]] 
double fimlEstimator_Gauss_group_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon,
    double n) {
  // Rf_PrintValue(wrap("USED"));
  
  double result = 0;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
    
    result += fimlEstimator_Gauss_group_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
      epsilon,
      n);
  }
  
  
  // Return
  return (1/n) * result;
}



// Outer function PER GROUP:
// [[Rcpp::export]] 
double fimlEstimator_Gauss_group_cpp_fullFIML(
    Rcpp::List sigma, 
    Rcpp::List kappa,
    Rcpp::List mu,
    Rcpp::List fimldata,
    double epsilon,
    double n) {
  // Rf_PrintValue(wrap("USED"));
  
  double result = 0;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
    
    
    result += fimlEstimator_Gauss_group_cpp_inner(
      sigma[i], 
      kappa[i],
      mu[i],
      fimldata[i],
              epsilon,
              n);
  }
  
  
  // Return
  return (1/n) * result;
}




// fMain function:
// [[Rcpp::export]]
double fimlestimator_Gauss_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  bool fullFIML = prep["fullFIML"];
  
  // Result:
  double fit = 0;
  
  for (int i=0; i<nGroup;i++){
    
    Rcpp::List grouplist = groupmodels[i];
    
    if (fullFIML){
      
      fit += (nPerGroup(i) / nTotal) * fimlEstimator_Gauss_group_cpp_fullFIML(grouplist["sigma"],grouplist["kappa"],grouplist["mu"],grouplist["fimldata"], 1.490116e-08,grouplist["fulln"]);
      
    } else {
      
      fit += (nPerGroup(i) / nTotal) * fimlEstimator_Gauss_group_cpp(grouplist["sigma"],grouplist["kappa"],grouplist["mu"],grouplist["fimldata"], 1.490116e-08, grouplist["fulln"]);
      
    }
    
  }
  
  return(fit);
}




