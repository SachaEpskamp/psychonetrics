// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// mat test(
//   mat X,
//   mat Y,
//   uvec ind
// ){
//   X(ind,ind) += Y;
//   return X;
// }


// innner function
// [[Rcpp::export]]
double logLikelihood_gaussian_subgroup_fiml_cpp_inner(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    Rcpp::List dat,
    double epsilon) {
  // Rf_PrintValue(wrap("USED"));
  
  double logdet = 0;
  double n_part;
  double log2pi = log(2*M_PI);
  
  // Loop over groups
  
  // Obtain list:
  // Rcpp::List dat = fimldata[i];
  
  // Observed indices:
  vec obs = dat["obs"];
  
  uvec inds = find(obs == true);
  
  // sample size:
  n_part = dat["n"];
  
  // Subset matrices:
  arma::mat sigma_p = sigma(inds,inds);
  arma::mat kappa_p(sigma_p);
  arma::vec mu_p = mu(inds);
  
  int nvar = mu_p.size();
  
  // Observed values:
  arma::mat S = dat["S"];
  arma::vec means = dat["means"];
  
  // log det:
  arma::vec ev = arma::eig_sym(sigma_p);
  bool ispos = true;
  int nEV = ev.size();
  for (int j = 0; j < nEV; j++){
    if (ev[j] < sqrt(epsilon)){
      ispos = false;
      break;
    }
  }
  if (ispos){
    kappa_p = inv(sigma_p);
    logdet = log(det(kappa_p));
    double logepsilon = log(epsilon);
    if (logdet < logepsilon){
      logdet = logepsilon;
    }
    if (logdet == R_PosInf){
      logdet = real(log_det(kappa_p));
    }
    
  } else {
    kappa_p = pinv(sigma_p);
    logdet = log(epsilon);
  }
  
  // Likelihood:
  double result = n_part * (
    logdet - nvar * log2pi - trace(S * kappa_p) -
      dot((means - mu_p).t(), kappa_p * (means - mu_p))
  );
  
  // nvar <- ncol(kappa)
  // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  
  // Return
  return result;
}


// Outer function
// [[Rcpp::export]]
double logLikelihood_gaussian_subgroup_fiml_cpp(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  double result = 0;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
    
    
    // Likelihood:
    result += logLikelihood_gaussian_subgroup_fiml_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
              epsilon) ;
    
    // nvar <- ncol(kappa)
    // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  }
  
  
  // Return
  return result;
}


// Outer function PER GORUP
// sigma, kappa and mu may each be supplied either as a per-row list (one
// element per row of the data) or as a single matrix/vector that applies to
// every row (e.g. the varcov family). Both cases are handled here without
// materializing per-row copies.
// [[Rcpp::export]]
double logLikelihood_gaussian_subgroup_fiml_cpp_fullFIML(
    SEXP sigma,
    SEXP kappa,
    SEXP mu,
    Rcpp::List fimldata,
    double epsilon) {


  double result = 0;

  // Detect per-row lists vs single matrices/vectors:
  bool sigmaIsList = Rf_isVectorList(sigma);
  bool kappaIsList = Rf_isVectorList(kappa);
  bool muIsList = Rf_isVectorList(mu);

  Rcpp::List sigmaList, kappaList, muList;
  arma::mat curSigma, curKappa;
  arma::vec curMu;

  if (sigmaIsList) sigmaList = sigma; else curSigma = Rcpp::as<arma::mat>(sigma);
  if (kappaIsList) kappaList = kappa; else curKappa = Rcpp::as<arma::mat>(kappa);
  if (muIsList) muList = mu; else curMu = Rcpp::as<arma::vec>(mu);

  // Loop over rows
  for (int i = 0; i < fimldata.size(); i++){

    if (sigmaIsList) curSigma = Rcpp::as<arma::mat>(sigmaList[i]);
    if (kappaIsList) curKappa = Rcpp::as<arma::mat>(kappaList[i]);
    if (muIsList) curMu = Rcpp::as<arma::vec>(muList[i]);

    // Likelihood:
    result += logLikelihood_gaussian_subgroup_fiml_cpp_inner(
      curSigma,
           curKappa,
                curMu,
                  fimldata[i],
                          epsilon) ;

    // nvar <- ncol(kappa)
    // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  }


  // Return
  return result;
}















// 
// 
// // Main function
// // [[Rcpp::export]]
// double logLikelihood_gaussian_cpp(
//     const Rcpp::List& prep
// ){
//   
//   Rcpp::List groupmodels = prep["groupModels"];
//   int nGroup = groupmodels.length();
//   arma::vec nPerGroup = prep["nPerGroup"];
//   double nTotal = prep["nTotal"];
//   bool fullFIML = groupmodels["fullFIML"];
//   
//   // Result:
//   double ll = 0;
//   
//   for (int i=0; i<nGroup;i++){
//     if (fullFIML){
//       ll += (nPerGroup(i) / 2) * logLikelihood_gaussian_group_cpp_fullFIML_outer(groupmodels[i]);
//     } else {
//       ll += (nPerGroup(i) / 2) * logLikelihood_gaussian_group_cpp_outer(groupmodels[i]);
//     }
//     
//   }
//   
//   return(ll);
// }
// 
// 
// 




