#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "02_algebrahelpers_RcppHelpers.h"


// Inner function
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    const Rcpp::List& dat,
    double epsilon) {
  // Rf_PrintValue(wrap("USED"));
  
  double n_part;
  
  // Number of parameters
  // int nmeans = mu.size();
  // int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  // arma::mat Jac = zeros(1, nvars + nmeans);
  
  // Integers
  // int i;
  // int j;
  
  // Loop over groups
  // for ( i = 0; i < fimldata.size(); i++){
  // Obtain list:
  // Rcpp::List dat = fimldata[i];
  
  // Observed indices:
  vec obs = dat["obs"];
  
  uvec inds = find(obs == true);
  
  // sample size:
  n_part = dat["n"];
  
  // Duplication matrix:
  arma::sp_mat D = dat["D"]; 
  arma::sp_mat L = dat["L"]; 
  
  // Subset matrices:
  arma::mat sigma_p = sigma(inds,inds);
  arma::mat kappa_p(sigma_p);
  arma::vec mu_p = mu(inds);
  
  
  // Observed values:
  arma::mat S = dat["S"];
  arma::vec means = dat["means"];
  
  // inverse:
  // arma::vec ev = arma::eig_sym(sigma_p);
  // bool ispos = true;
  // for ( j = 0; j < ev.size(); j++){
  //   if (ev[j] < sqrt(epsilon)){
  //     ispos = false;
  //     break;
  //   }
  // }
  // if (ispos){
  //   kappa_p = inv(sigma_p);
  // } else {
  //   kappa_p = pinv(sigma_p);
  // }
  // log det:
  // arma::vec ev = arma::eig_sym(sigma_p);
  // bool ispos = ev[0] > epsilon;
  // 
  // // bool ispos = true;
  // // for (int j = 0; j < ev.size(); j++){
  // //   if (ev[j] < sqrt(epsilon)){
  // //     ispos = false;
  // //     break;
  // //   }
  // // }
  // double logepsilon = log(epsilon);
  // 
  // if (ispos){
  //   kappa_p = inv(sigma_p);
  // } else {
  //   kappa_p = pinv(sigma_p);
  // }
  kappa_p = solve_symmetric_cpp_matrixonly(sigma_p, epsilon);
  
  
  // Mean part:
  arma::mat meanpart = -2 * (means - mu_p).t() * kappa_p;
  // Sigma part:
  arma::mat matpart = kappa_p * (S + (means - mu_p) * (means - mu_p).t() - sigma_p) * kappa_p;
  
  //Vectorize:
  arma::vec vecpart = vectorise(matpart);
  arma::mat sigmapart = -1.0 * (D.t() * vecpart).t();
  
  // Join:
  arma::mat Jac =  n_part * join_rows(meanpart, sigmapart) * L;
  
  // }
  
  // Return
  return Jac;
}


// Outer function
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Jac = zeros(1, nvars + nmeans);
  
  // Integers
  int i;

  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){

    // Join:
    Jac += jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
       epsilon);
    
  }
  
  // Return
  return Jac;
}


// Outer function
// sigma, kappa and mu may each be supplied either as a per-row list (one
// element per row of the data) or as a single matrix/vector that applies to
// every row (e.g. the varcov family). Both cases are handled here without
// materializing per-row copies.
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp_fullFIML(
    SEXP sigma,
    SEXP kappa,
    SEXP mu,
    const Rcpp::List& fimldata,
    double epsilon) {

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

  // Number of parameters
  arma::vec firstmu = muIsList ? Rcpp::as<arma::vec>(muList[0]) : curMu;
  int nmeans = firstmu.size();
  int nvars = nmeans * (nmeans + 1) / 2;

  // Empty Jacobian:
  arma::mat Jac = zeros(1, nvars + nmeans);

  // Integers
  int i;


  // Loop over rows
  for ( i = 0; i < fimldata.size(); i++){

    if (sigmaIsList) curSigma = Rcpp::as<arma::mat>(sigmaList[i]);
    if (kappaIsList) curKappa = Rcpp::as<arma::mat>(kappaList[i]);
    if (muIsList) curMu = Rcpp::as<arma::vec>(muList[i]);

    // Join:
    Jac += jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
      curSigma,
      curKappa,
      curMu,
      fimldata[i],
              epsilon);

  }

  // Return
  return Jac;
}



arma::mat jacobian_fiml_outer_cpp(
    const Rcpp::List& grouplist,
    bool fullFIML
){
  double fulln = grouplist["fulln"];
  
  bool corinput = false;
  if (grouplist.containsElementNamed("corinput")){
    corinput = grouplist["corinput"];
  }

  bool meanstructure = true;
  
  if (grouplist.containsElementNamed("meanstructure")){
    meanstructure = grouplist["meanstructure"];
  }
  
  arma::mat S = grouplist["S"];
  int nvar = S.n_rows;
  arma::mat Jac;
  
  // Subgroup models:
  if (fullFIML){
   Jac = 1.0/fulln * jacobian_fiml_gaussian_subgroup_sigma_cpp_fullFIML(grouplist["sigma"],grouplist["kappa"],grouplist["mu"],grouplist["fimldata"],1.490116e-08);
    
  } else {
    
    Jac = 1.0/fulln * jacobian_fiml_gaussian_subgroup_sigma_cpp(grouplist["sigma"],grouplist["kappa"],grouplist["mu"],grouplist["fimldata"],1.490116e-08);
  
  }

  
  
  // FIXME: Nicer to not bother with computing the parts not needed
  // Cut out variances if needed:
  if (corinput){
 
      arma::mat I = eye(nvar, nvar );
      arma::vec dummyvec = join_cols(
        zeros<vec>(nvar),
        vech(I)
      );
      
      uvec remove = find(dummyvec > 0);
      
      Jac.shed_cols(remove);      
    
  }
  
  if (!meanstructure){

    Jac = Jac.submat(0, nvar, Jac.n_rows-1, Jac.n_cols-1);
    
  } 
  
  return(Jac);
  
}



// full Jacobian function 
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_sigma_cpp(
    const Rcpp::List& prep
){
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  bool fullFIML = prep["fullFIML"];
  
  // JAcobian:
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    // if (fullFIML){
      
      arma::mat groupgrad =  (nPerGroup(i) / nTotal) * jacobian_fiml_outer_cpp(groupmodels[i], fullFIML);
      groupgradients[i]  = groupgrad;      
      

  }
  
  
  arma::mat res =  cbind_psychonetrics(groupgradients);
  
  
  return(res);
}






