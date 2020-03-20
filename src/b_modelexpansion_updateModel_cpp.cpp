// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// 
// // we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include <RcppArmadillo.h>
// #include <math.h>
// #include "02_algebragelpers_kronecker.h"
// #include "02_algebrahelpers_RcppHelpers.h"
// #include "03_modelformation_formModelMatrices_cpp.h"
// #include "03_modelformation_impliedcovstructures.h"
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
// 
// // [[Rcpp::export]]
// S4 updateModel_cpp(
//     arma::vec x,
//     S4& model,
//     bool updateMatrices
// ){
//   S4 newMod(model); // FIXME: this copies the entire model ...
//   
//   Rcpp::List parsList = newMod.slot("parameters");
//   
//   Rcpp::DataFrame parsDF = CheapDataFrameBuilder(parsList);
//   
//   int nPar = x.n_elem;
//   
//   
//   newMod.slot("parameters") = parsDF;
//   return(newMod);
// }
// 
