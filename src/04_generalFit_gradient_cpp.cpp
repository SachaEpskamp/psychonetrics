// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "05_MLestimator_gradient_Gauss_cpp.h"
#include "05_MLestimator_gradient_Ising_cpp.h"
#include "06_ULS_gradient_cpp.h"
#include "07_FIMLestimator_jacobian_gauss.h"
#include "14_varcov_derivatives_cpp.h"
#include "15_lvm_derivatives_cpp.h"
#include "16_var1_derivatives_cpp.h"
#include "18_dlvm1_derivatives_cpp.h"
#include "19_tsdlvm1_derivatives_cpp.h"
#include "20_meta_varcov_derivatives_cpp.h"
#include "21_Ising_derivatives_cpp.h"
#include "22_ml_lvm_derivatives_cpp.h"
#include "02_algebrahelpers_modelMatrix_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Dense - sparse - sparse
// [[Rcpp::export]]
arma::mat gradient_inner_cpp_DSS(
    const arma::mat& estimator,
    const arma::sp_mat& model,
    const arma::sp_mat& manual
) {
  // Sparse part first:
  arma::sp_mat sparse = model * manual;
  
  // Now second part:
  arma::vec grad = vectorise(estimator * sparse);
  
  // Return
  return grad;
}

// Dense - dense - sparse
// [[Rcpp::export]]
arma::mat gradient_inner_cpp_DDS(
    const arma::mat& estimator,
    const arma::mat& model,
    const arma::sp_mat& manual
) {
  arma::mat part1 = estimator * model;
  // Return
  arma::vec grad = vectorise(part1 * manual);
  
  return grad;
}





// [[Rcpp::export]]
void psychonetrics_gradient_cpp_inner(
    const arma::vec& x,
    arma::vec& grad,
    const S4& model,
    bool sparsemodel = false,
    bool useM = false
){
  // Prepare model:
  Rcpp::List prep = prepareModel_cpp(x, model);
  
  // Manual part:
  arma::sp_mat manualPart = Mmatrix_cpp_list(model.slot("parameters"));
  
  
  // What estimator:
  std::string estimator= prep["estimator"];
  
  
  // What distribution:
  std::string distribution = prep["distribution"];
  
  // What model:
  std::string usemodel = prep["model"];
  
  // Estimator part:
  arma::mat estimatorPart;
  
  if (estimator == "ML"){
    if (distribution == "Gaussian"){
      
      estimatorPart = jacobian_gaussian_sigma_cpp(prep);
      
    } else if (distribution == "Ising"){
      
      estimatorPart = jacobian_Ising_cpp(prep);
      
    } else {
      Rf_error("Distribution not supported for ML estimator.");
    }
    
  } else if (estimator == "ULS" || estimator == "WLS" || estimator == "DWLS"){
    
    estimatorPart = ULS_gradient_Gauss_cpp(prep);
    
  } else if (estimator == "FIML"){
    
    estimatorPart = jacobian_fiml_gaussian_sigma_cpp(prep);
    
  } else {
    Rf_error("Estimator not supported.");
  }
  
  // Model part:
  arma::mat modelPart;
  
  if (usemodel == "varcov"){
    
    modelPart = d_phi_theta_varcov_cpp(prep);
    
  } else   if (usemodel == "lvm"){
    
    modelPart = d_phi_theta_lvm_cpp(prep);
    
  } else  if (usemodel == "var1"){
    
    modelPart = d_phi_theta_var1_cpp(prep);
    
  }  else  if (usemodel == "dlvm1"){
    
    modelPart = d_phi_theta_dlvm1_cpp(prep);
    
  } else  if (usemodel == "tsdlvm1"){
    
    modelPart = d_phi_theta_tsdlvm1_cpp(prep);
    
  }  else if (usemodel == "meta_varcov"){
    
    modelPart = d_phi_theta_meta_varcov_cpp(prep);
    
  }  else if (usemodel == "Ising"){
    
    modelPart = d_phi_theta_Ising_cpp(prep);
    
  }  else if (usemodel == "ml_lvm"){
    
    modelPart = d_phi_theta_ml_lvm_cpp(prep);
    
  }
  
  // Compute the gradient
  // FIXME? Use sparse?
  
  // Get parameter table:
  Rcpp::List pars = model.slot("parameters");
  arma::vec parnum = pars["par"];
  int freePar = max(parnum);
  
  arma::mat Jac(1,freePar);
  
  
  if (useM){
    if (sparsemodel) {
      arma::sp_mat sparse_model = (arma::sp_mat)modelPart;
      
      Jac = gradient_inner_cpp_DSS(estimatorPart, sparse_model, manualPart);
      
    } else {
      
      Jac = gradient_inner_cpp_DDS(estimatorPart, modelPart, manualPart);
      
    }    
    
  }  else {
    
    
    // Avoid using the M matrix:
    
    // Compute inner part:
    // Fixme: Potentially sparse? Or use block structure in multigroup setting!!!
    arma::mat innerpart;
    
    if (sparsemodel) {
      arma::sp_mat sparse_model = (arma::sp_mat)modelPart;
      
      innerpart =  estimatorPart * sparse_model;
      
    } else {
      
      innerpart = estimatorPart * modelPart;
      
    }
    
    
    // Get parameter table:
    Rcpp::List pars = model.slot("parameters");
    arma::vec parnum = pars["par"];
    
    // ints to use:
    int i;
    
    // Number of free parameters:
    // int freePar = max(parnum);
    
    // Number of total parameters:
    int totalPar = parnum.n_elem;
    
    // Empty Fisher:
    Jac.fill(0);
    
    // Start looping:
    for (i=0;i<totalPar;i++){
        if (parnum(i)){
          
          Jac(0, parnum(i)-1) += innerpart(0,i);
          
        }
    }
    
    
    
  }
  
  
  grad = vectorise(Jac);
  // return(vectorise(Jac));
}




// [[Rcpp::export]]
arma::vec psychonetrics_gradient_cpp(
    arma::vec x,
    const S4& model,
    bool useM = false,
    bool sparsemodel = false
){
  
  
  // Create empty gradient:
  arma::vec grad(x.n_elem);
  
  // Run inner:
  psychonetrics_gradient_cpp_inner(x, grad, model, sparsemodel, false);
  
  
  return(grad);
}

