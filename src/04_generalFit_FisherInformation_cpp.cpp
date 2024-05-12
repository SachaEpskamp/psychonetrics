// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "05_MLestimator_expected_Hessian_Gauss_cpp.h"
#include "06_ULS_expectedHessian_cpp.h"
#include "07_FIMLestimator_expected_hessian_gauss_cppversion.h"
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
arma::mat FisherInformation_inner_cpp_DSS(
    const arma::mat& estimator,
    const arma::sp_mat& model,
    const arma::sp_mat& manual
) {
  // Sparse part first:
  arma::sp_mat sparse = model * manual;
  
  // Now second part:
  arma::mat Fis = 0.5 * sparse.t() * estimator * sparse;
  
  // Return
  return Fis;
}
// 0.5 * t(manualPart) %*% t(modelPart) %*% estimatorPartHessian %*% modelPart %*% manualPart

// Dense - dense - sparse
// [[Rcpp::export]]
arma::mat FisherInformation_inner_cpp_DDS(
    const arma::mat& estimator,
    const arma::mat& model,
    const arma::sp_mat& manual
) {
  // Dense part first:
  arma::mat dense = model.t() * estimator * model;
  
  // Now full:
  arma::mat Fis = 0.5 * manual.t() * dense * manual;
  
  return Fis;
}





// [[Rcpp::export]]
void psychonetrics_FisherInformation_cpp_inner(
    const arma::vec& x,
    arma::mat& Fisher,
    const S4& model,
    bool useM = true,
    bool sparsemodel = false
){
  // Prepare model:
  Rcpp::List prep = prepareModel_cpp(x, model);

  
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

      estimatorPart = expected_hessian_Gaussian_cpp(prep);

    } else if (distribution == "Ising"){

      // Obtain environment containing function
      Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
      
      // Make function callable from C++
      Rcpp::Function hesFun = base["expected_hessian_Ising"]; 
      
      estimatorPart = as<arma::mat>(hesFun(prep));

    } else {
      Rf_error("Distribution not supported for ML estimator.");
    }

  } else if (estimator == "ULS" || estimator == "WLS" || estimator == "DWLS"){

    estimatorPart = expected_hessian_ULS_Gaussian_cpp(prep);

  } else if (estimator == "FIML"){

    estimatorPart = expected_hessian_fiml_Gaussian_cppVersion(prep);

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
  

  // Compute the FisherInformation
  // FIXME? Use sparse?
  
  if (useM){
    // Manual part:
    arma::sp_mat manualPart = Mmatrix_cpp_list(model.slot("parameters"));
    
    if (sparsemodel) {
      
      Fisher = FisherInformation_inner_cpp_DSS(estimatorPart, (arma::sp_mat)modelPart, manualPart);
      
    } else {
    
    Fisher = FisherInformation_inner_cpp_DDS(estimatorPart, modelPart, manualPart);
      
    }
    
    
    
  } else {
    
    // Avoid using the M matrix:
    
    // Compute inner part:
    // Fixme: Potentially sparse? Or use block structure in multigroup setting!!!
    arma::mat innerpart;
    if (sparsemodel) {
      arma::sp_mat sparse_model = (arma::sp_mat)modelPart;
      
      innerpart =  sparse_model.t() * estimatorPart * sparse_model;
      
    } else {
      
      innerpart = modelPart.t() * estimatorPart * modelPart;
      
    }
      
      
      
    
    // Get parameter table:
    Rcpp::List pars = model.slot("parameters");
    arma::vec parnum = pars["par"];
    
    // ints to use:
    int i,j;
    
    // Number of free parameters:
    // int freePar = max(parnum);
    
    // Number of total parameters:
    int totalPar = parnum.n_elem;
    
    // Empty Fisher:
    Fisher.fill(0);
    
    // Start looping:
    for (i=0;i<totalPar;i++){
      for (j=0;j<totalPar;j++){
        if (parnum(i) > 0 && parnum(j) > 0){
          
          Fisher(parnum(i)-1, parnum(j)-1) += innerpart(i,j);
          
        }
      }
    }
    
    // Finally, mulyiply with 0.5:
    Fisher *= 0.5;
    
  }

}




// [[Rcpp::export]]
arma::mat psychonetrics_FisherInformation_cpp(
    const S4& model,
    bool useM = false,
    bool sparsemodel = false
){

  arma::vec x(parVector_cpp(model));
  
  // Create empty FisherInformation:
  arma::mat Fisher = zeros(x.n_elem, x.n_elem);

  // Run inner:
  psychonetrics_FisherInformation_cpp_inner(x, Fisher, model, useM, sparsemodel);

  
  return(Fisher);
}


