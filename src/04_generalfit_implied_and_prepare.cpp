// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "14_varcov_implied_cpp.h"
#include "14_varcov_prepare_cpp.h"
#include "15_lvm_prepare_cpp.h"
#include "15_lvm_implied_cpp.h"
#include "16_var1_implied_cpp.h"
#include "16_var1_prepare_cpp.h"
#include "18_dlvm1_implied_cpp.h"
#include "18_dlvm1_prepare_cpp.h"
#include "19_tsdlvm1_implied_cpp.h"
#include "19_tsdlvm1_prepare_cpp.h"
#include "20_meta_varcov_implied_cpp.h"
#include "20_meta_varcov_prepare_cpp.h"
#include "21_Ising_implied_cpp.h"
#include "21_Ising_prepare_cpp.h"
#include "22_ml_lvm_prepare_cpp.h"
#include "22_ml_lvm_implied_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List impliedModel_cpp(
    const S4& model,
    bool all = false
){
  // What model:
  std::string framework= model.slot("model");
  
  // Output list:
  Rcpp::List imp;
  
  
  // Run:
  if (framework == "varcov"){
    
    imp = implied_varcov_cpp(model, all); // = Updated!
    
  } else if (framework == "lvm"){
    
    imp = implied_lvm_cpp(model, all); // = Updated!
    
  } else if (framework == "var1"){
    
    imp = implied_var1_cpp(model, all); // = Updated!
    
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function impfun = base["implied_var1"]; 
    // 
    // imp = impfun(Rcpp::_["model"] = model, Rcpp::_["all"] =   all);
    
  } else if (framework == "dlvm1"){
    
    imp = implied_dlvm1_cpp(model, all); // = Updated!
    
  }  else if (framework == "tsdlvm1"){

    imp = implied_tsdlvm1_cpp(model, all); // = Updated!
    
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function impfun = base["implied_tsdlvm1"]; 
    // 
    // imp = impfun(Rcpp::_["model"] = model, Rcpp::_["all"] =   all);
    
  }   else if (framework == "meta_varcov"){
    
    imp = implied_meta_varcov_cpp(model, all); // = Updated!
    // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function impfun = base["implied_meta_varcov"]; 
    // 
    // imp = impfun(Rcpp::_["model"] = model, Rcpp::_["all"] =   all);
    
  }  else if (framework == "Ising"){
    
    imp = implied_Ising_cpp(model, all); // = Updated!
    
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function impfun = base["implied_Ising"]; 
    // 
    // imp = impfun(Rcpp::_["model"] = model, Rcpp::_["all"] =   all);
    
  }  else if (framework == "ml_lvm"){
    
    
    imp = implied_ml_lvm_cpp(model, all); // = Updated!
    // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function impfun = base["implied_ml_lvm"]; 
    // 
    // imp = impfun(Rcpp::_["model"] = model, Rcpp::_["all"] =   all);
    // 
  }
  
  // Return:
  return(imp);
}


// [[Rcpp::export]]
Rcpp::List prepareModel_cpp(
    arma::vec x,
    const S4& model
){
  int g;
  
  // What model:
  std::string framework= model.slot("model");
  
  // Output list:
  Rcpp::List prep;
  
  
  // Run:
  if (framework == "varcov"){
    
    prep = prepare_varcov_cpp(x, model); // = Updated!
    
  } else if (framework == "lvm"){
    
    prep = prepare_lvm_cpp(x, model); // = Updated!
    
  } else if (framework == "var1"){
    
    prep = prepare_var1_cpp(x, model); // = Updated!
    
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ;
    // 
    // // Make function callable from C++
    // Rcpp::Function prepfun = base["prepare_var1"];
    // 
    // prep = prepfun(Rcpp::_["x"] = x, Rcpp::_["model"] =   model);
    
  } else if (framework == "dlvm1"){
    
    prep = prepare_dlvm1_cpp(x, model); // = Updated!
    
    
  }  else if (framework == "tsdlvm1"){
    
    prep = prepare_tsdlvm1_cpp(x, model); // = Updated!
    
    // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function prepfun = base["prepare_tsdlvm1"]; 
    // 
    // prep = prepfun(Rcpp::_["x"] = x, Rcpp::_["model"] =   model);
    // 
  }   else if (framework == "meta_varcov"){
    
    prep = prepare_meta_varcov_cpp(x, model); // = Updated!
    // 
    // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function prepfun = base["prepare_meta_varcov"]; 
    // 
    // prep = prepfun(Rcpp::_["x"] = x, Rcpp::_["model"] =   model);
    
  }  else if (framework == "Ising"){
    
    prep = prepare_Ising_cpp(x, model); // = Updated!
    // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function prepfun = base["prepare_Ising"]; 
    // 
    // prep = prepfun(Rcpp::_["x"] = x, Rcpp::_["model"] =   model);
    
  }  else if (framework == "ml_lvm"){
    
    prep = prepare_ml_lvm_cpp(x, model); // = Updated!
    // // 
    // // Obtain environment containing function
    // Rcpp::Environment base = Environment::namespace_env( "psychonetrics" ) ; 
    // 
    // // Make function callable from C++
    // Rcpp::Function prepfun = base["prepare_ml_lvm"]; 
    // 
    // prep = prepfun(Rcpp::_["x"] = x, Rcpp::_["model"] =   model);
    // 
  }
  
  // Sample
  S4 sample = model.slot("sample");
  Rcpp::List fimldata = sample.slot("fimldata");
  Rcpp::List WLS_W = sample.slot("WLS.W");
  
  // PArameters:
  Rcpp::List parList = model.slot("parameters");
  arma::vec pars = parList["par"];
  
  // Groups
  Rcpp::List groups = sample.slot("groups");
  arma::vec nobspergroup = groups["nobs"];
  
  // Number of groups:
  Rcpp::List groupModels = prep["groupModels"];
  int nGroup = groupModels.length();
  
  // Estimator:
  std::string estimator = model.slot("estimator");
  
  // If the estimator is FIML, add the raw data:
  for (g=0;g<nGroup;g++){
    Rcpp::List groupmod = groupModels[g];
    
    groupmod["cpp"] = true;
    
    if (estimator == "FIML"){
      // Add the raw data to each group:
      groupmod["fimldata"] = fimldata[g];
      groupmod["fulln"] = nobspergroup[g];
      
    } else   if (estimator == "WLS" || estimator == "DWLS" || estimator == "ULS"){
      groupmod["WLS.W"] = WLS_W[g];
      
    }
    
    groupmod["estimator"] = estimator;
    
    // Return back:
    groupModels[g] = groupmod;
  }
  
  // Write back:
  prep["groupModels"] = groupModels;
  
  prep["fullFIML"] = sample.slot("fullFIML");
  
  // FIXME Add WLS.W:
  
  
  // Add number of parameters:
  prep["nParFull"] = pars.n_elem;
  prep["nParFree"] = max(pars);
    
    prep["estimator"] = model.slot("estimator");
    prep["distribution"] = model.slot("distribution");
    prep["model"] = model.slot("model");
    
    
    // Return:
    return(prep);
}
