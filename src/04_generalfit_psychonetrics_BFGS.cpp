// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "04_generalfit_fitfunction_cpp.h"
#include "04_generalFit_gradient_cpp.h"
#include "b_modelexpansion_updateModel_cpp.h"
#include "02_algebrahelpers_modelMatrix_cpp.h"
#include "roptim.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace roptim;



class psychonetrics_optim : public Functor {
public:
  psychonetrics_optim(const S4& model) : model_(model) {
    os.sann_use_custom_function_ = true;
  }
  
  double operator()(const arma::vec &x) override {

    double res = psychonetrics_fitfunction_cpp(x, model_);
    
    return res;
  }
  
  void Gradient(const arma::vec &x, arma::vec &gr) override {
    return psychonetrics_gradient_cpp_inner(x, gr, model_);
  }
  
private:
  S4 model_;
};




// [[Rcpp::export]]
S4 psychonetrics_optimizer(
    S4 model,
    const arma::vec& lower,
    const arma::vec& upper,
    std::string optimizer = "L-BFGS-B",
    bool bounded = false
) {
  int i;
  
  // Parameters:
  Rcpp::List pars = model.slot("parameters");
  
  // Total pars:
  arma::vec est = pars["est"];
  arma::vec parnum = pars["par"];
  int nPar = max(parnum);
  int totalPar = parnum.n_elem;
  
  // If 0, nothing to do!
  if (nPar == 0){
    model.slot("computed") = true;
    return(model);
  }
  
  // Form the start vector:
  arma::vec x(nPar);
  
  for (i=0;i<totalPar;i++){
    if (parnum(i) > 0){
      x(parnum(i)-1) = est(i);
    }
  }
  
  
  // Optimize:
  psychonetrics_optim rb(model);
  Roptim<psychonetrics_optim> opt(optimizer);
  
  // Control pars:
  opt.control.trace = 0; // <- no output
  opt.control.maxit = 20000L; 
  
  // 
  // opt.control.abstol = 1.490116e-08; // Same as nlminb
  // opt.control.reltol = 1.490116e-08; // Same as nlminb
  // //opt.control.ndeps = ones(nPar);
  // opt.control.pgtol = 1.490116e-08; // Same as nlminb
  // opt.control.factr = 0; // Same as nlminb
  // 
  // 
  
  
  opt.control.abstol = 1.490116e-08; // Same as nlminb
  opt.control.reltol = 1.490116e-08; // Same as nlminb
  opt.control.ndeps = ones(nPar);
  opt.control.pgtol = 1.490116e-08; // Same as nlminb
  opt.control.factr = 1.490116e-08; // Same as nlminb
  

  
  // Bounds:
  if (optimizer == "L-BFGS-B" && bounded){
    opt.set_lower(lower);
    opt.set_upper(upper);    
  }

  
  opt.set_hessian(false);

  
  // OPTIMIZE:
  opt.minimize(rb, x);

  // Output info:
  Rcpp::List optimout = List::create(
    Named("par") = opt.par() , 
    Named("convergence") = opt.convergence(),
    Named("message") = opt.message(),
      Named("value") = opt.value(),
      Named("fncount") = opt.fncount(),
      Named("grcount") = opt.grcount(),
      Named("optimizer") = optimizer
    );
    
    
    // Update model:
    model = updateModel_cpp(x,model,false);

    // Set computed:
    model.slot("computed") = true;
    
    // Add output:
    model.slot("optim") = optimout;
  
  return(model);
}
