// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <cstddef>
#include <algorithm>
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
    //return psychonetrics_fitfunction_cpp(x, model_);
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
    std::string optimizer = "L-BFGS-B"
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
  opt.control.maxit = 20000; 
  opt.control.abstol = 1.490116e-08 * 10;
  opt.control.reltol = 1e-3;
  opt.control.pgtol = 1e-5;
  // opt.control.reltol = 1.490116e-08;
  // opt.control.abstol = 1.490116e-08;
  
  // Bounds:
  if (optimizer == "L-BFGS-B"){
    opt.set_lower(lower);
    opt.set_upper(upper);    
  }

  
  opt.set_hessian(false);
//   
//   optim.control$control<- list(eval.max=20000L,
//                                iter.max=10000L,
//                                rel.tol=1e-5,
// #step.min=2.2e-14, # in =< 0.5-12
//                                step.min=1.0, # 1.0 in < 0.5-21
//                                  step.max=1.0,
//                                    x.tol=1.5e-8,
//                                    xf.tol=2.2e-14)
  
  // OPTIMIZE:
  opt.minimize(rb, x);

  
  // Rcpp::Rcout << "-------------------------" << std::endl;
  // opt.print();

  // Output info:
  Rcpp::List optimout = List::create(
    Named("par") = opt.par() , 
    Named("convergence") = opt.convergence(),
    Named("message") = opt.message(),
      Named("value") = opt.value(),
      Named("fncount") = opt.fncount(),
      Named("grcount") = opt.grcount()
    );
    
    // Update model:
    model = updateModel_cpp(x,model,false);

    // Set computed:
    model.slot("computed") = true;
    
    // Add output:
    model.slot("optim") = optimout;
  
  return(model);
}


// Old test:
// // [[Rcpp::export]]
// S4 psychonetrics_BFGS(
//     const S4& model,
//     arma::mat Hstart
// ){
//   int i;
//   
//   // Copy model:
//   S4 newMod(model);
//   
//   // Parameters:
//   Rcpp::List pars = newMod.slot("parameters");
//   
//   // Total pars:
//   arma::vec est = pars["est"];
//   arma::vec parnum = pars["par"];
//   int nPar = max(parnum);
//   
//   // If 0, nothing to do!
//   if (nPar == 0){
//     newMod.slot("computed") = true;
//     return(newMod);    
//   }
//   
//   // Form the start vector:
//   arma::vec x(nPar);
//   
//   for (i=0;i<nPar;i++){
//     if (parnum(i) > 0){
//       x(parnum(i)-1) = est(i);
//     }
//   }
//   
//   
//   // Manual part:
//   arma::sp_mat manualPart = Mmatrix_cpp_list(newMod.slot("parameters"));
//   
//   // Initial Hessian:
//   arma::mat Hes = Hstart; //eye(nPar, nPar);
//   
//   // Setup things I need:
//   arma::vec curGrad;
//   arma::vec newGrad;
//   arma::vec p;
//   arma::vec y;
//   double alpha_k ;
//   double dummyfit ;
//   double propfit;
//   double propalpha;
//   bool alphasearchcont;
//   arma::vec denominator;
//   
//   // Converged:
//   bool converged = false;
//   
//   int nIter = 0;
//   // START OPTIMIZER:
//   do{
//     nIter++;
//     // Current gradient:
//     curGrad = psychonetrics_gradient_cpp(x, newMod);
//     
//     // Step 1: Find initial direction:
//     p = solve(Hes, -curGrad);
//     
//     // STep 2:Find acceptable stepsize:
//     // FIXME: This is a silly optimizer ....
//     alpha_k = 1;
//     
//     // dummyfit = psychonetrics_fitfunction_cpp(x + alpha_k * p, newMod);
//     // propfit;
//     // propalpha;
//     // 
//     // alphasearchcont = true;
//     // 
//     // do{
//     //   propalpha = alpha_k / 2;
//     //   propfit = psychonetrics_fitfunction_cpp(x + propalpha * p, newMod);
//     // 
//     //   if (dummyfit < propfit){
//     //     alphasearchcont = false;
//     //   } else{
//     //     alpha_k = propalpha;
//     //     dummyfit = propfit;
//     //   }
//     // } while (alphasearchcont);
//     
//     
//     // Step 3: update x:
//     arma::vec s = alpha_k * p;
//     x = x + s;
//     
//     // Step 4, set y:
//     newGrad = psychonetrics_gradient_cpp(x, newMod);
//     y = newGrad - curGrad;
//     curGrad = newGrad;
//     
//     // Step 5: Update Hes estimate:
//     denominator = (s.t() * Hes * s);
//       Hes = Hes + (y * y.t()) / dot(y,s) - ( Hes * s * s.t() * Hes.t() )/ denominator(0);
//       
//     // Check if converged:
//     converged = sum(abs(curGrad)) < 0.01;
//     
//   } while (!converged);
//   
//   
//   // Update model:
//   newMod = updateModel_cpp(x,newMod,false);
//   
//   // Set computed:
//   newMod.slot("computed") = true;
//   
//   // Rf_PrintValue(wrap(nIter));
//   
//   // Return:
//   return(newMod);
// }