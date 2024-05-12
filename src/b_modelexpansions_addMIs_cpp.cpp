// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// 
// // we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include <RcppArmadillo.h>
// #include <math.h>
// #include "02_algebragelpers_kronecker.h"
// #include "02_algebrahelpers_RcppHelpers.h"
// 
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
// 
// 
// // [[Rcpp::export]]
// S4 addMIs_inner_full(
//     const S4& xOld,
//     std::string type
// ){
//   // Ints:
//   int i;
//   
//   
//   // clone model:
//   S4 x = clone(xOld); // FIXME: creates a copy, but it avoids a ton of weird stuff happening otherwise... Could do better.
//   
//   // Parameter table:
//   Rcpp::List parameters = x.slot("parameters");
//   
//   // Parameter indices and estimates:
//   NumericVector est = parameters["est"];
//   IntegerVector par = parameters["par"];
//   
//   NumericVector mi = parameters["mi"];
//   NumericVector pmi = parameters["pmi"];
//   NumericVector epc = parameters["epc"];
//   
//   
//   NumericVector mi_free = parameters["mi_free"];
//   NumericVector pmi_free = parameters["pmi_free"];
//   NumericVector epc_free = parameters["epc_free"];
//   
//   NumericVector mi_equal = parameters["mi_equal"];
//   NumericVector pmi_equal = parameters["pmi_equal"];
//   NumericVector epc_equal = parameters["epc_equal"];
//   
//   int nparTotal = par.n_elem;
//   int nparFree = max(par);
//   
//   
//   // If no constrained parameters, nothing to do!
//   if (nparTotal == nparFree){
//     
//     return(x);
//     
//   }
//   for (i=0;i>nparTotal;i++ ){
//     
//     // Clear old MIs:
//     if (type == "normal"){
//       mi[i] = 0;
//       pmi[i] = NA_REAL;
//       epc[i] = NA_REAL;
//     } else if (type == "free"){
//       mi_free[i] = 0
//       pmi_free[i] = NA_REAL;
//       epc_free[i] = NA_REAL;
//     } else {
//       mi_equal[i] = 0
//       pmi_equal[i] = NA_REAL;
//       epc_equal[i] = NA_REAL;
//     }
//   }
//   
//   
//   // Sample size:
//   sample = x.slot("sample");
//   groups = sample.slot("groups");
//   NumericVector nobs = groups["nobs"];
//   int n = sum(nobs);
//     
//     // Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but constrained per group #
//     // Fully free:
//     // Copy the model:
//     S4 modCopy = clone(x);
//     
//     // Obtain the full set of parameters that are constrained across all groups:
//     if (type == "equal"){
//       sum = modCopy@parameters %>% group_by_("matrix","row","col") %>% summarize_(anyConstrained = ~any(fixed)) %>% 
//         filter_(~anyConstrained)
//       // Add a unique number to each:
//       sum$par2 =  max(modCopy@parameters$par) + seq_len(nrow(sum))
//       
//       // Left join back:
//       modCopy@parameters = modCopy@parameters %>% left_join(sum,by=c("matrix","row","col")) %>% 
//         mutate_(par = ~ifelse(identified,0,ifelse(par==0,par2,par)))
//       
//     } else {
//       // Add free parameter numbers:
//       modCopy@parameters$par[modCopy@parameters$par==0 & !modCopy@parameters$identified] = max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$par==0 & !modCopy@parameters$identified))
//       
//       // For each group, free all parameters from equality constraints:
//       if (type == "free"){
//         if (nrow(modCopy@sample@groups)>1){
//           for (g in 2:nrow(modCopy@sample@groups)){
//             modCopy@parameters$par[modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified] = max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified))
//           }
//         }      
//       }
//       
//     }
//     
//     // Identify:
//     // modCopy = identify(modCopy)
//     
//     // Remake the model matrix:
//     // modCopy@extramatrices$M = Mmatrix(modCopy@parameters)
//     // 
//     if (modCopy@cpp){
//       modCopy@extramatrices$M   = Mmatrix_cpp(modCopy@parameters )
//     } else {
//       modCopy@extramatrices$M  = Mmatrix(modCopy@parameters)  
//     }
//     
//     
//     // Check if gradient and hessian are present:
//     // gradient = !is.null(x@fitfunctions$gradient)
//     // hessian = !is.null(x@fitfunctions$hessian)
//     // information = !is.null(x@fitfunctions$information)
//     
//     // Compute a gradient:
//     if (modCopy@cpp){
//       g = psychonetrics_gradient_cpp(parVector(modCopy), modCopy)
//     } else {
//       g = psychonetrics_gradient(parVector(modCopy), modCopy)
//       
//     }
//     // if (gradient){
//     //   g = x@fitfunctions$gradient(parVector(modCopy), modCopy)
//     // } else {
//     //   g = numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
//     // }
//     
//     // Compute the Fisher information:
//     // if (information){
//     //   H = x@fitfunctions$information(modCopy)
//     // } else if (hessian){
//     //   H = x@fitfunctions$hessian(parVector(modCopy), modCopy)
//     // } else if (gradient){
//     //   H = numDeriv::jacobian(x@fitfunctions$gradient,parVector(modCopy), model=modCopy) 
//     // } else {
//     //   H = numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
//     // }
//     // Total sample size:
//     nTotal = sum(x@sample@groups$nobs)
//       
//       // FIXME: 4 * n could be nicer here probably
//       if (modCopy@cpp){
//         H = 4 * nTotal * as(psychonetrics_FisherInformation_cpp(modCopy, analyticFisher), "matrix")
//       } else {
//         H = 4 * nTotal * as(psychonetrics_FisherInformation(modCopy, analyticFisher), "matrix")    
//       }
//       
//       
//       // For every new parameter:
//       curMax = max(x@parameters$par)
//         
// ##// NEW
//         curInds = seq_len(curMax)
//           newInds = curMax + seq_len(max(modCopy@parameters$par) - curMax)
//           V =  H[newInds,newInds] - H[newInds,curInds,drop=FALSE] %*% solve_symmetric(H[curInds,curInds]) %*% H[curInds,newInds,drop=FALSE]
//         V.diag = diag(V)
//           idx = which(V.diag < sqrt(.Machine$double.eps))
//           if(length(idx) > 0L) {
//             V.diag[idx] = as.numeric(NA)
//           }
//           // How many in total?
//           nTotalPars = length(c(curInds,newInds))
//             
//             // 
//             // 
//             // // Effective N:
//             // if (nrow(x@sample@groups) > 1){
//             //   par = modCopy@parameters
//             //   // par$id[!modCopy@parameters$identified] = seq_len(sum(!modCopy@parameters$identified))
//             //   par = par %>% left_join(modCopy@sample@groups, by = c("group_id" = "id")) %>%
//             //     group_by(par) %>% summarize(Neff = sum(nobs))
//             //   Neff = numeric(max(par$par))
//             //   Neff[par$par[par$par!=0]] = par$Neff[par$par!=0]
//             // } else {
//             //   Neff = rep(x@sample@groups$nobs[1],nTotalPars)
//             // }
//             
//             
//             // MIs:
//             // All MIs:
//             mi = numeric(nTotalPars)
//             
//             // mi[newInds] = ifelse(abs(V.diag) < 1e-10,0,((-Neff[newInds]*g[newInds])^2)/V.diag)
//             mi[newInds] = ifelse(abs(V.diag) < 1e-10,0,((-nTotal*g[newInds])^2)/V.diag)
//             if (length(curInds) > 0){
//               // mi[curInds] = ((-Neff[curInds]*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
//               mi[curInds] = ((-nTotal*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
//             }
//             p = pchisq(mi,df = 1,lower.tail = FALSE)     
//               
//               // Compute epc:
//               // d = 0.5 * (-1 * Neff) * g
//               d = 0.5 * (-1 * nTotal) * g
//               // needed? probably not; just in case
//               d[which(abs(d) < 1e-15)] = 1.0
//             
//             // Expected parameter change:
//             epc =   mi / d 
//               
//               // Which to fill:
//               fillInds = match(c(curInds,newInds),modCopy@parameters$par)
//               if (type == "normal"){
//                 x@parameters$mi[fillInds[!is.na(fillInds)]] = round(mi[!is.na(fillInds)],10) // round(mi, 3)
//                 x@parameters$pmi[fillInds[!is.na(fillInds)]] = round(p[!is.na(fillInds)],10)
//                 x@parameters$epc[fillInds[!is.na(fillInds)]] = round(epc[!is.na(fillInds)],10)
//               } else if (type == "free"){
//                 x@parameters$mi_free[fillInds[!is.na(fillInds)]] = round(mi[!is.na(fillInds)],10) // round(mi, 3)
//                 x@parameters$pmi_free[fillInds[!is.na(fillInds)]] = round(p[!is.na(fillInds)],10)
//                 x@parameters$epc[fillInds[!is.na(fillInds)]] = round(epc[!is.na(fillInds)],10)
//               } else {
//                 x@parameters$mi_equal[fillInds[!is.na(fillInds)]] = round(mi[!is.na(fillInds)],10) // round(mi,3)
//                 x@parameters$pmi_equal[fillInds[!is.na(fillInds)]] = round(p[!is.na(fillInds)], 10)
//                 x@parameters$epc[fillInds[!is.na(fillInds)]] = round(epc[!is.na(fillInds)],10)
//               }
//               
//               return(x)
//                 
// }
