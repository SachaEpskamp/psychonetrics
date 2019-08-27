// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include "10_ordinal_cppFuns.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// General function to do the following:
// 1. Compute means for continuous variables, thresholds for ordinal
// 2. Compute Pearson/polychoric/polyserial correlations/covariances
// 3. Compute WLS weights matrix
// [[Rcpp::export]]
List covPrepare_cpp(
    NumericMatrix data, // Data as data frame
    LogicalVector isOrdered,
    double tol = 0.000001,
    bool WLSweights = true
) { 
  // Iterators:
  int i, j, k, p;
  int nCase = data.nrow();
  int nUsed;
  double D1, D2;
  
  // Number of variables:
  int nVar = data.ncol();
  if (isOrdered.length() != nVar){
    Rf_error("'isOrdered' is not of the same length as the number of variables.");
  }
  
  // List that stores data as vectors:
  List DataList(nVar);
  
  // List that stores the missingness:
  List MissingList(nVar);
  
  // Integer that stores how many means or thresholds are stored per variable:
  IntegerVector nMeans_Thresh(nVar);
  
  // Fill the list (copying contiuous vectors as I don't want to mess with the input data):
  for (i=0;i<nVar;i++){
    // Dummy numeric:
    NumericVector curVar(nCase);
    
    // Fill the numeric with values:
    for (p=0;p<nCase;p++){
      curVar[p] = data(p,i);
    }
    
    
    // Dummy to store if missing:
    LogicalVector curNA = is_na(curVar);
    MissingList[i] = curNA;
    
    
    
    if (isOrdered[i]){
      // // All unique values sorted
      // NumericVector uniqueVals = unique(na_omit(curVar));
      // 
      // // Sort these:
      // std::sort(uniqueVals.begin(), uniqueVals.end());
      // 
      // // Get the size:
      // int nLevelCur = uniqueVals.size();
      // 
      // // Integervector to store data:
      // IntegerVector ordinalVar(nCase);
      // 
      // // Fill per person:
      // for (p=0;p<nCase;p++){
      //   if (curNA[p] == true){
      //     ordinalVar[p] = NA_INTEGER;
      //   } else {
      //     for (j=0;j<nLevelCur;j++){
      //       if (curVar[p] == uniqueVals[j]){
      //         ordinalVar[p] = j;
      //       }
      //     } 
      //   }
      // }
      // 
      // // Store in data list:
      // DataList[i] = ordinalVar;
      
      DataList[i] = toOrdinal(curVar);
      int nLevelCur = maxInt(DataList[i]) + 1;
      
      // Number of thresholds:
      nMeans_Thresh[i] = nLevelCur - 1;
    } else {
      // Store in data list:
      DataList[i] = curVar;
      
      // Number of means (1):
      nMeans_Thresh[i] = 1;
    }
  }
  
  // List that wil store a mean for every continuous variable or thresholds for every ordered variable:
  List meansAndThresholds(nVar);
  
  // FIXME: Maybe put this in the same loop as above?
  // Fill the list:
  for (i=0;i<nVar;i++){
    if (isOrdered[i] == false){
      // Compute mean:
      meansAndThresholds[i] = computeMean(DataList[i]);
      
    } else {
      // Compute thresholds:
      meansAndThresholds[i] = computeThresholds(DataList[i]);
      
    }
  }
  
  // Covariance matrix:
  NumericMatrix covMat(nVar, nVar);
  
  // Fill lower triangle columwise:
  for (j=0;j<nVar;j++){
    for(i=j;i<nVar;i++){
      // Variance?
      if (i == j){
        if (isOrdered[i]){
          covMat(i,j) = 1.0;
        } else {
          covMat(i,j) = pearsonCov(DataList[i], DataList[j], meansAndThresholds[i], meansAndThresholds[j]);
        }
      } else {
        // Else covariance or correlation:  
        if (isOrdered[i] && isOrdered[j]){
          // Rf_error("Polychoric correlation not yet supported");
          covMat(i,j) = covMat(j,i) = estimate_polychoric(DataList[i], DataList[j],  meansAndThresholds[i],  meansAndThresholds[j], tol = tol);
          
          
        } else if (isOrdered[i] && !isOrdered[j]){
          Rf_error("Polyserial correlation not yet supported");
        } else if (!isOrdered[i] && isOrdered[j]){
          Rf_error("Polyserial correlation not yet supported");
        } else {
          covMat(i,j) = covMat(j,i) = pearsonCov(DataList[i], DataList[j], meansAndThresholds[i], meansAndThresholds[j]);
        }
      }
    }
  }
  
  // Output list:
  List Result;
  
  // Add the raw data (probably should be removed):
  Result["data"] = DataList;
  
  // Add missings (probably should be removed):
  Result["missing"] = MissingList;
  
  // Add means and thresholds:
  Result["means_thresholds"] = meansAndThresholds;
  
  // Add var-cov matrix:
  Result["covmat"] = covMat;
  
  int nElements = 0;
  
  if (WLSweights){
    // Make the parameter vector { mean / threhsolds ; lower tri covariances }
    // First let's count how many elements I need:
    
    for (j=0;j<nVar;j++){
      
      if (isOrdered[j]){
        // Add thresholds
        nElements += nMeans_Thresh[j];
      } else {
        // Add one mean:
        nElements += 1;
      }
      // 
      // Var/cov:
      for(i=j;i<nVar;i++){
        if (i == j && !isOrdered[i]){
          nElements++;
        }
        
        if (i != j){
          nElements++;
        }
      }
    }
    
    // Parameter vector:
    NumericVector parVector(nElements, 0.0);
    
    
    // Integer vector to store which variable belongs to each element. -1 in var2 indicates a mean or threshold:
    IntegerVector var1(nElements);
    IntegerVector var2(nElements);
    // And a vector storing which threshold an element indicates:
    IntegerVector whichPar(nElements, 0);
    
    int curPar = 0;
    for (i=0;i<nVar;i++){
      for (j=0;j<nMeans_Thresh[i];j++){
        whichPar[curPar] = j;
        var1[curPar] = i;
        var2[curPar] = -1;
        curPar++;
      }
    }
    for (j=0;j<nVar;j++){
      for(i=j;i<nVar;i++){
        if (i == j && !isOrdered[i]){
          var1[curPar] = j;
          var2[curPar] = i;
          curPar++;
        }
        
        if (i != j){
          var1[curPar] = j;
          var2[curPar] = i;
          curPar++;
        }
      }
    }
    
    // Following Muthen 1984, I need two matrices:
    // DD: An outer product of the derivative vectors
    // B: A block matrix with elements of DD
    arma::mat B(nElements, nElements); // <- FIXME: I know B is SPARSE, but can't get it's inverse...
    arma::mat DD(nElements, nElements);
    
    bool var1mt = true;
    bool var2mt = true;
    
    
    // Fill the matrices:
    for (i=0;i<nElements;i++){
      for (j=0;j<=i;j++){
        var1mt = var2[i] == -1;
        var2mt = var2[j] == -1;
        DD(i,j) = DD(j,i) = B(i,j) = B(j,i) = 0.0;
        
        // Both a mean or threshold?
        if (var1mt && var2mt){
          // Missing patterns:
          LogicalVector mis1 = MissingList[var1[i]];
          LogicalVector mis2 = MissingList[var1[j]];
          nUsed = sum(!mis1 * !mis2);
          
          
          // For every subject:
          for (p=0; p<nCase; p++){
            
            if (!mis1[p] && !mis2[p]){
              
              // Var 1:
              if (isOrdered[var1[i]]){
                D1 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[i]])[p], whichPar[i], meansAndThresholds[var1[i]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              // Var 2:
              if (isOrdered[var1[j]]){
                D2 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[j]])[p], whichPar[j], meansAndThresholds[var1[j]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              
              // Fill in matrices:
              DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1* (-2.0/(double)nUsed)  * D2;
              // DD(i,j) = DD(j,i) = DD(i,j) +  D1 *  D2;
              if (var1[i] == var1[j]){
                B(i,j) = B(j,i) = DD(i,j);
              }
            }
          }
        }
        
        
        // part 1 is (co)variance and part 2 mean/threshold?
        if (!var1mt && var2mt){
          // Missing patterns:
          LogicalVector mis1a = MissingList[var1[i]];
          LogicalVector mis1b = MissingList[var2[i]];
          LogicalVector mis2 = MissingList[var1[j]];
          
          nUsed = sum(!mis1a * !mis1b * !mis2);
          
          
          // For every subject:
          for (p=0; p<nCase; p++){
            
            if (!mis1a[p] * !mis1b[p] * !mis2[p]){
              
              // part 1:
              if (isOrdered[var1[i]] && isOrdered[var2[i]]){
                D1 = polychor_grad_singlesubject(
                  ((IntegerVector)DataList[var1[i]])[p], 
                                                    ((IntegerVector)DataList[var2[i]])[p],
                                                                                      covMat(var1[i],var2[i]),
                                                                                      meansAndThresholds[var1[i]],
                                                                                                        meansAndThresholds[var2[i]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              // part 2:
              if (isOrdered[var1[j]]){
                D2 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[j]])[p], whichPar[j], meansAndThresholds[var1[j]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              
              // Fill in matrices:
              // DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
              // if (var1[i] == var1[j] || var1[i] == var2[j]){
              //   B(i,j) = B(j,i) = DD(i,j);
              // }
              // Not summetrical B?
              DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
              // DD(i,j) = DD(j,i) = DD(i,j) +  D1 *  D2;
              if (var1[i] == var1[j] || var1[i] == var2[j]){
                B(i,j) = DD(i,j);
              }
            }
          }
        }
        
        // Both a (co)variance?
        if (!var1mt && !var2mt){
          // Missing patterns:
          LogicalVector mis1a = MissingList[var1[i]];
          LogicalVector mis1b = MissingList[var2[i]];
          LogicalVector mis2a = MissingList[var1[j]];
          LogicalVector mis2b = MissingList[var2[j]];
          
          nUsed = sum(!mis1a * !mis1b * !mis2a * !mis2b);
          
          
          // For every subject:
          for (p=0; p<nCase; p++){
            
            if (!mis1a[p] * !mis1b[p] * !mis2a[p] * !mis2b[p]){
              
              // part 1:
              if (isOrdered[var1[i]] && isOrdered[var2[i]]){
                D1 = polychor_grad_singlesubject(
                  ((IntegerVector)DataList[var1[i]])[p], 
                  ((IntegerVector)DataList[var2[i]])[p],
                  covMat(var1[i],var2[i]),
                 meansAndThresholds[var1[i]],
                  meansAndThresholds[var2[i]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              // part 2:
              if (isOrdered[var1[j]] && isOrdered[var2[j]]){
                D2 = polychor_grad_singlesubject(
                  ((IntegerVector)DataList[var1[j]])[p], 
                  ((IntegerVector)DataList[var2[j]])[p],
                 covMat(var1[j],var2[j]),
                 meansAndThresholds[var1[j]],
                 meansAndThresholds[var2[j]]);
              } else {
                Rf_error("Only ordinal data supported now...");
              }
              
              
              // Fill in matrices:
              DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
              // DD(i,j) = DD(j,i) = DD(i,j) +  D1 *  D2;
              if (i == j){
                B(i,j) = B(j,i) = DD(i,j);
              }
            }
          }
        }
        
        // Dummy for now:
        // if (i == j && var2[i] > -1 && var2[j] > -1){
        //   DD(i,j) = 1.0;
        //   B(i,j) = 1.0;
        // }
      }
    }
    
    
    // Rf_PrintValue(wrap(B));
    // Rf_PrintValue(wrap(DD));
    
    // Compute the final matrix:
    arma::mat Binv = inv(B);
    // sp_mat I = diag_ones(nElements);
    // arma::sp_mat Binv = spsolve( B, I );
    arma::mat WLS_V = Binv * DD * Binv;
    
    
    // Store in output:
    Result["parameter_vector"] = parVector;
    Result["parameter_index"] = whichPar;
    Result["pars_var1"] = var1;
    Result["pars_var2"] = var2;
    Result["WLS_V"] = WLS_V;
    Result["DD"] = DD;
    Result["B"] = B;
  }
  
  // Return output:
  return(Result);
}
