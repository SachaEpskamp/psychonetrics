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
    double tol = 0.000001
) { 
  // Iterators:
  int i, j, p;
  int nCase = data.nrow();
  
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
  
  // Return output:
  return(Result);
}
