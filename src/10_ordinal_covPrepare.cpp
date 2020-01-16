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
    bool WLSweights = true,
    bool verbose = true
) { 
  // Iterators:
  int i, j, k, p;
  int nCase = data.nrow();
  // int nUsed;
  double D2;
  
  // Number of variables:
  int nVar = data.ncol();
  if (isOrdered.length() != nVar){
    Rf_error("'isOrdered' is not of the same length as the number of variables.");
  }
  
  // List that stores data as vectors:
  List DataList(nVar);
  
  // List that stores the missingness:
  // List MissingList(nVar);
  // Matrix that stores missingness:
  LogicalMatrix Missings(nCase, nVar);
  
  // Integer that stores how many means or thresholds are stored per variable:
  IntegerVector nMeans_Thresh(nVar);
  
  
  // List that wil store a mean for every continuous variable or thresholds for every ordered variable:
  List meansAndThresholds(nVar);
  
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
    Missings(_,i) = curNA;
    
    
    
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
  // }
  

  
  // FIXME: Maybe put this in the same loop as above?
  // Fill the list:
  // for (i=0;i<nVar;i++){
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
          covMat(i,j)  = estimate_polychoric(DataList[i], DataList[j],  meansAndThresholds[i],  meansAndThresholds[j], tol);
          if (i != j){
            covMat(j,i) =  covMat(i,j);            
          }

          
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
  Result["missing"] = Missings;
  
  // Add means and thresholds:
  Result["means_thresholds"] = meansAndThresholds;
  
  // Add var-cov matrix:
  Result["covmat"] = covMat;
  
    
  if (WLSweights){
    int nElements = 0;
    int ntotalMT = sum(nMeans_Thresh);
    
    // Augment thresholds for all gradients already:
    List meansAndThresholds_aug(nVar);
    for (i=0;i<nVar;i++){
      if (isOrdered[i]){
        NumericVector curthresh = meansAndThresholds[i];
        NumericVector augthresh(nMeans_Thresh[i] + 2);
        augthresh[0] = -999;
        augthresh[nMeans_Thresh[i] + 1] = 999;
        for (j=0;j<nMeans_Thresh[i];j++){
          augthresh[j+1] = curthresh[j];
        }
        meansAndThresholds_aug[i] = augthresh;        
      } else {
        meansAndThresholds_aug[i] = meansAndThresholds[i];
      }

    }
    
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
    // NumericVector parVector(nElements, 0.0);
    
    
    
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
          var1[curPar] = i;
          var2[curPar] = j;
          curPar++;
        }
        
        if (i != j){
          var1[curPar] = i;
          var2[curPar] = j;
          curPar++;
        }
      }
    }
    
    // Following Muthen 1984, I need two matrices:
    // DD: An outer product of the derivative vectors
    // B: A block matrix with elements of DD
    arma::mat B(nElements, nElements); // <- FIXME: I know B is SPARSE, but can't get it's inverse...
    arma::mat DD(nElements, nElements);
    arma::mat partA21(nElements-ntotalMT,ntotalMT);
    
    
    // Set to zeroes:
    B.zeros();
    DD.zeros();
    
    // I also want to store the effictive N per element:
    // IntegerMatrix effectiveSampleSize(nElements, nElements);
    // arma::mat effectiveSampleSize(nElements, nElements);
    // effectiveSampleSize.zeros();
    
    // I also need a vector per subject:
    arma::vec obsElement(nElements);
    
    // I need a D vector per subject:
    arma::vec D(nElements);
    
    
    // I also want to compute the bivariate likelihood per subject only once, as this is costly. This matrix will store the bivariate probabilities:
    NumericMatrix bivariatLikelihood(nVar, nVar);
    
    
    // bool var1mt = true;
    // bool var2mt = true;
    bool mis1, mis2;
    
    // Loop over all subjects to fill matrices:
    for (p=0; p<nCase;p++){

      // Reset D:
      D.zeros();
      
      // Reset missings:
      obsElement.zeros();

      
      // Loop over all elemenets:
      for (i=0;i<nElements;i++){

        // Check if the relevant variables are missing, or skip otherwise
        mis1 = Missings(p,var1[i]);
        
        // Check if var 2 is missing (not relevant for mean/threshold)
        mis2 = false;
        if (var2[i]>-1){
          mis2 = Missings(p,var2[i]);
        }

        // If missing, nothing to do, else do stuff:
        if (!mis1 & !mis2){

          // Add to the sample size:
          obsElement[i] = 1;

          
          // Is i a mean or threshold?
          if (var2[i] == -1){
            
            // Check for ordinal (nothing else implemented yet)
          
            if (isOrdered[var1[i]]){
              D[i] = threshold_grad_singlesubject(((IntegerVector)DataList[var1[i]])[p], whichPar[i], meansAndThresholds_aug[var1[i]]);
            } else {
              Rf_error("Only ordinal data supported now...");
            }
          } else {

            // Variance or covariance:
            if (isOrdered[var1[i]] & isOrdered[var2[i]]){
            
              // Compute pi:
              bivariatLikelihood(var1[i],var2[i]) = bivariatLikelihood(var2[i],var1[i]) = 
                ordered_bivariate_likelihood(((IntegerVector)DataList[var1[i]])[p], 
                                             ((IntegerVector)DataList[var2[i]])[p],
                                             covMat(var1[i],var2[i]),
                                              meansAndThresholds_aug[var1[i]],
                                             meansAndThresholds_aug[var2[i]]);
            
              D[i] = polychor_grad_singlesubject(
                ((IntegerVector)DataList[var1[i]])[p], 
                ((IntegerVector)DataList[var2[i]])[p],
                covMat(var1[i],var2[i]),
                meansAndThresholds_aug[var1[i]],
                meansAndThresholds_aug[var2[i]],
                bivariatLikelihood(var1[i],var2[i]));

              
            } else {
              Rf_error("Only ordinal data supported now...");
            }         
          }
        }
      }
      

      // Now that I looped over all elements, I can update DD:
      DD += D * D.t();
      
      // Also update sample size:
      // effectiveSampleSize += obsElement * obsElement.t();
      
     
      
      // TODO: Part A21
      // I need to loop over all means/thresholds:
      // FIXME: Long and ugly loop... sorry...
      for (j=0;j<ntotalMT;j++){
        for (i=ntotalMT;i<nElements;i++){
          if ( !Missings(p,var1[i]) && !Missings(p,var2[i])){
            if (var1[j] == var1[i]){
              D2 = bthreshold_grad_singlesubject(
                ((IntegerVector)DataList[var1[j]])[p],
                ((IntegerVector)DataList[var2[i]])[p],
                covMat(var1[j],var2[i]),
                whichPar[j],
                meansAndThresholds_aug[var1[j]],
                meansAndThresholds_aug[var2[i]],
               bivariatLikelihood(var1[j],var2[i]));
              B(i,j) +=  D[i] * D2; 
            }
            
            if (var1[j] == var2[i]){
              D2 = bthreshold_grad_singlesubject(
                ((IntegerVector)DataList[var1[j]])[p],
                ((IntegerVector)DataList[var1[i]])[p],
               covMat(var1[j],var1[i]),
                whichPar[j],
                 meansAndThresholds_aug[var1[j]],
                  meansAndThresholds_aug[var1[i]],
                  bivariatLikelihood(var1[j],var1[i]));
              
              
              B(i,j) +=  D[i] * D2;
            }
            
          }

          
        }
      }
      
  
    }
    
    // As for B, I can directly copy A11 and A22. First looping over all means and thresholds:
    curPar=0;
    for (i=0;i<nVar;i++){
      for (j=0;j<nMeans_Thresh[i];j++){
        for (k=0;k<nMeans_Thresh[i]-j;k++){
          B(curPar + k, curPar) =  B(curPar, curPar + k) = DD(curPar, curPar+k);
        }
        curPar++;
      }
    }
    
    // Next for all variances/covariances I just need to copy the diagonal::
    for (j=0;j<nVar;j++){
      for(i=j;i<nVar;i++){
        if (i == j && !isOrdered[i]){
          B(curPar, curPar) = DD(curPar, curPar);
          curPar++;
        }
        
        if (i != j){
          B(curPar, curPar) = DD(curPar, curPar);
          curPar++;
        }
      }
    }
    
    
    // Finally, adjust for effective sample size:
    // B = B % pow(effectiveSampleSize, -1.0);
    // DD = DD % pow(effectiveSampleSize, -1.0);
    // 
    // B = B * pow((double)nCase, -1.0);
    // DD = DD * pow((double)nCase, -1.0);
    // 
    // 
    // 
    // // Fill the matrices:
    // for (j=0;j<nElements;j++){
    //   for (i=j;i<nElements;i++){
    //     // Initial setup:
    //     var1mt = var2[i] == -1; // Is i a mean/threshold?
    //     var2mt = var2[j] == -1; // Is j a mean/threshold?
    //     
    //     // Initialize to 0:
    //     DD(i,j) = DD(j,i) = B(i,j) = B(j,i) = 0.0;
    //     effectiveSampleSize(i,j) = effectiveSampleSize(j,i) = 0;
    //     
    //     
    //     // Both a mean or threshold?
    //     if (var1mt && var2mt){
    //       // For every subject:
    //       for (p=0; p<nCase; p++){
    //         
    //         if (!Missings(p,var1[i]) * !Missings(p,var1[j])){
    //           
    //           // Var 1:
    //           if (isOrdered[var1[i]]){
    //             D1 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[i]])[p], whichPar[i], meansAndThresholds[var1[i]]);
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           // Var 2:
    //           if (isOrdered[var1[j]]){
    //             D2 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[j]])[p], whichPar[j], meansAndThresholds[var1[j]]);
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           
    //           // Fill in matrices:
    //           // DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1* (-2.0/(double)nUsed)  * D2;
    //           DD(i,j) = DD(i,j) +  D1 *  D2;
    //           effectiveSampleSize(i,j)++;
    //           
    //           
    //           if (i != j){
    //             DD(j,i) = DD(j,i) +  D1 *  D2;                
    //             effectiveSampleSize(j,i)++;
    //           }
    //         }
    //       }
    //       if (var1[i] == var1[j]){
    //         B(i,j) = DD(i,j);
    //         if (i != j){
    //           B(j,i) = DD(j,i);
    //         }
    //       }
    //     }
    //   
    //     
    //     // part 1 is (co)variance and part 2 mean/threshold?
    //     if (!var1mt && var2mt){
    //       
    //       // For every subject:
    //       for (p=0; p<nCase; p++){
    //         
    //         if (!Missings(p,var1[i]) * !Missings(p,var2[i]) * !Missings(p,var1[j])){
    //           
    //           // part 1:
    //           if (isOrdered[var1[i]] && isOrdered[var2[i]]){
    //             D1 = polychor_grad_singlesubject(
    //               ((IntegerVector)DataList[var1[i]])[p], 
    //               ((IntegerVector)DataList[var2[i]])[p],
    //               covMat(var1[i],var2[i]),
    //               meansAndThresholds[var1[i]],
    //               meansAndThresholds[var2[i]]);
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           // part 2:
    //           if (isOrdered[var1[j]]){
    //             
    //             if (var1[i] == var1[j] || var2[i] == var2[j]){
    //               // For B I need the derivative to the bivariate function:
    //               if (var1[i] == var1[j]){
    //                 D2bpart = bthreshold_grad_singlesubject(
    //                   ((IntegerVector)DataList[var1[j]])[p], 
    //                   ((IntegerVector)DataList[var2[i]])[p],
    //                   covMat(var1[j],var2[i]),
    //                   whichPar[j],
    //                    meansAndThresholds[var1[j]],
    //                    meansAndThresholds[var2[i]]);                    
    //               } else {
    //                 D2bpart = bthreshold_grad_singlesubject(
    //                   ((IntegerVector)DataList[var1[j]])[p], 
    //                     ((IntegerVector)DataList[var1[i]])[p],
    //                     covMat(var1[j],var1[i]),
    //                      whichPar[j],
    //                     meansAndThresholds[var1[j]],
    //                    meansAndThresholds[var1[i]]);   
    //               }
    // 
    //             }
    // 
    //               
    //             // For DD instead I need the one to the univariate function:
    //             D2 = threshold_grad_singlesubject(((IntegerVector)DataList[var1[j]])[p], whichPar[j], meansAndThresholds[var1[j]]);
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           
    //           // Fill in matrices:
    //           // DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
    //           // if (var1[i] == var1[j] || var1[i] == var2[j]){
    //           //   B(i,j) = B(j,i) = DD(i,j);
    //           // }
    //           // Not summetrical B?
    //           // DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
    //           DD(i,j) = DD(i,j) +  D1 *  D2;
    //           effectiveSampleSize(i,j)++;
    //           
    //           if (var1[i] == var1[j] || var2[i] == var1[j]){
    //             B(i,j) = B(i,j) + D1 * D2bpart;
    //           }
    //           
    //           if (i != j){
    //             DD(j,i) = DD(j,i) +  D1 *  D2;                
    //             effectiveSampleSize(j,i)++;
    //           }
    //          
    //         }
    //       }
    //     }
    //     
    //     // Both a (co)variance?
    //     if (!var1mt && !var2mt){
    //       
    //       // For every subject:
    //       for (p=0; p<nCase; p++){
    //         
    //         if (!Missings(p,var1[i]) * !Missings(p,var2[i]) * !Missings(p,var1[j]) * !Missings(p,var2[j])){
    //     
    //           
    //           // part 1:
    //           if (isOrdered[var1[i]] && isOrdered[var2[i]]){
    //             D1 = polychor_grad_singlesubject(
    //               ((IntegerVector)DataList[var1[i]])[p], 
    //             ((IntegerVector)DataList[var2[i]])[p],
    //             covMat(var1[i],var2[i]),
    //                                                                                   meansAndThresholds[var1[i]],
    //                                                                                                     meansAndThresholds[var2[i]]);
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           // part 2:
    //           if (isOrdered[var1[j]] && isOrdered[var2[j]]){
    //             D2 = polychor_grad_singlesubject(
    //               ((IntegerVector)DataList[var1[j]])[p], 
    //               ((IntegerVector)DataList[var2[j]])[p],
    //               covMat(var1[j],var2[j]),
    //               meansAndThresholds[var1[j]],
    //               meansAndThresholds[var2[j]]);
    //             
    //             
    //             // 
    //             // if (std::isnan(D2)){
    //             //   Rf_PrintValue(wrap(((IntegerVector)DataList[var1[j]])[p]));
    //             //   Rf_PrintValue(wrap(((IntegerVector)DataList[var2[j]])[p]));
    //             //   Rf_PrintValue(wrap(covMat(var1[j],var2[j])));
    //             //   Rf_PrintValue(wrap(meansAndThresholds[var1[j]]));
    //             //   Rf_PrintValue(wrap( meansAndThresholds[var2[j]]));
    //             // }
    //           } else {
    //             Rf_error("Only ordinal data supported now...");
    //           }
    //           
    //           
    //           // Fill in matrices:
    //           // DD(i,j) = DD(j,i) = DD(i,j) + (-2.0/(double)nUsed) * D1 * (-2.0/(double)nUsed) * D2;
    //           DD(i,j) = DD(i,j) +  D1 *  D2;
    //           effectiveSampleSize(i,j)++;
    //           
    //           
    //           if (i != j){
    //             DD(j,i) = DD(j,i) +  D1 *  D2;                
    //             effectiveSampleSize(j,i)++;
    //           }
    //         }
    //       }
    //       if (i == j){
    //         B(i,j) = DD(i,j);
    //       }
    //     }
    //     
    //     // Dummy for now:
    //     // if (i == j && var2[i] > -1 && var2[j] > -1){
    //     //   DD(i,j) = 1.0;
    //     //   B(i,j) = 1.0;
    //     // }
    //     // 
    //     // // Multiply by -2/n:
    //     // B(i,j) *= -2.0 / ((double)effectiveSampleSize(i,j)) * -2.0 / ((double)effectiveSampleSize(i,j));
    //     // DD(i,j) *= -2.0 / ((double)effectiveSampleSize(i,j)) * -2.0 / ((double)effectiveSampleSize(i,j));
    //     // 
    //     // if (i != j){
    //     //   B(j,i) *= -2.0 / ((double)effectiveSampleSize(j,i)) * -2.0 / ((double)effectiveSampleSize(j,i));
    //     //   DD(j,i) *= -2.0 / ((double)effectiveSampleSize(j,i)) * -2.0 / ((double)effectiveSampleSize(j,i));
    //     // }
    //     
    //     
    //     // Multiply by 1/n FIXME: Why?!:
    //     B(i,j) *= 1.0 / ((double)effectiveSampleSize(i,j));
    //     DD(i,j) *= 1.0 / ((double)effectiveSampleSize(i,j));
    //     
    //     if (i != j){
    //       B(j,i) *= 1.0 / ((double)effectiveSampleSize(j,i));
    //       DD(j,i) *= 1.0 / ((double)effectiveSampleSize(j,i));
    //     }
    // 
    // 
    //   }
    // }
    // 
    

    // Compute the final matrix:
    arma::mat Binv = inv(B);
    arma::mat WLS_V = pow((double)nCase, -1.0) * inv(Binv * DD * Binv.t());
    
    // Store in output:
    Result["parameter_index"] = whichPar;
    Result["pars_var1"] = var1;
    Result["pars_var2"] = var2;
    Result["WLS_V"] = WLS_V;
    Result["DD"] = DD;
    Result["B"] = B;
    // Result["effectiveSampleSize"] = effectiveSampleSize;
  }
  
  // Return output:
  return(Result);
}
