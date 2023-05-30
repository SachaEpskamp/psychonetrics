// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// sparseordense
// [[Rcpp::export]]
bool is_sparse_cpp(
  const arma::mat& X
){
  int i, j;
  
  // Size of X:
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  int nNotZero = 0;
  int nZero = 0;
  
  // Threshold to determine matrix is sparse:
  int thresh_nonZero = 0.25 * nrow * ncol;
  int thresh_Zero = 0.75 * nrow * ncol;
  
  // Result value:
  bool sparse = true;
  
  // Look over all values:
  for (i=0; i<nrow && sparse; i++){
    for (j=0; i<ncol && sparse; i++){
      if (X(i,j) != 0){
        nNotZero++;
      } else {
        nZero++;
      }
      if (nNotZero >= thresh_nonZero){
        sparse = false;
      } else if (nZero >= thresh_Zero){
        sparse = false;
      }
    }
  }

  return(sparse);
}

// diag_sparse_dense: returns 0 for diag, 1 for sparse, 2 for dense
// [[Rcpp::export]]
int diag_sparse_dense_cpp(
    const arma::mat& X
){
  int i, j;
  
  // Size of X:
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  int nNotZero = 0;
  int nZero = 0;
  
  // Threshold to determine matrix is sparse:
  int thresh_nonZero = 0.25 * nrow * ncol;
  int thresh_Zero = 0.75 * nrow * ncol;

  // Result value:
  bool sparse = true;
  bool diag = nrow == ncol;
  
  // Look over all values:
  // Run while we think the matrix might be sparse or diagonal
  for (i=0; i<nrow && (sparse||diag); i++){
    for (j=0; j<ncol && (sparse||diag); j++){
      
      // Is the element nonzero?
      if (X(i,j) != 0){
        
        // Add to counter:
        nNotZero++;
        
        // If this is an offdiagonal element then we know matrix is not diagonal:
        if (i != j){
          diag = false;
        }
      } else {
        
        // Add to zeroes counter:
        nZero++;
      }
      
      // if we can determine matrix is not sparse we can stop:
      if (nNotZero >= thresh_nonZero){
        sparse = false;
      } else if (nZero >= thresh_Zero){
        sparse = false;
      }
    }
  }
  
  // Return value:
  int retval;
  if (diag){
    retval = 0;
  } else if (sparse){
    retval = 1;
  } else {
    retval = 2;
  }
  
  return(retval);
}
