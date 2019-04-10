// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Asymptotic cov matrix
// [[Rcpp::export]]
arma::mat WLS_wmat(
    arma::mat data,
    arma::vec means,
    int ncase,
    int nvar) {
  int i, j, g, h, p;
  
  // Asyptotic cov matrix:
  double acov[nvar][nvar][nvar][nvar];
  
  // Fill acov:
  for (i = 0; i < nvar; i++){
    for (j = i; j < nvar; j ++){
      for (g = 0; g < nvar; g++){
        for (h = g; h < nvar; h ++){
          
          acov[i][j][g][h] = 0;
          for (p=0; p < ncase; p++){
            acov[i][j][g][h] += (data(p,i) - means(i)) * (data(p,j) - means(j)) * (data(p,g) - means(g)) * (data(p,h) - means(h));
          }
          acov[i][j][g][h] *= 1/ncase;
        }
      }
    }
  }
  
  // Weights matrix:
  arma::mat Wmat = zeros(nvar*nvar,nvar*nvar);

  
  // Return
  return Wmat;
}
