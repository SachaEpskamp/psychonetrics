// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Asymptotic cov matrix
// [[Rcpp::export]]
arma::sp_mat DWLS_wmat(
    arma::mat data,
    arma::vec means,
    int ncase,
    int nvar) {
  int i, j, g, h, p;
  
  // Asyptotic 2nd, 3rd and 4th order cov matrix:
  double secondorder[nvar][nvar];
  // double thirdorder[nvar][nvar][nvar];
  double fourthorder[nvar][nvar][nvar][nvar];
  for (p=0; p < ncase; p++){
    for (g = 0; g < nvar; g++){
      for (h = g; h < nvar; h ++){
        if (p==0){
          secondorder[g][h] = 0;
        }
        secondorder[g][h] += pow(ncase,-1.0) * (data(p,g) - means(g)) * (data(p,h) - means(h));
        for (i = 0; i < nvar; i++){
          // if (p==0){
          //   thirdorder[g][h][i] = 0;
          // }
          // thirdorder[g][h][i] += pow(ncase,-1.0) * (data(p,g) - means(g)) * (data(p,h) - means(h)) * (data(p,i) - means(i)) ;
          for (j = i; j < nvar; j ++){
            if (p==0){
              fourthorder[i][j][g][h] = 0;
            }
            fourthorder[i][j][g][h] +=  pow(ncase,-1.0) *(data(p,i) - means(i)) * (data(p,j) - means(j)) * (data(p,g) - means(g)) * (data(p,h) - means(h));
          }
        }
      }
    }
  }
  
  // second-order cov mat
  // double S[nvar][nvar];
  // for (i = 0; i < nvar; i++){
  //   for (j = i; j < nvar; j ++){
  //     S[i][j] = 0;
  //     for (p=0; p < ncase; p++){
  //       S[i][j] += (data(p,i) - means(i)) * (data(p,j) - means(j));
  //     }
  //     S[i][j] *= pow(ncase,-1.0);
  //     S[j][i] = S[i][j];
  //   }
  // }
  
  
  // Weights matrix:
  arma::sp_mat Wmat(nvar + nvar*(nvar+1)/2,nvar + nvar*(nvar+1)/2);
  // First fill the mean part:
  // for (i=0;i<nvar;i++){
  //   for (j=0;j<=i;j++){
  //     Wmat(i,j) = Wmat(j,i) = secondorder[j][i];
  //   }
  // }
  
  // Fill the mean by variance part:
  
  
  // Fill the variance part:
  int rowind = nvar;
  int colind = nvar;
  for (g = 0; g < nvar; g++){
    for (h = g; h < nvar; h ++){
      // Fill the mean part:
      if (g == h){
        Wmat(g,h) = Wmat(h,g) = secondorder[g][h]; 
      }
      
      for (i = 0; i < nvar; i++){
        // // Fill the mean by var part:
        // Wmat(i,colind) = Wmat(colind,i) = thirdorder[g][h][i]; // - (means(i) * secondorder[g][h]);
        // 
        for (j = i; j < nvar; j ++){
          // Fill the var by var part:
          if (rowind == colind){
            Wmat(rowind,colind) =  Wmat(colind,rowind) = fourthorder[i][j][g][h] - (secondorder[i][j] * secondorder[g][h]);
          }
          rowind++;
        }
      }
      colind++;
      rowind = nvar;
    }
  }
  
  
  // Return
  return Wmat;
}
