// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Eigenvalues of symmetric matrix:
// [[Rcpp::export]]
arma::vec eig_sym_cpp(
    arma::mat X
){
  return(eig_sym(arma::symmatl(X)));
}

// Check symmetric pd:
// [[Rcpp::export]]
bool sympd_cpp(
    arma::mat X
){
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // Check if posdef:
  double epsilon = 1.490116e-08;
  
  bool posdef = eig_sym(arma::symmatl(X))[0] > -epsilon; //X.is_sympd();
  
  // return:
  return(posdef);
}
// 
// // Nearest PD matrix:
// // [[Rcpp::export]]
// void psychonetrics_nearPD(
//     arma::mat &X
// ){
//   // Number of variables:
//   int nvar = X.n_cols;
//   
//   // Counter:
//   int i;
//   
//   // Initiate eigen values and vectors:
//   arma::vec eigval;
//   arma::mat eigvec;
//   
//   // Check if symmetric:
//   bool issym = X.is_symmetric();
//   
//   // if not, make symmetric:
//   if (!issym){
//     X = 0.5* (X + X.t());
//   }
//   
//   // Compute:
//   arma::eig_sym(eigval, eigvec, arma::symmatl(X));
//   
//   // FIXME: hardcoded epsilon:
//   double sqrt_epsilon = 1.490116e-08;
//   
//   // If lowest eigenvalue is > sqrt(epsilon), do nothing:
//   if (eigval[0] < sqrt_epsilon){
//     i = 0;
//     while (eigval[i] < 0){
//       eigval[i] = 0;
//       i++;
//     }
//     arma::mat Lambda = diagmat(eigval);
//     X = eigvec * Lambda * eigvec.t();
//     
//     // Now do a small spectral shift:
//     for (i=0;i<nvar;i++){
//       //X(i,i) -= lowestEV;
//       X(i,i) = X(i,i)  + sqrt_epsilon;
//     }
//     
//   } 
// }
// 
// // Symmetric solve:
// // [[Rcpp::export]]
// Rcpp::List solve_symmetric_cpp(
//     arma::mat X,
//     bool logdet = false,
//     double sqrt_epsilon  = 1.490116e-08
// ){
//   
//   
//   double logdetval = R_NegInf;
//   Rcpp::List res;
//   int i;
//   int nvar = X.n_cols;
//   
//   // Check if symmetric:
//   bool issym = X.is_symmetric();
//   
//   // if not, make symmetric:
//   if (!issym){
//     X = 0.5* (X + X.t());
//   }
//   
//   // // Check if posdef:
//   // bool posdef = X.is_sympd();
//   // Check if posdef:
//   double lowestEV = eig_sym(arma::symmatl(X))[0];
//   
//   bool posdef = lowestEV > -sqrt_epsilon; //X.is_sympd();
//   
//   // If not, pseudoinverse:
//   if (!posdef){
//     
//     
//     arma::mat inv = pinv(X);
//     res["inv"] = inv;
//     
//     if (logdet){
//       logdetval = log(sqrt_epsilon);
//       res["logdet"] = logdetval;
//     }
//     
//   } else {
//     
//     
//     if (lowestEV < sqrt_epsilon){
//       // Nearest PD matrix and small spectral shift if needed:
//       psychonetrics_nearPD(X);
//       // 
//       // for (i=0;i<nvar;i++){
//       //   //X(i,i) -= lowestEV;
//       //   X(i,i) = X(i,i) -lowestEV + sqrt(epsilon);
//       // }
//     }
//     // Rf_PrintValue(wrap(lowestEV));
//     // invert:
//     arma::mat inv = inv_sympd(X); // FIXME
//     // arma::mat inv = inv(X);
//     res["inv"] = inv;
//     
//     if (logdet){
//       double logepsilon = log(sqrt_epsilon);
//       logdetval =  log(det(inv));
//       if (logdetval == R_PosInf){
//         logdetval = real(log_det(inv));
//       }
//       if (logdetval < logepsilon){
//         logdetval = logepsilon;
//       }
//       
//       res["logdet"] = logdetval;
//     }
//   }
//   
//   return(res);
// }
// 
// 
// // [[Rcpp::export]]
// arma::mat solve_symmetric_cpp_matrixonly(
//     arma::mat X,
//     double sqrt_epsilon  = 1.490116e-08
// ){
//   int i;
//   int nvar = X.n_cols;
//   
//   // Check if symmetric:
//   bool issym = X.is_symmetric();
//   
//   // if not, make symmetric:
//   if (!issym){
//     X = 0.5* (X + X.t());
//   }
//   
//   // // Check if posdef:
//   // bool posdef = X.is_sympd();
//   // Check if posdef:
//   // Rf_PrintValue(wrap(X));
//   
//   double lowestEV = eig_sym(arma::symmatl(X))[0];
//   
//   
//   bool posdef = lowestEV > - sqrt_epsilon; //X.is_sympd();
//   
//   
//   // If not, pseudoinverse:
//   if (!posdef){
//     
//     
//     arma::mat res = pinv(X);
//     return(res);
//     
//   } else {
//     // Small spectral shift:
//     if (lowestEV < sqrt_epsilon){
//       // Nearest PD matrix and small spectral shift if needed:
//       psychonetrics_nearPD(X);
//       
//     }
//     
//     // invert:
//     arma::mat res = inv_sympd(X); // FIXME
//     return(res);
//   }
//   
//   
// }
// 
// 
// 
// // Same as above, but with extra check:
// 
// // [[Rcpp::export]]
// arma::mat solve_symmetric_cpp_matrixonly_withcheck(
//     arma::mat X,
//     bool& proper
// ){
//   double sqrt_epsilon  = 1.490116e-08;
//   int i,j;
//   int nvar = X.n_cols;
//   
//   // Check if symmetric:
//   bool issym = X.is_symmetric();
//   
//   // if not, make symmetric:
//   if (!issym){
//     X = 0.5* (X + X.t());
//   }
//   
//   // Check again:
//   issym = X.is_symmetric();
//   
//   // if not, loop over to fix:
//   if (!issym){
//     proper = false;
//     for (i=0;i<nvar;i++){
//       for (j=0;j<nvar;j++){
//         if (!is_finite(X(i,j))){
//           if (i==j){
//             X(i,j) = 1;
//           } else {
//             X(i,j) = 0;
//           }
//         }
//         
//       }
//     }
//   }
//   
//   // // Check if posdef:
//   // bool posdef = X.is_sympd();
//   // Check if posdef:
//   // Rf_PrintValue(wrap(X));
//   
//   
//   double lowestEV = eig_sym(arma::symmatl(X))[0];
//   
//   
//   bool posdef = lowestEV > - sqrt_epsilon; //X.is_sympd();
//   
//   // If not, pseudoinverse:
//   if (!posdef){
//     proper = false;
//     
//     arma::mat res = pinv(X);
//     return(res);
//     
//   } else {
//     // Small spectral shift:
//     // if (lowestEV < sqrt(epsilon)){
//     //   proper = false;
//     //   
//     //   for (i=0;i<nvar;i++){
//     //     X(i,i) = X(i,i) - lowestEV + sqrt(epsilon);
//     //   }
//     if (lowestEV < sqrt_epsilon){
//       
//       proper = false;
//       
//       // Nearest PD matrix and small spectral shift if needed:
//       psychonetrics_nearPD(X);
//       
//     }
//     
//     // invert:
//     arma::mat res = inv_sympd(X); // FIXME
//     return(res);
//   }
//   
//   
// }


// New ways of inverting matrices (0.8.2+):

// Symmetric solve:
// [[Rcpp::export]]
Rcpp::List solve_symmetric_cpp(
    arma::mat X,
    bool logdet = false,
    double sqrt_epsilon  = 1.490116e-08
){
  

  
  double logdetval = R_NegInf;
  Rcpp::List res;
  int i;
  int nvar = X.n_cols;
  
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // Return value matrix:
  mat inv(nvar,nvar);
  
  // First try a normal inverse:
  bool success = inv_sympd(inv,X);
  
  // If this failed, do pseudo-inverse:
  if (!success){
    
    inv = pinv(X);
    res["inv"] = inv;
    
    if (logdet){
      logdetval = log(sqrt_epsilon);
      res["logdet"] = logdetval;
    }
    
  } else {
    
    res["inv"] = inv;
    
    if (logdet){
      double logepsilon = log(sqrt_epsilon);
      logdetval =  log(det(inv));
      if (logdetval == R_PosInf){
        logdetval = real(log_det(inv));
      }
      if (logdetval < logepsilon){
        logdetval = logepsilon;
      }
      
      res["logdet"] = logdetval;
    }
  }
  
  return(res);
}


// [[Rcpp::export]]
arma::mat solve_symmetric_cpp_matrixonly(
    arma::mat X,
    double sqrt_epsilon  = 1.490116e-08
){
  int i;
  int nvar = X.n_cols;
  
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // Return value matrix:
  mat inv(nvar,nvar);
  
  // First try a normal inverse:
  bool success = inv_sympd(inv,X);
  
  // If this failed, do pseudo-inverse:
  if (!success){
    
    inv = pinv(X);
    
  }
  
  return(inv);
}



// Same as above, but with extra check:

// [[Rcpp::export]]
arma::mat solve_symmetric_cpp_matrixonly_withcheck(
    arma::mat X,
    bool& proper
){
  double sqrt_epsilon  = 1.490116e-08;
  int i,j;
  int nvar = X.n_cols;
  
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // Check again:
  issym = X.is_symmetric();
  
  // if not, loop over to fix:
  if (!issym){
    proper = false;
    for (i=0;i<nvar;i++){
      for (j=0;j<nvar;j++){
        if (!is_finite(X(i,j))){
          if (i==j){
            X(i,j) = 1;
          } else {
            X(i,j) = 0;
          }
        }
        
      }
    }
  }
  
  // Return value matrix:
  mat inv(nvar,nvar);
  
  // First try a normal inverse:
  bool success = inv_sympd(inv,X);
  
  // If this failed, do pseudo-inverse:
  if (!success){
    
    inv = pinv(X);
    
  }
  
  // proper:
  proper = !success;
  
  return(inv);
  
  
}


// Block diag:
// [[Rcpp::export]]
arma::mat bdiag_psychonetrics(
    const Rcpp::List  mats
){
  int nMat = mats.length();
  arma::vec cols(nMat);
  arma::vec rows(nMat);
  int totalrow = 0;
  int totalcol = 0;
  int i, j, k;
  
  // Compute size:
  for (i=0;i<nMat;i++){
    arma::mat curmat = mats[i];
    cols[i] = curmat.n_cols;
    totalcol += cols[i];
    rows[i] = curmat.n_rows;
    totalrow += rows[i];
  }
  
  // Form zero matrix:
  arma::mat diagmat = zeros(totalrow, totalcol);
  
  // Fill the matrix:
  int startrow = 0;
  int startcol = 0;
  
  for (k =0; k<nMat; k++){
    arma::mat curmat = mats[k];
    
    for (i=0;i<rows[k];i++){
      
      for (j=0;j<cols[k];j++){
        
        diagmat(startrow+i, startcol + j)  = curmat(i,j);
        
      }
    }
    startrow += rows[k];
    startcol += cols[k];
  }
  
  
  return(diagmat);
}

// [[Rcpp::export]]
arma::mat cbind_psychonetrics(
    const Rcpp::List  mats
){
  int nMat = mats.length();
  arma::vec cols(nMat);
  arma::vec rows(nMat);
  
  arma::mat curmat = mats[0];
  
  int totalrow = curmat.n_rows;
  int totalcol = 0;
  int i, j, k;
  
  // Compute size:
  for (i=0;i<nMat;i++){
    arma::mat curmat = mats[i];
    cols[i] = curmat.n_cols;
    totalcol += cols[i];
    rows[i] = curmat.n_rows;
    if (rows[i] != totalrow){
      Rf_error("Number of rows are not consistent");
    }
    // totalrow += rows[i];
  }
  
  // Form zero matrix:
  arma::mat resmat = zeros(totalrow, totalcol);
  
  // Fill the matrix:
  // int startrow = 0;
  int startcol = 0;
  
  for (k =0; k<nMat; k++){
    arma::mat curmat = mats[k];
    
    for (i=0;i<rows[k];i++){
      
      for (j=0;j<cols[k];j++){
        
        resmat(i, startcol + j)  = curmat(i,j);
        
      }
    }
    // startrow += rows[k];
    startcol += cols[k];
  }
  
  
  return(resmat);
}

// Half vectorization (lower triangle):
// [[Rcpp::export]]
arma::vec vech(
    arma::mat X,
    bool diag = true
){
  int nvar = X.n_rows;
  int nelement = nvar * (nvar - 1 + 2 * diag) / 2;
  int curel = 0;
  arma::vec out(nelement);
  for (int j=0;j<nvar;j++){
    for (int i=j;i<nvar;i++){
      
      if (diag || (i != j)){
        
        // if (diag){
        out(curel) = X(i,j);  
        curel++;
        // }
        //   
        // } else {
        //   out(curel) = X(i,j);
        //   curel++;
        // }
        
      }
    }
  }
  
  
  return(out);
}


// Indices (start end) for seq_len - 1:
// [[Rcpp::export]]
arma::vec seq_len_inds(
    int start,
    int n
){
  arma::vec res(2);
  res[0] = start;
  res[1] = start + n - 1;
  return(res);
}


// [[Rcpp::export]]
arma::mat cov2cor_cpp(
    arma::mat X
){
  int i, j;
  int n = X.n_rows;
  
  // // Check diagonal:
  // for (i=0;i<n;i++){
  //   if (X(i,i)){
  //     if (X(i,i) < 1.490116e-08){
  //       X(i,i) = 1.490116e-08;
  //     }
  //   }
  // }
  
  
  arma::mat Y = eye(n,n);
  
  
  for (i=0;i<n;i++){
    for (j=0;j<=i;j++){
      Y(i,j) = Y(j,i) = 
        X(i,j) / sqrt(X(i,i) * X(j,j));
    }  
  }
  
  
  return(Y);
}

// [[Rcpp::export]]
arma::mat wi2net_cpp(
    const arma::mat& X
){
  int n = X.n_rows;
  arma::mat Y = zeros(n,n);
  
  int i, j;
  for (i=0;i<n;i++){
    for (j=0;j<i;j++){
      Y(i,j) = Y(j,i) = 
        - X(i,j) / sqrt(X(i,i) * X(j,j));
    }  
  }
  
  return(Y);
}

// [[Rcpp::export]]
arma::mat SDmat(
    const arma::mat& X
){
  int n = X.n_rows;
  arma::mat Y = zeros(n,n);
  
  int i;
  for (i=0;i<n;i++){
    Y(i,i) = sqrt(X(i,i));
  }
  
  return(Y);
}

// [[Rcpp::export]]
arma::mat invSDmat(
    const arma::mat& X
){
  int n = X.n_rows;
  arma::mat Y = zeros(n,n);
  
  int i;
  for (i=0;i<n;i++){
    if (X(i,i) > 0){
      Y(i,i) = pow(X(i,i), -0.5);  
    }
    
  }
  
  return(Y);
}

// [[Rcpp::export]]
bool anyNon0(
    const arma::mat& X
){
  bool anyNon0 = false;
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  for (int i=0; i<nrow; i++){
    for (int j=0; j<ncol; j++){
      
      if (X(i,j) != 0){
        anyNon0 = true;
      }
      
    }
  }
  
  return(anyNon0);
}


// [[Rcpp::export]]
void growlist(
    Rcpp::List& X,
    const Rcpp::List Y
){
  // Names of Y:
  CharacterVector names = Y.names();
  
  // Length of Y:
  int n = Y.length();
  
  // var:
  std::string var;
  
  // Loop
  for (int i=0; i < n; i++){
    var = names[i];
    X[var] = Y[i];
  }
  
}

// [[Rcpp::export]]
arma::vec parVector_cpp(
    const S4& model
){
  Rcpp::List pars = model.slot("parameters");
  arma::vec parnum = pars["par"];
  arma::vec parest = pars["est"];
  
  // Number of free parameters:
  int freePar = max(parnum);
  
  // Number of total parameters:
  int totalPar = parnum.n_elem;
  
  // output:
  arma::vec parvec(freePar);
  
  for (int i=0;i<totalPar;i++){
    if (parnum(i)>0){
      parvec(parnum(i)-1) = parest(i);
    }
  }
  
  return(parvec);
}

// [[Rcpp::export]]
arma::mat computePDC_cpp(
    const arma::mat& beta,
    const arma::mat& kappa,
    const arma::mat& sigma
){
  
  arma::vec sigmaDiag = sigma.diag();
  arma::vec kappaDiag = kappa.diag();
  arma::mat PDCt = beta / sqrt(sigmaDiag * kappaDiag.t() + beta % beta);
  return(PDCt.t());
}


// [[Rcpp::export]]
arma::mat blockToeplitz_cpp(
    const Rcpp::List& X
){
  int b,r, start_col, start_row, end_col, end_row;
  
  // length:
  int nBlock = X.length();
  
  // first block:
  arma::mat first = X[0];
  int n = first.n_rows; 
  
  
  // Assume symmetric and all equal
  arma::mat Toeplitz = zeros(n*nBlock,n*nBlock);
  
  // Fill for each block:
  for (b=0;b<nBlock;b++){
    
    // Block 1 (0 in cpp) is repeated n times, block  n (n-1 in cpp) is repeated 1 time...
    for (r=0;r<(nBlock-b);r++){
      arma::mat block = X[b];
      
      // The r-th repition starts on the r*nth column and the 
      start_row = r * n + b * n;
      start_col = start_row - b*n;
      end_row = start_row + n - 1;
      end_col = start_col + n - 1;
      
      Toeplitz.submat(start_row,start_col,end_row,end_col) = block;
      
      if (b > 0){
        
        Toeplitz.submat(start_col,start_row,end_col,end_row) = block.t();
        
      }
      
      
    }
    
    
  }
  
  return(Toeplitz);
}



// [[Rcpp::export]]
arma::mat matrixform(
    const arma::vec& x
){
  int n = sqrt((double)x.n_elem);
  arma::mat out(n,n);
  int curelem=0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      
      out(i,j) = x(curelem);
      curelem++;
    }
  }
  
  return(out);
}
