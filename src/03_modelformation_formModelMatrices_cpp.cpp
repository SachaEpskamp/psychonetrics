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

// Eigenvalues of symmetric matrix:
// [[Rcpp::export]]
Rcpp::List formModelMatrices_cpp(
    const S4& model
){
  int i,j, m, g, p;
  
  // Obtain from S4:
  Rcpp::List mats = model.slot("matrices");
  Rcpp::List pars = model.slot("parameters");
  S4 sample = model.slot("sample");
  Rcpp::List groups = sample.slot("groups");
  
  // All groups:
  arma::vec allGroups = groups["id"];
  Rcpp::StringVector groupnames = groups["label"];
  
  // Things I need from the list of matrices:
  Rcpp::StringVector matname = mats["name"];
  Rcpp::NumericVector matnrow = mats["nrow"];
  Rcpp::NumericVector matncol = mats["ncol"];
  Rcpp::LogicalVector matsym = mats["symmetrical"];
  Rcpp::LogicalVector matlowtry = mats["lowertri"];
  Rcpp::LogicalVector matdiag = mats["diagonal"];
  Rcpp::LogicalVector matincomplete = mats["incomplete"];
  
  // Things I need from the parameter df:
  Rcpp::StringVector par_mat = pars["matrix"];
  Rcpp::IntegerVector par_row = pars["row"];
  Rcpp::IntegerVector par_col = pars["col"];
  Rcpp::IntegerVector group = pars["group_id"];
  Rcpp::NumericVector par_est = pars["est"];
  
  int nParTotal = par_mat.length();
  int nMat = matname.length();
  int nGroup = allGroups.n_elem;
  int nrow, ncol;
  std::string name;
  
  
  // Empty dummy matrices list:
  Rcpp::List matrices;

  // For each group:
  for (g=0; g < nGroup; g++){
    Rcpp::List grouplist;
    std::string groupname = (std::string)groupnames(g);

    
      // For every matrix:
    for (m=0; m<nMat; m++){
      nrow = matnrow(m);
      ncol = matncol(m);
      name = matname(m);
      
      
      // Empty matrix:
      arma::mat curmat(nrow,ncol);
      
      // If incomplete, fill with na. Otheriwse, fill with 0:
      for (i = 0; i < nrow; i++){
        for (j=0;j<ncol;j++){
          if (matincomplete(m)){
            curmat(i,j) = NA_REAL;
          } else {
            curmat(i,j) = 0;
          }
        }
      }
      
      // FIXME: loop over all rows in parameter table. this can be nicer...
      for (p=0;p<nParTotal;p++){
        if (par_mat[p] == name && group[p] == (g+1)){
          curmat(par_row[p]-1,par_col[p]-1) = par_est[p];
          if (matsym(m)){
            curmat(par_col[p]-1,par_row[p]-1) = par_est[p];
          }
        }
        
        
      }
  
      
      grouplist[name] = curmat;
    }
    
    
    matrices[groupname] = grouplist;
  }
  
  
  
 return(matrices); 
}
