#ifndef RCPPHELPERS_H
#define RCPPHELPERS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::vec eig_sym_cpp(
        arma::mat X
);

bool sympd_cpp(
                arma::mat X
);

Rcpp::List solve_symmetric_cpp(
        arma::mat X,
        bool logdet = false,
        double epsilon  = 1.490116e-08
);

arma::mat solve_symmetric_cpp_matrixonly(
        arma::mat X,
        double epsilon  = 1.490116e-08
);

arma::mat solve_symmetric_cpp_matrixonly_withcheck(
                arma::mat X,
                bool& proper
);


arma::mat bdiag_psychonetrics(
                Rcpp::List mats
);

arma::vec vech(
                arma::mat X,
                bool diag = true
);

arma::vec seq_len_inds(
                int start,
                int n
);

arma::mat cov2cor_cpp(
                arma::mat X
);

arma::mat wi2net_cpp(
                const arma::mat& X
);

arma::mat SDmat(
                const arma::mat& X
);

bool anyNon0(
                const arma::mat& X
);

arma::mat invSDmat(
                const arma::mat& X
);

arma::mat cbind_psychonetrics(
                const Rcpp::List  mats
);

void growlist(
                Rcpp::List& X,
                const Rcpp::List Y
);

arma::vec parVector_cpp(
                const S4& model
);

arma::mat computePDC_cpp(
                const arma::mat& beta,
                const arma::mat& kappa,
                const arma::mat& sigma
);

arma::mat blockToeplitz_cpp(
                const Rcpp::List& X
);

arma::mat matrixform(
                const arma::vec& x
);

#endif
