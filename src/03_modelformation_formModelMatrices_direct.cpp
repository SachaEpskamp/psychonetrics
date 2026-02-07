// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <string>
#include "03_modelformation_formModelMatrices_direct.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Build the mapping from an S4 model (call once before optimization)
MatrixMapping buildMatrixMapping(const S4& model) {
    MatrixMapping mapping;

    // Extract from S4:
    Rcpp::List mats = model.slot("matrices");
    Rcpp::List pars = model.slot("parameters");
    S4 sample = model.slot("sample");
    Rcpp::List groups = sample.slot("groups");

    // Groups:
    arma::vec allGroups = groups["id"];
    Rcpp::StringVector groupnames = groups["label"];
    mapping.nGroup = allGroups.n_elem;
    mapping.groupNames.resize(mapping.nGroup);
    for (int g = 0; g < mapping.nGroup; g++) {
        mapping.groupNames[g] = (std::string)groupnames(g);
    }

    // Matrices:
    Rcpp::StringVector matname = mats["name"];
    Rcpp::NumericVector matnrow = mats["nrow"];
    Rcpp::NumericVector matncol = mats["ncol"];
    Rcpp::LogicalVector matsym = mats["symmetrical"];
    Rcpp::LogicalVector matincomplete = mats["incomplete"];

    mapping.nMat = matname.length();
    mapping.matNames.resize(mapping.nMat);
    mapping.matDims.resize(mapping.nMat);
    mapping.matIncomplete.resize(mapping.nMat);

    for (int m = 0; m < mapping.nMat; m++) {
        mapping.matNames[m] = (std::string)matname(m);
        mapping.matDims[m] = std::make_pair((int)matnrow(m), (int)matncol(m));
        mapping.matIncomplete[m] = (bool)matincomplete(m);
    }

    // Build name-to-index map for matrices:
    std::map<std::string, int> matNameToIndex;
    for (int m = 0; m < mapping.nMat; m++) {
        matNameToIndex[mapping.matNames[m]] = m;
    }

    // Parameter table:
    Rcpp::StringVector par_mat = pars["matrix"];
    Rcpp::IntegerVector par_row = pars["row"];
    Rcpp::IntegerVector par_col = pars["col"];
    Rcpp::IntegerVector group_id = pars["group_id"];
    Rcpp::NumericVector par_est = pars["est"];
    arma::vec par_num_vec = pars["par"];
    Rcpp::IntegerVector par_num(par_num_vec.n_elem);
    for (int i = 0; i < (int)par_num_vec.n_elem; i++) {
        par_num[i] = (int)par_num_vec(i);
    }

    mapping.nParTotal = par_mat.length();
    mapping.nFreePar = (int)max(par_num_vec);

    // Build allEntries:
    mapping.allEntries.resize(mapping.nParTotal);
    for (int p = 0; p < mapping.nParTotal; p++) {
        std::string mname = (std::string)par_mat(p);
        int midx = matNameToIndex[mname];
        mapping.allEntries[p].group = group_id(p) - 1;  // 0-based
        mapping.allEntries[p].matrix = midx;
        mapping.allEntries[p].row = par_row(p) - 1;      // 0-based
        mapping.allEntries[p].col = par_col(p) - 1;      // 0-based
        mapping.allEntries[p].symmetric = (bool)matsym(midx);
    }

    // Build freeToTableRows:
    mapping.freeToTableRows.resize(mapping.nFreePar);
    for (int p = 0; p < mapping.nParTotal; p++) {
        if (par_num[p] > 0) {
            mapping.freeToTableRows[par_num[p] - 1].push_back(p);
        }
    }

    // Build base matrices (with fixed parameter values pre-filled):
    mapping.baseMatrices.resize(mapping.nGroup);
    for (int g = 0; g < mapping.nGroup; g++) {
        mapping.baseMatrices[g].resize(mapping.nMat);
        for (int m = 0; m < mapping.nMat; m++) {
            int nr = mapping.matDims[m].first;
            int nc = mapping.matDims[m].second;
            if (mapping.matIncomplete[m]) {
                mapping.baseMatrices[g][m] = arma::mat(nr, nc, arma::fill::none);
                mapping.baseMatrices[g][m].fill(NA_REAL);
            } else {
                mapping.baseMatrices[g][m] = arma::mat(nr, nc, arma::fill::zeros);
            }
        }
    }

    // Fill base matrices with fixed parameter values (par == 0):
    for (int p = 0; p < mapping.nParTotal; p++) {
        if (par_num[p] == 0) {
            const ParamEntry& e = mapping.allEntries[p];
            mapping.baseMatrices[e.group][e.matrix](e.row, e.col) = par_est(p);
            if (e.symmetric) {
                mapping.baseMatrices[e.group][e.matrix](e.col, e.row) = par_est(p);
            }
        }
    }

    return mapping;
}


// Use the pre-computed mapping to form model matrices from parameter vector x
// This is O(nFreePar) instead of O(nMat * nParTotal) with string comparisons
Rcpp::List formModelMatrices_direct(
    const arma::vec& x,
    const MatrixMapping& mapping
) {
    // Start with copies of the base matrices (fixed values already in place):
    std::vector<std::vector<arma::mat>> mats(mapping.nGroup);
    for (int g = 0; g < mapping.nGroup; g++) {
        mats[g].resize(mapping.nMat);
        for (int m = 0; m < mapping.nMat; m++) {
            mats[g][m] = mapping.baseMatrices[g][m]; // copy
        }
    }

    // Direct update: for each free parameter, write its value to all mapped cells
    for (int k = 0; k < mapping.nFreePar; k++) {
        double val = x(k);
        for (int idx : mapping.freeToTableRows[k]) {
            const ParamEntry& e = mapping.allEntries[idx];
            mats[e.group][e.matrix](e.row, e.col) = val;
            if (e.symmetric) {
                mats[e.group][e.matrix](e.col, e.row) = val;
            }
        }
    }

    // Convert to Rcpp::List format matching formModelMatrices_cpp output:
    Rcpp::List result;
    for (int g = 0; g < mapping.nGroup; g++) {
        Rcpp::List grouplist;
        for (int m = 0; m < mapping.nMat; m++) {
            grouplist[mapping.matNames[m]] = mats[g][m];
        }
        result[mapping.groupNames[g]] = grouplist;
    }

    return result;
}


// Convenience wrapper: builds mapping from S4 model, then forms matrices
// This is exported to R for testing purposes
// [[Rcpp::export]]
Rcpp::List formModelMatrices_cpp_direct(
    const arma::vec& x,
    const S4& model
) {
    MatrixMapping mapping = buildMatrixMapping(model);
    return formModelMatrices_direct(x, mapping);
}

// Benchmark helper: build mapping, then call formModelMatrices_direct N times
// This measures the amortized cost of the direct method with mapping cached
// [[Rcpp::export]]
double benchmark_formModelMatrices_direct(
    const arma::vec& x,
    const S4& model,
    int n_iter = 1000
) {
    // Build mapping once:
    MatrixMapping mapping = buildMatrixMapping(model);

    // Time N iterations of the direct method only:
    for (int i = 0; i < n_iter; i++) {
        Rcpp::List result = formModelMatrices_direct(x, mapping);
    }

    return (double)n_iter; // just return count; timing done in R
}
