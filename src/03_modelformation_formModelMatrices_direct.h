#ifndef MODELMATS_DIRECT_H
#define MODELMATS_DIRECT_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Structure describing where a single parameter entry maps to in the matrix system
struct ParamEntry {
    int group;       // group index (0-based)
    int matrix;      // matrix index within group (0-based)
    int row;         // row in the matrix (0-based)
    int col;         // col in the matrix (0-based)
    bool symmetric;  // if true, also set (col, row)
};

// Pre-computed mapping for fast matrix formation
struct MatrixMapping {
    // Per-group, per-matrix: base matrix with fixed (par==0) values filled in
    // baseMatrices[group][matrix] = arma::mat
    std::vector<std::vector<arma::mat>> baseMatrices;

    // Matrix dimensions: dims[matrix] = {nrow, ncol}
    std::vector<std::pair<int,int>> matDims;

    // Matrix names in order
    std::vector<std::string> matNames;

    // Group names in order
    std::vector<std::string> groupNames;

    // Number of groups and matrices
    int nGroup;
    int nMat;

    // For each row in the parameter table (total params), store its mapping
    // totalParamMap[i] = ParamEntry for parameter table row i
    // Only entries with par > 0 OR par == 0 (fixed) are stored
    // We store ALL entries so we can rebuild from est vector
    std::vector<ParamEntry> allEntries;

    // For the free parameters: freeParamMap[free_par_index] = list of param table rows
    // that map to this free parameter. This lets us go from x[k] to all cells.
    std::vector<std::vector<int>> freeToTableRows;

    // Total number of parameters in the table
    int nParTotal;

    // Number of free parameters
    int nFreePar;

    // Whether each matrix is incomplete (fill with NA vs 0)
    std::vector<bool> matIncomplete;
};

// Build the mapping from an S4 model (call once)
MatrixMapping buildMatrixMapping(const S4& model);

// Use the mapping to form model matrices from a parameter vector x
// This is the fast replacement for formModelMatrices_cpp
Rcpp::List formModelMatrices_direct(
    const arma::vec& x,
    const MatrixMapping& mapping
);

// Wrapper that takes the S4 model (for backward compatibility / first call)
// Builds the mapping internally and then forms matrices
Rcpp::List formModelMatrices_cpp_direct(
    const arma::vec& x,
    const S4& model
);

#endif
