#ifndef OPTIMWORKSPACE_H
#define OPTIMWORKSPACE_H

#include <RcppArmadillo.h>
#include <memory>
#include "03_modelformation_formModelMatrices_direct.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Workspace struct that caches constant data extracted from the S4 model.
// Built once before optimization, reused across all fn/grad calls.
struct OptimWorkspace {
    // Pre-built matrix mapping (eliminates buildMatrixMapping per call)
    MatrixMapping mapping;

    // Pre-built M matrix for gradient (eliminates Mmatrix_cpp_list per call)
    arma::sp_mat Mmatrix;

    // Pre-extracted parameter info
    arma::vec parnum;       // parameters["par"] vector
    int nParTotal;          // parnum.n_elem
    int nFreePar;           // max(parnum)

    // Cache key: raw SEXP pointer of the model object
    SEXP modelSEXP;
};

// Build a workspace from an S4 model (call once before optimization)
OptimWorkspace buildOptimWorkspace(const S4& model);

// Get or build workspace using static cache keyed by SEXP pointer.
// If the same model is passed again, returns the cached workspace.
// If a different model is passed, rebuilds and caches the new one.
const OptimWorkspace& getOrBuildWorkspace(const S4& model);

// Manually invalidate the cache (for testing or cleanup)
void invalidateWorkspaceCache();

#endif
