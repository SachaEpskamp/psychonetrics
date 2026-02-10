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

    // --- Constant model slot data (Opt 3) ---
    // From model slots:
    std::string framework;        // model.slot("model")
    std::string estimator;        // model.slot("estimator")
    std::string distribution;     // model.slot("distribution")

    // Penalty data for PML (constant during optimization)
    double penalty_alpha;            // elastic net mixing (1=LASSO, 0=ridge)
    arma::vec penalty_lambda_vec;    // per-free-parameter penalty weights (length nFreePar)
    bool meanstructure;           // model.slot("meanstructure")
    Rcpp::List extramatrices;     // model.slot("extramatrices")
    Rcpp::List types;             // model.slot("types")

    // From sample slots:
    bool corinput;                // sample.slot("corinput")
    bool fullFIML;                // sample.slot("fullFIML")
    Rcpp::List sampleCovs;        // sample.slot("covs")
    Rcpp::List sampleMeans;       // sample.slot("means")
    Rcpp::List sampleThresholds;  // sample.slot("thresholds")
    Rcpp::List sampleSquares;     // sample.slot("squares") - Ising only
    Rcpp::List fimldata;          // sample.slot("fimldata")
    Rcpp::List WLS_W;             // sample.slot("WLS.W")

    // Pre-computed from groups/variables:
    int nGroup;                   // groups["id"].n_elem
    arma::vec nPerGroup;          // groups["nobs"]
    double nTotal;                // sum(nPerGroup)
    int nVar;                     // variables["id"].n_elem

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
