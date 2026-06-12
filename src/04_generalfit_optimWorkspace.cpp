// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <memory>
#include "04_generalfit_optimWorkspace.h"
#include "04_generalFit_implied_and_prepare.h"
#include "03_modelformation_formModelMatrices_direct.h"
#include "02_algebrahelpers_modelMatrix_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Build a workspace from an S4 model
OptimWorkspace buildOptimWorkspace(const S4& model) {
    OptimWorkspace ws;

    // Build the matrix mapping (the most expensive part to cache)
    ws.mapping = buildMatrixMapping(model);

    // Build the M matrix (parameter-to-free-parameter mapping)
    ws.Mmatrix = Mmatrix_cpp_list(model.slot("parameters"));

    // Extract parameter info
    Rcpp::List pars = model.slot("parameters");
    ws.parnum = Rcpp::as<arma::vec>(pars["par"]);
    ws.nParTotal = ws.parnum.n_elem;
    ws.nFreePar = (int)max(ws.parnum);

    // --- Constant model slot data (Opt 3) ---
    // From model slots:
    ws.framework = Rcpp::as<std::string>(model.slot("model"));
    ws.estimator = Rcpp::as<std::string>(model.slot("estimator"));
    ws.distribution = Rcpp::as<std::string>(model.slot("distribution"));

    // Build penalty data for PML/PFIML (replicates penaltyVector() in C++)
    if (ws.estimator == "PML" || ws.estimator == "PFIML") {
        Rcpp::List penalty = Rcpp::as<Rcpp::List>(model.slot("penalty"));
        ws.penalty_alpha = Rcpp::as<double>(penalty["alpha"]);
        arma::vec penalty_lambda_col = Rcpp::as<arma::vec>(pars["penalty_lambda"]);
        ws.penalty_lambda_vec.zeros(ws.nFreePar);
        for (int i = 0; i < ws.nParTotal; i++) {
            int p = (int)ws.parnum(i);
            if (p > 0 && ws.penalty_lambda_vec(p - 1) == 0.0) {
                double lam_val = penalty_lambda_col(i);
                // NA in R becomes NaN in C++; treat as 0 (auto-select not yet resolved)
                ws.penalty_lambda_vec(p - 1) = std::isnan(lam_val) ? 0.0 : lam_val;
            }
        }
    } else {
        ws.penalty_alpha = 1.0;
        ws.penalty_lambda_vec.zeros(ws.nFreePar);
    }
    ws.meanstructure = Rcpp::as<bool>(model.slot("meanstructure"));
    ws.extramatrices = Rcpp::as<Rcpp::List>(model.slot("extramatrices"));
    ws.types = Rcpp::as<Rcpp::List>(model.slot("types"));

    // From sample slots:
    S4 sample = model.slot("sample");
    ws.corinput = Rcpp::as<bool>(sample.slot("corinput"));
    ws.fullFIML = Rcpp::as<bool>(sample.slot("fullFIML"));
    ws.sampleCovs = Rcpp::as<Rcpp::List>(sample.slot("covs"));
    ws.sampleMeans = Rcpp::as<Rcpp::List>(sample.slot("means"));
    ws.sampleThresholds = Rcpp::as<Rcpp::List>(sample.slot("thresholds"));
    if (sample.hasSlot("squares")) {
        ws.sampleSquares = Rcpp::as<Rcpp::List>(sample.slot("squares"));
    }
    ws.fimldata = Rcpp::as<Rcpp::List>(sample.slot("fimldata"));
    // Two-level sufficient statistics (ml_lvm estimator = "ML"). Objects saved
    // before 0.15.31 lack the slot, so guard with hasSlot (twin of
    // get_twolevel_stats in R):
    if (sample.hasSlot("twolevel")) {
        ws.sampleTwolevel = Rcpp::as<Rcpp::List>(sample.slot("twolevel"));
    }
    ws.WLS_W = Rcpp::as<Rcpp::List>(sample.slot("WLS.W"));

    // Pre-computed from groups/variables:
    Rcpp::List groups = Rcpp::as<Rcpp::List>(sample.slot("groups"));
    arma::vec groupId = Rcpp::as<arma::vec>(groups["id"]);
    ws.nGroup = groupId.n_elem;
    ws.nPerGroup = Rcpp::as<arma::vec>(groups["nobs"]);
    ws.nTotal = arma::sum(ws.nPerGroup);
    Rcpp::List variables = Rcpp::as<Rcpp::List>(sample.slot("variables"));
    arma::vec varId = Rcpp::as<arma::vec>(variables["id"]);
    ws.nVar = varId.n_elem;

    // Store cache key
    ws.modelSEXP = (SEXP)model;

    return ws;
}

// Static cache: single workspace + its SEXP key.
//
// The key SEXP is protected with R_PreserveObject while it is the cache key.
// This makes the address-based keying sound:
//  - the garbage collector cannot free the keyed model, so no *different*
//    object can ever be allocated at the cached address while it is cached
//    (previously, after rm() + gc(), a new model at the same address would
//    silently receive the stale workspace);
//  - preservation also bumps the object's reference count, so any R-level
//    slot modification of the keyed model duplicates it (copy-on-write),
//    yielding a new address and therefore a cache miss + rebuild.
static std::unique_ptr<OptimWorkspace> s_cachedWorkspace;
static SEXP s_cachedModelSEXP = R_NilValue;

// Replace the cached key, releasing the previously preserved one:
static void setCachedModelKey(SEXP key) {
    if (s_cachedModelSEXP != R_NilValue) {
        R_ReleaseObject(s_cachedModelSEXP);
    }
    s_cachedModelSEXP = key;
    if (key != R_NilValue) {
        R_PreserveObject(key);
    }
}

// Get or build workspace using static cache
const OptimWorkspace& getOrBuildWorkspace(const S4& model) {
    SEXP currentSEXP = (SEXP)model;

    // Cache hit: same (live, preserved) model object, return cached workspace
    if (s_cachedWorkspace && s_cachedModelSEXP == currentSEXP) {
        return *s_cachedWorkspace;
    }

    // Cache miss: build new workspace and cache it
    s_cachedWorkspace = std::make_unique<OptimWorkspace>(buildOptimWorkspace(model));
    setCachedModelKey(currentSEXP);

    // The prepareModel_cpp result cache is keyed on the same model SEXP; it
    // is only sound while its key equals the preserved workspace key, so
    // clear it whenever the workspace key changes:
    invalidatePrepCache();

    return *s_cachedWorkspace;
}

// Manually invalidate the cache
// [[Rcpp::export]]
void invalidateWorkspaceCache() {
    s_cachedWorkspace.reset();
    setCachedModelKey(R_NilValue);
    // Also clear the prepareModel_cpp result cache (Optimization 4):
    invalidatePrepCache();
}

// Update only the penalty lambda vector in the cached workspace
// without rebuilding the entire workspace. Used during lambda grid search
// where only the penalty changes between iterations.
// [[Rcpp::export]]
void updateWorkspacePenaltyLambda(const arma::vec& new_lambda_vec, SEXP modelSEXP) {
    if (s_cachedWorkspace && (int)new_lambda_vec.n_elem == s_cachedWorkspace->nFreePar) {
        s_cachedWorkspace->penalty_lambda_vec = new_lambda_vec;
        s_cachedWorkspace->modelSEXP = modelSEXP;
        setCachedModelKey(modelSEXP);
    }
    // Invalidate prep cache since model parameters changed
    invalidatePrepCache();
}
