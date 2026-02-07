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

// Static cache: single workspace + its SEXP key
static std::unique_ptr<OptimWorkspace> s_cachedWorkspace;
static SEXP s_cachedModelSEXP = R_NilValue;

// Get or build workspace using static cache
const OptimWorkspace& getOrBuildWorkspace(const S4& model) {
    SEXP currentSEXP = (SEXP)model;

    // Cache hit: same model object, return cached workspace
    if (s_cachedWorkspace && s_cachedModelSEXP == currentSEXP) {
        return *s_cachedWorkspace;
    }

    // Cache miss: build new workspace and cache it
    s_cachedWorkspace = std::make_unique<OptimWorkspace>(buildOptimWorkspace(model));
    s_cachedModelSEXP = currentSEXP;

    return *s_cachedWorkspace;
}

// Manually invalidate the cache
// [[Rcpp::export]]
void invalidateWorkspaceCache() {
    s_cachedWorkspace.reset();
    s_cachedModelSEXP = R_NilValue;
    // Also clear the prepareModel_cpp result cache (Optimization 4):
    invalidatePrepCache();
}
