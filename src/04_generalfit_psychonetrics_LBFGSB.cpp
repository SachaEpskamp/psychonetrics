// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// LBFGS++ (LBFGSpp) based optimizer for psychonetrics
// Uses the L-BFGS-B algorithm from https://github.com/yixuan/LBFGSpp
// MIT license — see inst/include/LBFGSpp-LICENSE.md

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h>
#include <LBFGSB.h>

#include "04_generalFit_implied_and_prepare.h"
#include "04_generalfit_fitfunction_cpp.h"
#include "04_generalFit_gradient_cpp.h"
#include "b_modelexpansion_updateModel_cpp.h"
#include "02_algebrahelpers_modelMatrix_cpp.h"
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Functor for LBFGS++ : combined fn + grad in one operator()
// ---------------------------------------------------------------------------
class PsychonetricsFunctor {
public:
    PsychonetricsFunctor(const S4& model)
        : model_(model), fncount_(0), grcount_(0)
    {
        // Pre-warm the static workspace cache for this model:
        getOrBuildWorkspace(model);
    }

    // LBFGS++ calls  f(x_eigen, grad_eigen)  and expects:
    //   - return value = objective function at x
    //   - grad_eigen overwritten with gradient at x
    double operator()(const Eigen::VectorXd& x_eigen, Eigen::VectorXd& grad_eigen)
    {
        const int n = x_eigen.size();

        // --- Eigen → Armadillo (zero-copy for x, allocate for grad) ---
        // const_cast is safe here: arma::vec with copy_aux_mem=false will not
        // modify the data; psychonetrics_fitfunction_cpp takes const ref.
        arma::vec x_arma(const_cast<double*>(x_eigen.data()), n, false, true);
        arma::vec grad_arma(n);

        // --- Fit function ---
        double fx = psychonetrics_fitfunction_cpp(x_arma, model_);
        fncount_++;

        // NaN / Inf protection: return a large but finite value so the
        // line-search can back-track instead of crashing.
        if (!std::isfinite(fx)) {
            grad_eigen.setZero();
            return 1e20;
        }

        // --- Gradient (reuses the cached prepareModel from the fit call) ---
        psychonetrics_gradient_cpp_inner(x_arma, grad_arma, model_);
        grcount_++;

        // --- Armadillo grad → Eigen grad (memcpy, ~nPar*8 bytes) ---
        std::memcpy(grad_eigen.data(), grad_arma.memptr(), n * sizeof(double));

        // Allow R user interrupts (Ctrl-C / ESC)
        Rcpp::checkUserInterrupt();

        return fx;
    }

    int fncount() const { return fncount_; }
    int grcount() const { return grcount_; }

private:
    S4 model_;
    int fncount_;
    int grcount_;
};

// ---------------------------------------------------------------------------
// Rcpp-exported entry point
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
S4 psychonetrics_lbfgsb_optimizer(
    S4 model,
    const arma::vec& lower,
    const arma::vec& upper,
    bool bounded = false
) {
    int i;

    // ---- Extract start vector from parameter table (same as existing optimizer) ----
    Rcpp::List pars = model.slot("parameters");
    arma::vec est    = pars["est"];
    arma::vec parnum = pars["par"];
    int nPar      = arma::max(parnum);
    int totalPar  = parnum.n_elem;

    // If 0 free parameters, nothing to optimise
    if (nPar == 0) {
        model.slot("computed") = true;
        return model;
    }

    // Build start vector (free-parameter ordering)
    arma::vec x_arma(nPar);
    for (i = 0; i < totalPar; i++) {
        if (parnum(i) > 0) {
            x_arma(parnum(i) - 1) = est(i);
        }
    }

    // ---- Convert start vector & bounds to Eigen ----
    Eigen::VectorXd x(nPar);
    std::memcpy(x.data(), x_arma.memptr(), nPar * sizeof(double));

    Eigen::VectorXd lb(nPar), ub(nPar);
    if (bounded) {
        std::memcpy(lb.data(), lower.memptr(), nPar * sizeof(double));
        std::memcpy(ub.data(), upper.memptr(), nPar * sizeof(double));
    } else {
        lb.setConstant(-std::numeric_limits<double>::infinity());
        ub.setConstant( std::numeric_limits<double>::infinity());
    }

    // ---- Configure LBFGS++ parameters ----
    // Key insight: nlminb converges primarily on relative change in objective
    // (rel.tol ~ 1.49e-8) and x-change (x.tol = 1.5e-8), NOT on gradient norm.
    // LBFGS++ has two convergence criteria:
    //   1) Gradient: ||Pg||_inf < max(epsilon, epsilon_rel * ||x||)
    //   2) Delta:    |f_{k-d} - f_k| < delta * max(1, |f_k|, |f_{k-d}|)
    // For SEM/FIML, criterion 2 (function change) is more appropriate because
    // the gradient can remain non-negligible at the optimum for ill-conditioned
    // problems, while the objective stabilises early.
    LBFGSpp::LBFGSBParam<double> param;
    param.m              = 10;       // more history for better Hessian approx on larger models
    param.epsilon        = 1e-5;     // gradient tolerance (less aggressive than before)
    param.epsilon_rel    = 1e-5;     // relative gradient tolerance
    param.past           = 1;        // look back 1 iteration for function change test
    param.delta          = 1e-8;     // stop when relative f-change < 1e-8 (≈ nlminb rel.tol)
    param.max_iterations = 10000;
    param.max_linesearch = 40;       // generous line-search budget
    param.max_submin     = 10;
    param.ftol           = 1e-4;     // sufficient decrease (Armijo condition)
    param.wolfe          = 0.9;      // curvature condition

    // ---- Create solver and functor ----
    LBFGSpp::LBFGSBSolver<double> solver(param);
    PsychonetricsFunctor fun(model);

    // ---- Run optimisation ----
    double fx;
    int convergence = 0;
    std::string message = "converged";

    try {
        solver.minimize(fun, x, fx, lb, ub);
    } catch (std::exception& e) {
        convergence = 1;
        message = std::string("LBFGS++ error: ") + e.what();
    }

    // ---- Copy optimal x back to Armadillo ----
    arma::vec x_result(nPar);
    std::memcpy(x_result.memptr(), x.data(), nPar * sizeof(double));

    // ---- Build output list (same structure as existing optimizer) ----
    Rcpp::List optimout = Rcpp::List::create(
        Rcpp::Named("par")         = x_result,
        Rcpp::Named("convergence") = convergence,
        Rcpp::Named("message")     = message,
        Rcpp::Named("value")       = fx,
        Rcpp::Named("fncount")     = fun.fncount(),
        Rcpp::Named("grcount")     = fun.grcount(),
        Rcpp::Named("optimizer")   = "LBFGS++"
    );

    // ---- Update model with optimal parameters ----
    model = updateModel_cpp(x_result, model, false);
    model.slot("computed") = true;
    model.slot("optim")    = optimout;

    return model;
}
