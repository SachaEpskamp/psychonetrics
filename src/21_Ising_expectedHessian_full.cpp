// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Single-pass, pure-C++ expected Hessian for Ising group.
//
// Replaces the R callback expected_hessian_Ising_group, which itself called
// IsingSampler::IsingLikelihood (full state matrix materialized in R) and
// then expHessianCpp (which made three more passes through all 2^N states:
// expHcpp + expH2cpp + the main moment-accumulation loop).
//
// This implementation walks all 2^N states EXACTLY ONCE, accumulating Z,
// expH, expH2, and the first/second/third/fourth-order moments needed by
// the Hessian assembly. The fourth-order tensor is filled only on its
// canonical lower-tri-of-pairs sub-shape (the assembly itself only reads
// those entries; the upper triangle of the Hessian is mirrored at the end).

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Per-group expected Hessian. Drop-in replacement for the R function
// expected_hessian_Ising_group + expHessianCpp combo.
//
// Parameter ordering of the Hessian (matching expHessianCpp):
//   [tau_0 .. tau_{N-1},
//    omega(1,0), omega(2,0), .., omega(N-1,0),
//    omega(2,1), omega(3,1), .., omega(N-1,1),
//    .., omega(N-1, N-2),
//    beta]
// i.e. column-major lower triangle of omega between thresholds and beta.
arma::mat expected_hessian_Ising_group_full_cpp(
    const arma::mat& omega,
    const arma::vec& tau,
    double beta,
    const arma::vec& responses,
    double min_sum
){
    const int nvar = omega.n_rows;
    const int E    = nvar * (nvar - 1) / 2;          // number of edges (lower-tri of omega)
    const int nElement = nvar + E + 1;

    const int nResp = responses.n_elem;
    const bool alwaysUpdate = (min_sum == R_NegInf);

    // ---- Edge index helper: (i,j) with i>j -> column-major lower-tri ordering
    //   edge_idx(j=0, i=1)=0, (j=0,i=2)=1, .., (j=0,i=N-1)=N-2,
    //   (j=1, i=2)=N-1, ..
    // matches expHessianCpp's currow traversal order in the graph block.
    std::vector<int> edge_i(E), edge_j(E);
    {
        int k = 0;
        for (int j = 0; j < nvar; j++) {
            for (int i = j + 1; i < nvar; i++) {
                edge_i[k] = i;
                edge_j[k] = j;
                k++;
            }
        }
    }
    // Quick lookup for edge index from (i,j) pair (only for i>j):
    // not needed at runtime — moments are indexed by edge_idx directly.

    // ---- Accumulators (un-normalized: weight = pot, divide by Z at end) ----
    double Z = 0.0;
    double expH_u  = 0.0;   // sum pot * H
    double expH2_u = 0.0;   // sum pot * H^2

    arma::vec fir_u(nvar, fill::zeros);
    arma::vec fir_H_u(nvar, fill::zeros);

    // sec stored as full N x N (we'll only fill lower tri including diagonal)
    arma::mat sec_u(nvar, nvar, fill::zeros);
    arma::mat sec_H_u(nvar, nvar, fill::zeros);

    // thi: edge x N (third-order moments needed for graph-threshold block)
    arma::mat thi_u(E, nvar, fill::zeros);

    // four: edge-pair lower triangle (size E*(E+1)/2)
    // index by  fp_idx(a,b) = a*(a+1)/2 + b   for a >= b
    std::vector<double> four_u((size_t)E * (E + 1) / 2, 0.0);

    // ---- Enumerate all nResp^N states via a mixed-radix counter ----
    // (digit 0 fastest; reduces exactly to the binary r1/r2 enumeration when
    // nResp == 2, so the two-outcome expected Hessian is unchanged).
    std::vector<int> idx(nvar, 0);
    arma::vec curstate(nvar, fill::value(responses(0)));
    bool done = false;

    do {
        bool update = alwaysUpdate || (sum(curstate) >= min_sum);

        if (update) {
            // Hamiltonian for current state
            double H = 0.0;
            for (int i = 0; i < nvar; i++) {
                H -= tau(i) * curstate(i);
                for (int j = 0; j < i; j++) {
                    H -= omega(i, j) * curstate(i) * curstate(j);
                }
            }

            const double pot = std::exp(-1.0 * beta * H);
            const double potH  = pot * H;
            const double potH2 = pot * H * H;

            Z       += pot;
            expH_u  += potH;
            expH2_u += potH2;

            // First / second order moments and their H-weighted variants
            for (int i = 0; i < nvar; i++) {
                const double xi = curstate(i);
                const double pot_xi  = pot  * xi;
                const double potH_xi = potH * xi;
                fir_u(i)   += pot_xi;
                fir_H_u(i) += potH_xi;

                for (int j = 0; j <= i; j++) {
                    const double xj = curstate(j);
                    sec_u(i, j)   += pot_xi  * xj;
                    sec_H_u(i, j) += potH_xi * xj;
                }
            }

            // Third / fourth order moments — restricted to the canonical
            // shapes the Hessian assembly reads.
            for (int a = 0; a < E; a++) {
                const int i = edge_i[a];
                const int j = edge_j[a];
                const double xij = curstate(i) * curstate(j);
                const double pot_xij = pot * xij;

                // thi(i,j,g) for g = 0..nvar-1
                for (int g = 0; g < nvar; g++) {
                    thi_u(a, g) += pot_xij * curstate(g);
                }

                // four(edge_a, edge_b) only for b <= a (lower tri of edge x edge)
                size_t base = (size_t)a * (a + 1) / 2;
                for (int b = 0; b <= a; b++) {
                    const int ii = edge_i[b];
                    const int jj = edge_j[b];
                    four_u[base + b] += pot_xij * curstate(ii) * curstate(jj);
                }
            }
        }

        // Increment the mixed-radix counter (digit 0 fastest):
        int c = 0;
        while (c < nvar && idx[c] == nResp - 1) {
            idx[c] = 0;
            curstate(c) = responses(0);
            c++;
        }
        if (c == nvar) {
            done = true;
        } else {
            idx[c]++;
            curstate(c) = responses(idx[c]);
        }
    } while (!done);

    // ---- Normalize ----
    const double invZ = 1.0 / Z;
    const double expH  = expH_u  * invZ;
    const double expH2 = expH2_u * invZ;

    arma::vec fir   = fir_u   * invZ;
    arma::vec fir_H = fir_H_u * invZ;

    // sec/sec_H: only lower tri populated, but assembly reads sec(i,j) for j<=i
    // and beta-graph reads sec_H(i,j) for i>j (also lower tri). No need to mirror.
    arma::mat sec   = sec_u   * invZ;
    arma::mat sec_H = sec_H_u * invZ;

    arma::mat thi   = thi_u   * invZ;
    // four_u in-place divide:
    {
        const size_t fsz = four_u.size();
        for (size_t k = 0; k < fsz; k++) four_u[k] *= invZ;
    }
    auto four_idx = [&](int a, int b) -> double {
        // assembly only ever reads with a >= b; if a < b swap (won't happen
        // in normal assembly, but keep for safety):
        if (a < b) std::swap(a, b);
        return four_u[(size_t)a * (a + 1) / 2 + b];
    };

    // ---- Assemble Hessian ----
    const double beta2 = beta * beta;
    arma::mat Hessian(nElement, nElement, fill::zeros);

    // tau-tau block (lower tri)
    for (int i = 0; i < nvar; i++) {
        for (int j = 0; j <= i; j++) {
            Hessian(i, j) = 2.0 * beta2 * (sec(i, j) - fir(i) * fir(j));
        }
    }

    // omega-tau block: rows = nvar..nvar+E-1, cols = 0..nvar-1
    for (int a = 0; a < E; a++) {
        const int i = edge_i[a];
        const int j = edge_j[a];
        const int row = nvar + a;
        for (int g = 0; g < nvar; g++) {
            Hessian(row, g) = 2.0 * beta2 * (thi(a, g) - sec(i, j) * fir(g));
        }
    }

    // omega-omega block: rows/cols = nvar..nvar+E-1, only fill row >= col
    for (int b = 0; b < E; b++) {
        const int ii = edge_i[b];
        const int jj = edge_j[b];
        const int col = nvar + b;
        for (int a = b; a < E; a++) {
            const int i = edge_i[a];
            const int j = edge_j[a];
            const int row = nvar + a;
            Hessian(row, col) = 2.0 * beta2 *
                (four_idx(a, b) - sec(i, j) * sec(ii, jj));
        }
    }

    // beta-tau block (last row)
    for (int i = 0; i < nvar; i++) {
        Hessian(nElement - 1, i) = 2.0 * beta * (fir(i) * expH - fir_H(i));
    }

    // beta-omega block (last row)
    for (int a = 0; a < E; a++) {
        const int i = edge_i[a];
        const int j = edge_j[a];
        Hessian(nElement - 1, nvar + a) =
            2.0 * beta * (sec(i, j) * expH - sec_H(i, j));
    }

    // beta-beta diagonal
    Hessian(nElement - 1, nElement - 1) = 2.0 * (expH2 - expH * expH);

    // Mirror lower triangle to upper
    for (int j = 0; j < nElement; j++) {
        for (int i = 0; i < j; i++) {
            Hessian(i, j) = Hessian(j, i);
        }
    }

    return Hessian;
}


// Top-level entry: takes prep (the prepared Ising model) and returns the
// block-diagonal expected Hessian over groups, weighted by group share.
// [[Rcpp::export]]
arma::mat expected_hessian_Ising_full_cpp(const Rcpp::List& prep) {
    Rcpp::List groupModels = prep["groupModels"];
    arma::vec  nPerGroup   = prep["nPerGroup"];
    double     nTotal      = prep["nTotal"];
    int nGroup = groupModels.length();

    Rcpp::List perGroup(nGroup);

    for (int g = 0; g < nGroup; g++) {
        Rcpp::List grp = groupModels[g];
        arma::mat omega    = grp["omega"];
        arma::vec tau      = grp["tau"];
        double    beta     = grp["beta"];
        arma::vec responses = grp["responses"];

        // min_sum lives in extramatrices (already merged into the grouplist
        // by prepare_Ising_cpp via growlist):
        double min_sum;
        if (grp.containsElementNamed("min_sum")) {
            min_sum = grp["min_sum"];
        } else {
            min_sum = R_NegInf;
        }

        arma::mat H = expected_hessian_Ising_group_full_cpp(
            omega, tau, beta, responses, min_sum);

        // Weight by group size
        H *= (nPerGroup(g) / nTotal);
        perGroup[g] = H;
    }

    return bdiag_psychonetrics(perGroup);
}
