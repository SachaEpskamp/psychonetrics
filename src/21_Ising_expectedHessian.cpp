// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Expected Hamiltonian for the Spin distribution. The Hamiltonian includes the
// per-node quadratic (Blume-Capel) term delta_i x_i^2 (delta = 0 -> classical
// Ising):
//   H(x) = - sum_i tau_i x_i + sum_i delta_i x_i^2 - sum_{i<j} omega_ij x_i x_j
// [[Rcpp::export]]
double expHcpp(
    const arma::mat& states,
    const arma::vec& probabilities,
    const arma::mat& omega,
    const arma::vec& tau,
    const arma::vec& delta,
    const int nstate,
    const int nvar) {
  int s, i, j;

  // Value to store:
  double expH = 0;

  // Start looping:
  for (s=0; s < nstate; s++){
    // row:
    for (i=0; i < nvar; i++){
      expH -= probabilities[s] * tau[i] * states(s,i);
      expH += probabilities[s] * delta[i] * states(s,i) * states(s,i);

      //col:
      for (j=0; j<i; j++){
        expH -= probabilities[s] * omega(i,j) * states(s,i) * states(s,j);
      }
    }
  }

  // Return
  return expH;
}

// Same but for the sum of squares:
double expH2cpp(
    const arma::mat& states,
    const arma::vec& probabilities,
    const arma::mat& omega,
    const arma::vec& tau,
    const arma::vec& delta,
    const int nstate,
    const int nvar) {
  int s, i, j;

  // Value to store:
  double expH = 0;

  // Start looping:
  for (s=0; s < nstate; s++){
    // Dummy for current state:
    double stateH = 0;

    // row:
    for (i=0; i < nvar; i++){
      stateH -= tau[i] * states(s,i);
      stateH += delta[i] * states(s,i) * states(s,i);

      //col:
      for (j=0; j<i; j++){
        stateH -= omega(i,j) * states(s,i) * states(s,j);
      }
    }
    expH += probabilities[s] *  stateH * stateH;
  }

  // Return
  return expH;
}

// Asymptotic cov matrix
// [[Rcpp::export]]
arma::mat expHessianCpp(
    const arma::mat& states,
    const arma::vec& probabilities,
    const arma::mat& omega,
    const arma::vec& tau,
    const arma::vec& delta,
    double beta,
    const int nstate,
    const int nvar) {
  // Order of parameters is (thresholds, omega[lower.tri], delta, beta). The
  // delta block is inserted between omega and beta so that beta remains the
  // last parameter (the beta-* fills below index it as nElement-1).
  int i, j, g, h, s, ii, jj;



  // Compute expected hamiltonian:
  double expH = expHcpp(
    states,
    probabilities,
    omega,
    tau,
    delta,
    nstate,
    nvar);

  // And its square:
  double expH2 = expH2cpp(
    states,
    probabilities,
    omega,
    tau,
    delta,
    nstate,
    nvar);
  
  // Asyptotic 1st, 2nd, 3rd and 4th order cov matrix:
  std::vector<double> firstoder(nvar);
  std::vector<double> secondorder(nvar * nvar);
  std::vector<double> thirdorder(nvar * nvar * nvar);
  std::vector<double> fourthorder(nvar * nvar * nvar * nvar);
  
  // Trick from Joris 
  auto fir = [&](size_t s1)->double& { return firstoder[s1]; };
  auto sec = [&](size_t s1, size_t s2)->double& { return secondorder[s1 + s2 * nvar]; };
  auto thi = [&](size_t s1, size_t s2, size_t s3)->double& { return thirdorder[s1 + s2 * nvar + s3 * nvar * nvar]; };
  auto four = [&](size_t s1, size_t s2, size_t s3, size_t s4)->double& { return fourthorder[s1 + s2 * nvar+ s3 * nvar * nvar + s4 * nvar * nvar * nvar]; };
  
  // I also need the first and second order times hamiltonian:
  std::vector<double> firs_times_H(nvar);
  std::vector<double> second_times_H(nvar * nvar);
  // Trick from Joris 
  auto fir_H = [&](size_t s1)->double& { return firs_times_H[s1]; };
  auto sec_H = [&](size_t s1, size_t s2)->double& { return second_times_H[s1 + s2 * nvar]; };
  
  // Vector to store Hamiltonians:
  std::vector<double> allHs(nstate);
  
  // First compute all moments:
  // Loop over all states: FIXME: Currently filling all elements
  for (s=0; s < nstate; s++){
    // first compute the Hamiltonian:
    allHs[s] = 0;
    for (ii=0;ii<nvar;ii++){

      // Update threshold part:
      allHs[s] -= tau[ii] * states(s,ii);

      // Update quadratic (Blume-Capel) part:
      allHs[s] += delta[ii] * states(s,ii) * states(s,ii);

      for (jj=0;jj<ii;jj++){
        // Update network part:
        allHs[s] -= omega(ii,jj) * states(s,ii)  * states(s,jj);
      }
    }
    
    
    // row:
    for (i=0; i < nvar; i++){
      // First order moment:
      if (s==0){
        fir(i) = 0; // Initialize
        fir_H(i) = 0;
      }
      // Update first order moment:
      fir(i) += probabilities[s] * states(s,i);
      // Update Hamiltonian:
      fir_H(i) +=  probabilities[s] * states(s,i) * allHs[s];
      
      // Column (ignore diagonal):
      // for (j=0; j<=i; j++){
      for (j=0; j<nvar; j++){
        // Initialize second order moment:
        if (s==0){
          sec(i,j) = 0; // Initialize
          sec_H(i,j) = 0;
        }
        // Update second order moment:
        sec(i,j) += probabilities[s] * states(s,i) * states(s,j);
        // Update Hamiltonian:
        sec_H(i,j) +=  probabilities[s] * states(s,i) * states(s,j) * allHs[s];
        
        
        // Third slice:
        // for (g=0; g<=j; g++){
        for (g=0; g<nvar; g++){
          // Initialize second order moment:
          if (s==0){
            thi(i,j,g) = 0; // Initialize
          }
          // Update second order moment:
          thi(i,j,g) += probabilities[s] * states(s,i) * states(s,j) * states(s,g);
          
          // Fourth slice:
          // for (h=0; h<=g; h++){
          for (h=0; h<nvar; h++){
            // Initialize second order moment:
            if (s==0){
              four(i,j,g,h) = 0; // Initialize
            }
            // Update second order moment:
            four(i,j,g,h) += probabilities[s] * states(s,i) * states(s,j) * states(s,g) * states(s,h);
            
          }
        }
        
      }
    } 
    
    
    
    
  }
  
  // Now it is time to form the expected Hessian matrix. Parameter order is
  // (thresholds, omega[lower.tri], delta, beta): nvar thresholds, nvar*(nvar-1)/2
  // network parameters, nvar quadratic (delta) parameters, and 1 beta.
  int nE_omega = nvar * (nvar - 1) / 2;
  int deltaOffset = nvar + nE_omega;            // first index of the delta block
  int nElement = nvar + nE_omega + nvar + 1;    // beta is the last parameter
  arma::mat Hessian = arma::mat(nElement,nElement,fill::zeros);
  
  // Fill threshold-threshold part:
  for (i=0;i<nvar;i++){
    for (j=0;j<=i;j++){
      Hessian(i,j) = 2.0 * pow(beta, 2.0) * (
        sec(i,j) - fir(i) * fir(j)
      );
    }  
  }
  
  // Fill the graph-threshold part:
  
  // Loop over lower triangle of the graph:
  int currow = nvar;
  int curcol = 0;
  for (j=0;j<nvar;j++){
    for (i=j+1;i<nvar;i++){
      // Now loop over thresholds:
      for (g=0;g<nvar;g++){
        Hessian(currow,curcol) = 
          2.0 * pow(beta, 2.0) * (
          thi(i,j,g) - sec(i,j) * fir(g)
        );
        curcol++;
      }
      currow++;
      curcol = 0;
    }  
  }
  
  // Fill graph-graph part:
  // Loop over lower triangle of the graph:
  currow = nvar;
  curcol = nvar;
  
  for (j=0;j<nvar;j++){
    for (i=j+1;i<nvar;i++){
      
      // Loop over lower triangle of the graph again:
      for (jj=0;jj<nvar;jj++){
        for (ii=jj+1;ii<nvar;ii++){
          
          // FIXME: This is just pure lazyness....
          if (currow >= curcol){
            Hessian(currow,curcol) = 
              2.0 * pow(beta, 2.0) * (
              four(i,j,ii,jj) - sec(i,j) * sec(ii,jj)   
            );            
          }
          
          currow++;
        }
      }
      currow=nvar;
      curcol++;
    }  
  }
  
  // ---- delta (Blume-Capel quadratic) blocks ----
  // The expected-Hessian entry for parameters (a, b) is 2 * Cov(g_a, g_b) of the
  // per-observation score features: g_tau_i = beta x_i, g_omega_ab = beta x_a x_b,
  // g_delta_i = -beta x_i^2, g_beta = -H. The delta features carry a minus sign,
  // which flips the sign relative to the corresponding tau/omega blocks.

  // Fill the delta - threshold part:
  //   2 Cov(-beta x_i^2, beta x_j) = -2 beta^2 ( E[x_i^2 x_j] - E[x_i^2] E[x_j] )
  for (i=0;i<nvar;i++){
    for (j=0;j<nvar;j++){
      Hessian(deltaOffset + i, j) =
        -2.0 * pow(beta, 2.0) * ( thi(i,i,j) - sec(i,i) * fir(j) );
    }
  }

  // Fill the delta - omega part:
  //   2 Cov(-beta x_i^2, beta x_a x_b) = -2 beta^2 ( E[x_i^2 x_a x_b] - E[x_i^2] E[x_a x_b] )
  // omega edges (a,b) with a>b in the same column-major lower-triangle order as above.
  for (i=0;i<nvar;i++){
    int curcol_d = nvar;
    for (j=0;j<nvar;j++){
      for (ii=j+1;ii<nvar;ii++){
        Hessian(deltaOffset + i, curcol_d) =
          -2.0 * pow(beta, 2.0) * ( four(i,i,ii,j) - sec(i,i) * sec(ii,j) );
        curcol_d++;
      }
    }
  }

  // Fill the delta - delta part (lower triangle):
  //   2 Cov(-beta x_i^2, -beta x_j^2) = 2 beta^2 ( E[x_i^2 x_j^2] - E[x_i^2] E[x_j^2] )
  for (i=0;i<nvar;i++){
    for (j=0;j<=i;j++){
      Hessian(deltaOffset + i, deltaOffset + j) =
        2.0 * pow(beta, 2.0) * ( four(i,i,j,j) - sec(i,i) * sec(j,j) );
    }
  }

  // Fill the beta - delta part (beta is the last parameter, row nElement-1):
  //   2 Cov(-beta x_i^2, -H) = 2 beta ( E[x_i^2 H] - E[x_i^2] E[H] )
  for (i=0;i<nvar;i++){
    Hessian(nElement-1, deltaOffset + i) = 2.0 * beta * ( sec_H(i,i) - sec(i,i) * expH );
  }

  // Fill the beta - threshold part:
  for (i=0;i<nvar;i++){
    Hessian(nElement-1,i) =  2.0 * beta * (fir(i) * expH - fir_H(i));
  }

  // Fill the beta - graph part:
  curcol = nvar;
  for (j=0;j<nvar;j++){
    for (i=j+1;i<nvar;i++){
      Hessian(nElement-1,curcol) =  2.0 * beta * (sec(i,j) * expH - sec_H(i,j));
      curcol++;
    }
  }
  

  // Fil the beta - beta part:
  Hessian(nElement-1,nElement-1) =  2.0 * (expH2 - expH * expH);

  // Fill the upper triangle:
  for (j=0;j<nElement;j++){
    for (i=0;i<j;i++){
      Hessian(i,j) = Hessian(j,i);
    }
  }
  

  // Return
  return Hessian; // Dummy return omega
}
