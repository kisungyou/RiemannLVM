// =============================================================================
// Riemannian Gaussian Latent Position Model
// (1) rglpm_loglkd : compute the log-likelihood given pars + distance matrix
// =============================================================================

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// (1) rglpm_loglkd 
// [[Rcpp::export]]
double rglpm_loglkd(arma::mat &A, const double tau, const double sig2, arma::mat &D){
  // parameters and auxiliaries
  int n = A.n_rows;
  double scaler = 2.0*sig2;
  
  // compute the loglkd directly
  double output = 0.0;
  double Pij    = 0.0; // probability
  double Lnow   = 0.0; // current entry
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i!=j){
        // compute an entry for Probability matrix P
        Pij  = tau*std::exp(-(D(i,j)*D(i,j))/scaler);
        Lnow = A(i,j)*std::log(Pij) + (1-A(i,j))*std::log(1.0-Pij);
        output += Lnow;
      }
    }
  }
  
  // return
  return(output);
}