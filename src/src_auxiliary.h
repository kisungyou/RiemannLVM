#ifndef SRC_AUXILIARY_H_
#define SRC_AUXILIARY_H_

#include "RcppArmadillo.h"
#include <string>

using namespace Rcpp;
using namespace arma;
using namespace std;

/* MISCELLANEOUS ROUTINES 
 * aux_equalstr       : boolean for equal strings
 * lorentzian_product : Lorentzian product on hyperboloid
 * 
 * GEOMETRIC ROUTINES - 'euclidean', 'spherical', and 'hyperbolic'
 * 
 * POINTERS
 * XPtr<distPtr> : distance
 * 
 * SHARED RESOURCES 
 * aux_pdist : compute pairwise distance matrix
 * 
 */

// MISCELLANEOUS ROUTINES  -----------------------------------------------------
bool aux_equalstr(std::string string1, std::string string2);
double lorentzian(arma::rowvec x, arma::rowvec y);




// GEOMETRIC_ROUTINES : EUCLIDEAN ----------------------------------------------
double euclidean_dist(arma::rowvec x, arma::rowvec y);
arma::rowvec euclidean_log(arma::rowvec x, arma::rowvec y);
arma::rowvec euclidean_exp(arma::rowvec x, arma::rowvec u, double t);
double euclidean_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v);
double euclidean_tnorm(arma::rowvec x, arma::rowvec u);
arma::rowvec euclidean_normalize(arma::rowvec x);




// GEOMETRIC_ROUTINES : SPHERICAL ----------------------------------------------
arma::rowvec spherical_proj(arma::rowvec x, arma::rowvec u);
arma::rowvec spherical_normalize(arma::rowvec x);

double spherical_dist(arma::rowvec x, arma::rowvec y);
arma::rowvec spherical_log(arma::rowvec x, arma::rowvec y);
arma::rowvec spherical_exp(arma::rowvec x, arma::rowvec u, double t);
double spherical_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v);
double spherical_tnorm(arma::rowvec x, arma::rowvec u);


// GEOMETRIC_ROUTINES : HYPERBOLIC(HYPERBOLOID) --------------------------------
// Note that we use the convention of the last element to be negated.
arma::rowvec hyperbolic_proj(arma::rowvec x, arma::rowvec u);
arma::rowvec hyperbolic_normalize(arma::rowvec x);

double hyperbolic_dist(arma::rowvec x, arma::rowvec y);
arma::rowvec hyperbolic_log(arma::rowvec x, arma::rowvec y);
arma::rowvec hyperbolic_exp(arma::rowvec x, arma::rowvec u, double t);
double hyperbolic_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v);
double hyperbolic_tnorm(arma::rowvec x, arma::rowvec u);


// POINTERS --------------------------------------------------------------------
typedef double (*distPtr)(arma::rowvec x, arma::rowvec y);
typedef double (*innerPtr)(arma::rowvec x, arma::rowvec u, arma::rowvec v);
typedef double (*tnormPtr)(arma::rowvec x, arma::rowvec u);
typedef arma::rowvec (*logPtr)(arma::rowvec x, arma::rowvec y);
typedef arma::rowvec (*expPtr)(arma::rowvec x, arma::rowvec u, double t);

XPtr<distPtr> decideDistPtr(std::string geometry);   // pairwise distance
XPtr<innerPtr> decideInnerPtr(std::string geometry); // inner product on TM
XPtr<tnormPtr> decideTnormPtr(std::string geometry); // norm of tangent vector
XPtr<logPtr> decideLogPtr(std::string geometry);     // logarithmic map
XPtr<expPtr> decideExpPtr(std::string geometry);     // exponential map

  

// SHARED RESOURCES ------------------------------------------------------------
arma::mat aux_pdist(arma::mat &X, std::string geometry);
arma::mat aux_pdist2(arma::mat &X, arma::mat &Y, std::string geometry);
arma::vec aux_distVec2Mat(arma::rowvec x, arma::mat Y, std::string geometry);
arma::mat aux_initembed(int n, int p, std::string geometry, double sd); // wrapped normal approach

#endif