// =============================================================================
// Dissimilarity : Generalized Riemannian Multidimensional Scaling 
// =============================================================================

#include "RcppArmadillo.h"
#include "src_auxiliary.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// update a single vector by an index 'i'
arma::rowvec src_germ_update(int i, arma::rowvec Drow, arma::mat X, std::string geometry, 
                             int maxiter, double abstol, 
                             bool exact_search, 
                             double learning_rate, double weight_order){
  // parameters
  int n = X.n_rows;
  double n_div = static_cast<double>(n-1);
  int p = X.n_cols;
  int n_now = 0;
  arma::uword sub_idx = 0;
  arma::rowvec Dweight(n,fill::zeros);
  if (weight_order < 0.5){
    Dweight.fill(1.0);
  } else {
    for (int i=0; i<n; i++){
      Dweight(i) = std::pow(Drow(i), weight_order);
    }
  }

  // prepare : arrays
  arma::rowvec x_old = X.row(i);
  arma::rowvec x_new(p,fill::zeros);
  arma::rowvec x_grad(p,fill::zeros);
  arma::rowvec x_tmp(p,fill::zeros);
  
  arma::vec dist2mat(n,fill::zeros);
  double cutoff = 100.0*arma::datum::eps;
  
  arma::uvec nz_true = arma::find(Drow > cutoff);
  arma::uvec nz_current(2,fill::zeros);
  arma::uvec nz_nonzero(2,fill::zeros);
  
  arma::vec grid_alpha = arma::regspace<arma::vec>(0.1, 0.1, 2.0);
  int n_alpha = grid_alpha.n_elem;
  arma::vec grid_cost(n_alpha, fill::zeros);
  
  // prepare : function pointers
  XPtr<expPtr> fun_exp = decideExpPtr(geometry);
  XPtr<logPtr> fun_log = decideLogPtr(geometry);
  XPtr<distPtr> fun_dist = decideDistPtr(geometry);
  
  expPtr general_exp = *fun_exp;
  logPtr general_log = *fun_log;
  distPtr general_dist = *fun_dist;
  
// 
//   typedef arma::rowvec (*logPtr)(arma::rowvec x, arma::rowvec y);
//   typedef arma::rowvec (*expPtr)(arma::rowvec x, arma::rowvec u, double t);
//   
//   XPtr<distPtr> decideDistPtr(std::string geometry);   // pairwise distance
//   XPtr<innerPtr> decideInnerPtr(std::string geometry); // inner product on TM
//   XPtr<tnormPtr> decideTnormPtr(std::string geometry); // norm of tangent vector
//   XPtr<logPtr> decideLogPtr(std::string geometry);     // logarithmic map
//   XPtr<expPtr> decideExpPtr(std::string geometry);     // exponential map
//   
  
  // iterate
  double increment = 1000.0;
  int itcount = 0;
  while (increment > abstol){
    // compute all the distances
    dist2mat = aux_distVec2Mat(x_old, X, geometry);
    
    // non-zero distances
    nz_current.reset();
    nz_current = arma::find(dist2mat > cutoff);
    
    // current index to update
    nz_nonzero.reset();
    nz_nonzero = arma::intersect(nz_current, nz_true);
    if (nz_nonzero.is_empty()){
      break;
    }
    n_now = nz_nonzero.n_elem;
    
    // update a gradient
    x_grad.fill(0.0);
    for (int i=0; i<n_now; i++){
      // current index for the entire ones
      sub_idx = nz_nonzero(i);
      
      // directly update
      x_grad += (4.0/n_div)*((Drow(sub_idx)-dist2mat(sub_idx))/(Dweight(sub_idx)*dist2mat(sub_idx)))*general_log(x_old,X.row(sub_idx)); 
    }
    
    // search
    if (exact_search){
      // compute 20 objective values (marginal)
      for (int i=0; i<n_alpha; i++){
        // temporary element
        x_tmp = general_exp(x_old, -grid_alpha(i)*x_grad, 1.0);
        
        // cost function calculation
        double tmpcost = 0.0;
        for (int j=0; j<n; j++){
          if (Dweight(j) > 0){
            tmpcost += (1.0/Dweight(j))*std::pow((general_dist(X.row(j), x_tmp)-Drow(j)), 2.0);
          }
        }
        grid_cost(i) = tmpcost;
      }
      
      // find the minimal
      arma::uword id_min = grid_cost.index_min();
      
      // update
      x_tmp = -grid_alpha(id_min)*x_grad;
      x_new = general_exp(x_old, x_tmp, 1.0);
    } else {
      // update with a fixed learning rate (alpha)
      x_tmp = -learning_rate*x_grad; 
      x_new = general_exp(x_old, x_tmp, 1.0);
    }
    
    // updating rule
    increment = arma::norm(x_old - x_new, 2)/static_cast<double>(n); // average
    itcount += 1;
    x_old = x_new;
    
    if (itcount >= maxiter){
      break;
    }
  }

  // return
  return(x_old);
} 

// overall cost computation
double src_germ_cost(arma::mat D, arma::mat X, std::string geometry, double weight_order){
  // params
  int n = X.n_rows;
  
  // pointer
  XPtr<distPtr> fun_dist = decideDistPtr(geometry);
  distPtr general_dist = *fun_dist;
  
  // iterate
  double output = 0.0;
  if (weight_order < 0.5){
    for (int i=0; i<(n-1); i++){
      for (int j=(i+1); j<n; j++){
        if (D(i,j) > 10*arma::datum::eps){
          output += (std::pow(D(i,j)-general_dist(X.row(i), X.row(j)), 2.0));
        }
      }
    } 
  } else if (weight_order < 1.5){
    for (int i=0; i<(n-1); i++){
      for (int j=(i+1); j<n; j++){
        if (D(i,j) > 10*arma::datum::eps){
          output += (std::pow(D(i,j)-general_dist(X.row(i), X.row(j)), 2.0))/D(i,j);
        }
      }
    } 
  } else {
    for (int i=0; i<(n-1); i++){
      for (int j=(i+1); j<n; j++){
        if (D(i,j) > 10*arma::datum::eps){
          output += (std::pow(D(i,j)-general_dist(X.row(i), X.row(j)), 2.0))/(D(i,j)*D(i,j));
        }
      }
    } 
  }
  return(output);
}

// MAIN FUNCTION FOR GERM ------------------------------------------------------
// D             : pairwise distance matrix
// p             : dimensionality
// geometry      : manifold geometry
// iter_outer    : maximum iteration for outer loop
// iter_inner    : maximum iteration for cyclic optimization
// tol_outer     : stopping criterion for outer loop
// tol_inner     : stopping criterion for cyclic optimization
// exact_search  : TRUE for exact search, FALSE for fixed search
// learning_rate : increment for inexact line search
// weight_order  : 0 for metric MDS, 1 for sammon, or 2 for kamada-kawai
// [[Rcpp::export]]
arma::mat src_germ(arma::mat &D, int p, std::string geometry, 
                   int iter_outer, int iter_inner,
                   double tol_outer, double tol_inner, 
                   bool exact_search, double learning_rate,
                   double weight_order){
  // param
  int n = D.n_rows; // number of elements

  // initialization
  arma::mat old_X = aux_initembed(n,p,geometry,0.5); // sd=0.5
  arma::mat new_X = old_X;
  arma::rowvec tmpvec(old_X.n_cols, fill::zeros);
  
  double old_cost = src_germ_cost(D, old_X, geometry, weight_order);
  double new_cost = 100.0;
    
  // iteration
  double increment = 1000.0;
  int iter_count = 0;
  while (increment > tol_outer){
    // iterate over each item
    for (int i=0; i<n; i++){
      // compute
      tmpvec = src_germ_update(i, D.row(i), new_X, geometry, 
                               iter_inner, tol_inner, 
                               exact_search, learning_rate, weight_order);
      new_X.row(i) = tmpvec;
    }
    
    // updating rule
    increment = arma::norm(old_X-new_X,"fro");
    new_cost  = src_germ_cost(D, new_X, geometry, weight_order);
    if (new_cost >= old_cost){ // the update does not do any good.
      break;
    }
    iter_count += 1;
    old_X = new_X;
    if (iter_count >= iter_outer){
      break;
    }
    
    // print
    Rcpp::Rcout << "iteration " << iter_count << " complete : increment=" << increment << endl;
  }
  
  
  // return
  return(old_X);
}