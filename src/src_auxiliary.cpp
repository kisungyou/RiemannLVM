// =============================================================================
// All the routines
// =============================================================================


#include "src_auxiliary.h"

// MISCELLANEOUS ROUTINES  -----------------------------------------------------
// [[Rcpp::export]]
bool aux_equalstr(std::string string1, std::string string2){
  return string1.compare(string2) == 0;
}
double lorentzian(arma::rowvec x, arma::rowvec y){
  int n = x.n_elem;
  double output = 0.0;
  for (int i=0; i<(n-1); i++){
    output += x(i)*y(i);
  }
  output -= x(n-1)*y(n-1);
  return(output);
}


// GEOMETRIC_ROUTINES : EUCLIDEAN ----------------------------------------------
double euclidean_dist(arma::rowvec x, arma::rowvec y){
  return(arma::norm(x-y,2));
}
arma::rowvec euclidean_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec logvec = y-x;
  return(logvec);
}
arma::rowvec euclidean_exp(arma::rowvec x, arma::rowvec u, double t){
  arma::rowvec output = x + t*u;
  return(output);
}
double euclidean_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v){
  double output = arma::cdot(u,v);
  return(output);
}
double euclidean_tnorm(arma::rowvec x, arma::rowvec u){
  double gxuu = euclidean_inner(x,u,u);
  double output = std::sqrt(gxuu);
  return(output);
}
arma::rowvec euclidean_normalize(arma::rowvec x){
  return(x);
}

// GEOMETRIC_ROUTINES : SPHERICAL ----------------------------------------------
arma::rowvec spherical_proj(arma::rowvec x, arma::rowvec u){
  return(u-x*(arma::dot(x,u)));
}
arma::rowvec spherical_normalize(arma::rowvec x){
  double x_norm = arma::norm(x, 2);
  arma::rowvec output = x/x_norm;
  return(output);
}
double spherical_dist(arma::rowvec x, arma::rowvec y){
  arma::rowvec vecxy = x-y;
  double dotxy = arma::dot(x,y);
  
  if (arma::norm(vecxy, 2) < arma::datum::eps){
    return(0.0);
  } else if (std::sqrt(dotxy*dotxy) >= (1.0-arma::datum::eps)){
    return(arma::datum::pi);
  } else {
    return(std::acos(arma::dot(x, y)));  
  }
}
// [[Rcpp::export]]
arma::rowvec spherical_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec v = spherical_proj(x,y-x);
  double di = spherical_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, "fro");
    v = v*(di/nv);
  }
  return(v);
}
// [[Rcpp::export]]
arma::rowvec spherical_exp(arma::rowvec x, arma::rowvec u, double t){
  double nrm_td = arma::norm(t*u, 2); // theta
  arma::rowvec out(x.n_elem, fill::zeros);
  if (nrm_td < 1e-15){ // very close
    out = x;
  } else {
    out = cos(nrm_td)*x + ((sin(nrm_td))/nrm_td)*t*u;
    out /= arma::norm(out, 2);
  }
  return(out);
}
double spherical_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v){
  double output = arma::dot(u,v);
  return(output);
}
double spherical_tnorm(arma::rowvec x, arma::rowvec u){
  double gxuu   = spherical_inner(x,u,u);
  double output = std::sqrt(gxuu);
  return(output);
}




// GEOMETRIC_ROUTINES : HYPERBOLIC/HYPERBOLOID ---------------------------------
arma::rowvec hyperbolic_proj(arma::rowvec x, arma::rowvec u){
  arma::rowvec output = u + x*lorentzian(x,u);
  return(output);
}
arma::rowvec hyperbolic_normalize(arma::rowvec x){
  arma::rowvec output = x/std::sqrt(std::abs(lorentzian(x,x)));
  return(output);
}
double hyperbolic_dist(arma::rowvec x, arma::rowvec y){
  if (euclidean_dist(x,y) < 10*arma::datum::eps){
    return(0.0);
  } else {
    return(std::acosh(-lorentzian(x,y)));
  }
}
// [[Rcpp::export]]
arma::rowvec hyperbolic_log(arma::rowvec x, arma::rowvec y){
  // method 1 : direct formula
  // double alpha = -lorentzian(x,y);
  // arma::rowvec output = (std::acosh(alpha)/sqrt((alpha*alpha)-1))*(y-alpha*x);
  // return(output);
  
  // method 2 : numerically safer option by manopt
  double d = hyperbolic_dist(x,y);
  double a = d/std::sinh(d);
  if (std::isnan(a)){
    a = 1.0;
  }
  arma::rowvec u = hyperbolic_proj(x, y*a);
  return(u);
}
// [[Rcpp::export]]
arma::rowvec hyperbolic_exp(arma::rowvec x, arma::rowvec u, double t){
  arma::rowvec tu = t*u;
  double val_inner = lorentzian(tu,tu);
  double val_norm = 0;
  if (val_inner > 0){
    val_norm = std::sqrt(val_inner);
  }
  
  double a = std::cosh(val_norm);
  double b = std::sinh(val_norm)/val_norm;
  if (std::isnan(b)){
    b = 1.0;
  }
  arma::rowvec y = x*a + tu*b;
  arma::rowvec output = hyperbolic_normalize(y);
  return(output);
}
double hyperbolic_inner(arma::rowvec x, arma::rowvec u, arma::rowvec v){
  double output = lorentzian(u,v);
  return(output);
}
double hyperbolic_tnorm(arma::rowvec x, arma::rowvec u){
  double gxuu   = hyperbolic_inner(x,u,u);
  double output = 0.0;
  if (gxuu > 0){
    output = std::sqrt(gxuu);
  }
  return(output); 
}

// POINTERS --------------------------------------------------------------------
typedef double (*distPtr)(arma::rowvec x, arma::rowvec y);
typedef double (*innerPtr)(arma::rowvec x, arma::rowvec u, arma::rowvec v);
typedef double (*tnormPtr)(arma::rowvec x, arma::rowvec u);
typedef arma::rowvec (*logPtr)(arma::rowvec x, arma::rowvec y);
typedef arma::rowvec (*expPtr)(arma::rowvec x, arma::rowvec u, double t);

XPtr<distPtr> decideDistPtr(std::string geometry){
  if (aux_equalstr(geometry, "euclidean")){
    return(XPtr<distPtr>(new distPtr(&euclidean_dist)));
  } else if (aux_equalstr(geometry, "spherical")){
    return(XPtr<distPtr>(new distPtr(&spherical_dist)));
  } else if (aux_equalstr(geometry, "hyperbolic")){
    return(XPtr<distPtr>(new distPtr(&hyperbolic_dist)));
  }else {
    return XPtr<distPtr>(R_NilValue);
  }
}
XPtr<innerPtr> decideInnerPtr(std::string geometry){
  if (aux_equalstr(geometry, "euclidean")){
    return(XPtr<innerPtr>(new innerPtr(&euclidean_inner)));
  } else if (aux_equalstr(geometry, "spherical")){
    return(XPtr<innerPtr>(new innerPtr(&spherical_inner)));
  } else if (aux_equalstr(geometry, "hyperbolic")){
    return(XPtr<innerPtr>(new innerPtr(&hyperbolic_inner)));
  }else {
    return XPtr<innerPtr>(R_NilValue);
  }
}
XPtr<tnormPtr> decideTnormPtr(std::string geometry){
  if (aux_equalstr(geometry, "euclidean")){
    return(XPtr<tnormPtr>(new tnormPtr(&euclidean_tnorm)));
  } else if (aux_equalstr(geometry, "spherical")){
    return(XPtr<tnormPtr>(new tnormPtr(&spherical_tnorm)));
  } else if (aux_equalstr(geometry, "hyperbolic")){
    return(XPtr<tnormPtr>(new tnormPtr(&hyperbolic_tnorm)));
  }else {
    return XPtr<tnormPtr>(R_NilValue);
  }
}
XPtr<logPtr> decideLogPtr(std::string geometry){
  if (aux_equalstr(geometry, "euclidean")){
    return(XPtr<logPtr>(new logPtr(&euclidean_log)));
  } else if (aux_equalstr(geometry, "spherical")){
    return(XPtr<logPtr>(new logPtr(&spherical_log)));
  } else if (aux_equalstr(geometry, "hyperbolic")){
    return(XPtr<logPtr>(new logPtr(&hyperbolic_log)));
  }else {
    return XPtr<logPtr>(R_NilValue);
  }
}
XPtr<expPtr> decideExpPtr(std::string geometry){
  if (aux_equalstr(geometry, "euclidean")){
    return(XPtr<expPtr>(new expPtr(&euclidean_exp)));
  } else if (aux_equalstr(geometry, "spherical")){
    return(XPtr<expPtr>(new expPtr(&spherical_exp)));
  } else if (aux_equalstr(geometry, "hyperbolic")){
    return(XPtr<expPtr>(new expPtr(&hyperbolic_exp)));
  }else {
    return XPtr<expPtr>(R_NilValue);
  }
}

// SHARED RESOURCES ------------------------------------------------------------
// [[Rcpp::export]]
arma::mat aux_pdist(arma::mat &X, std::string geometry){
  // parameters
  int N = X.n_rows;
  
  // get a pointer
  XPtr<distPtr> xpfun = decideDistPtr(geometry);
  distPtr pairwise_distance = *xpfun;
  
  // compute
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = pairwise_distance(X.row(i), X.row(j));
      output(j,i) = output(i,j);
    }
  }
  
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat aux_pdist2(arma::mat &X, arma::mat &Y, std::string geometry){
  // parameters
  int M = X.n_rows;
  int N = Y.n_rows;
  
  // get a pointer
  XPtr<distPtr> xpfun = decideDistPtr(geometry);
  distPtr pairwise_distance = *xpfun;
  
  // compute
  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      output(m,n) = pairwise_distance(X.row(m), X.row(n));
    }
  }
  // return
  return(output);
}
// [[Rcpp::export]]
arma::vec aux_distVec2Mat(arma::rowvec x, arma::mat Y, std::string geometry){
  // get a pointer
  XPtr<distPtr> xpfun = decideDistPtr(geometry);
  distPtr pairwise_distance = *xpfun;
  
  // parameter
  int N = Y.n_rows;
  
  // compute
  arma::vec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = pairwise_distance(x, Y.row(n));
  }
  
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat aux_initembed(int n, int p, std::string geometry, double sd){
  arma::mat output;
  if (aux_equalstr(geometry,"euclidean")){
    output = arma::randn(n,p,arma::distr_param(0.0, sd));
  } else if (aux_equalstr(geometry,"spherical")){
    // centered at the north pole
    arma::rowvec origin(p+1,fill::zeros);
    origin(p) = 1.0;
    
    // tangent vectors
    arma::mat tangents = arma::randn(n,p+1,arma::distr_param(0.0, sd));
    tangents.col(p).fill(0.0);
    
    // embed
    output = arma::zeros(n,p+1);
    for (int i=0; i<n; i++){
      output.row(i) = spherical_exp(origin, tangents.row(i), 1.0);
    }
  } else if (aux_equalstr(geometry,"hyperbolic")){
    // an origin at [0,...,0,1]
    arma::rowvec origin(p+1,fill::zeros);
    origin(p) = 1.0;
    
    // tanget vectors
    arma::mat tangents = arma::randn(n,p+1,arma::distr_param(0.0, sd));
    tangents.col(p).fill(0.0);
    
    // embed
    output = arma::zeros(n,p+1);
    for (int i=0; i<n; i++){
      output.row(i) = hyperbolic_exp(origin, tangents.row(i), 1.0);
    }
  }
  
  // return
  return(output);
}