#' Generalized Riemannian Multidimensional Scaling 
#' 
#' @param diss a S3 \code{\link[stats]{dist}} object of distance/dissimilarity among \eqn{n} objects.
#' @param p dimension of the target low-dimensional manifold (default: 2).
#' @param geometry name of the target geometry; one of \code{"euclidean"},\code{"spherical"}, and \code{"hyperbolic"} (default: \code{"euclidean"}).
#' @param ... extra parameters including \describe{
#' \item{main.iter}{maximum number of iterations for the main problem (default: 496).}
#' \item{main.tol}{stopping criterion for the main problem (default: \code{1e-5}).}
#' \item{sub.iter}{maximum number of iterations for subproblem (default: 10).}
#' \item{sub.tol}{stopping criterion for subproblems (default: \code{1e-5}.)}
#' \item{learning.fixed}{logical on how to solve subproblems. \code{TRUE} means that each iteration will proceed with a fixed step size and \code{FALSE} runs a line search, which may be expensive (default: \code{FALSE}).}
#' \item{learning.rate}{learning rate is a stepsize for subproblem optimization. It is only needed when \code{learning.fixed=TRUE} (default: 0.1).}
#' \item{weight.scheme}{exponent to control weighing scheme. It should be one of \eqn{0},\eqn{1}, or \eqn{2}. For more detailed explanation, see the section below (default: 0).}
#' \item{scale.spherical}{boolean to scale the \code{diss} for the spherical geometry case for better visualization (default: \code{TRUE}).}
#' \item{scale.factor}{scaling factor for \code{scale.spherical=TRUE} case. It should be a number in \eqn{(0,1)} (default: 0.99).}
#' }
#' 
#' @return a named list containing \describe{
#' \item{embedding}{low-dimensional embedding of \eqn{n} objects. If \code{geometry="euclidean"}, it returns an \eqn{(n\times p)} matrix. Otherwise, it returns an \eqn{(n\times (p+1))} matrix. In both cases, each row corresponds to an embedded vector onto the chosen manifold.}
#' }
#' 
#' @concept dissimilarity
#' @export
GERM <- function(diss, p=2, geometry=c("euclidean","spherical","hyperbolic"), ...){
  #=============================================================================
  # INPUT : EXPLICIT
  # diss
  if (!inherits(diss, "dist")){
    stop("* GERM : 'diss' should be a S3 'dist' class object.")
  }
  par_p    = max(1, round(p))
  par_geom = match.arg(geometry)

  # INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("learning.fixed"%in%pnames){
    exact_search = as.logical(params[["learning.fixed"]])
  } else {
    exact_search = TRUE
  }
  if ("learning.rate"%in%pnames){
    learning_rate = max(1e-10, as.double(params[["learning.rate"]]))
  } else {
    learning_rate = 0.1
  }
  if ("sub.iter"%in%pnames){
    iter_inner = max(2, round(as.numeric(params[["sub.iter"]])))
  } else {
    iter_inner = 10
  }
  if ("sub.tol"%in%pnames){
    tol_inner = max(10*.Machine$double.eps, as.double(params[["sub.tol"]]))
  } else {
    tol_inner = 1e-5
  }
  if ("main.iter"%in%pnames){
    iter_outer = max(3, round(as.numeric(params[["main.iter"]])))
  } else {
    iter_outer = 496
  }
  if ("main.tol"%in%pnames){
    tol_outer = max(10*.Machine$double.eps, as.double(params[["main.tol"]]))
  } else {
    tol_outer = 1e-5
  }
  if ("weight.scheme"%in%pnames){
    weight_order = round(as.numeric(params[["weight.scheme"]]))
    if (!(weight_order%in%c(0,1,2))){
      stop("* GERM : 'weight.scheme' should a numeric value of 0, 1, or 2.")
    }
  } else {
    weight_order = 0
  }
  if ("scale.spherical"%in%pnames){
    scale_sphere = as.logical(params[["scale.spherical"]])
  } else {
    scale_sphere = TRUE
  }
  if ("scale.factor"%in%pnames){
    scale_factor = as.double(params[["scale.factor"]])
    if ((scale_factor <= .Machine$double.eps)||(scale_factor >= (1-.Machine$double.eps))){
      stop("* GERM : 'scale.factor' should be a number in (0,1).")
    }
  } else {
    scale_factor = 0.99
  }

  #=============================================================================
  # COMPUTE
  # coerce the pairwise distance as a regular square matrix
  D = as.matrix(diss)

  # scaling for the spherical case
  if (aux_equalstr(par_geom, "spherical")){
    if (scale_sphere){
      D = D*((scale_factor*base::pi)/base::max(D))
    }
    if (any(D >= (2*base::pi))){
      stop("* GERM : 'diss' contains a value larger than 2*pi. Use 'scale.spherical=TRUE' or other geometries.")
    }
  }
  
  # run
  embedding = src_germ(D, par_p, par_geom, 
                       iter_outer, iter_inner, tol_outer, tol_inner, 
                       exact_search, learning_rate, weight_order)
  
  #=============================================================================
  # RETURN
  output = list()
  output$embedding = embedding
  return(output)
}

# X = as.matrix(iris[,1:4])
# dX = stats::dist(X)
# Y1 <- stats::cmdscale(dX)
# Y2 <- GERM(dX, main.tol=1e-5)$embedding
# par(mfrow=c(1,2))
# plot(Y1[,1], Y1[,2], pch=19, col=factor(iris[,5]), main="CMDS")
# plot(Y2[,1], Y2[,2], pch=19, col=factor(iris[,5]), main="MMDS")
# # 
# Z0 <- GERM(dX, main.tol=1e-5, weight.scheme=0)
# Z1 <- GERM(dX, main.tol=1e-5, weight.scheme=1)
# Z2 <- GERM(dX, main.tol=1e-5, weight.scheme=2)
# 
# par(mfrow=c(1,3))
# plot(Z0$embedding[,1], Z0$embedding[,2], pch=19, col=factor(iris[,5]), main="MMDS-weight=0")
# plot(Z1$embedding[,1], Z1$embedding[,2], pch=19, col=factor(iris[,5]), main="MMDS-weight=1")
# plot(Z2$embedding[,1], Z2$embedding[,2], pch=19, col=factor(iris[,5]), main="MMDS-weight=2")
# # 
# # G1 <- GERM(dX, main.tol=1e-4, geometry="euclidean")
# G1 <- cmdscale(dX, k=3)
# G2 <- GERM(dX, main.tol=1e-4, geometry="spherical")
# G3 <- GERM(dX, main.tol=1e-10, geometry="hyperbolic")
# colvec <- round(as.numeric(as.factor(iris[,5])))
# # 
# graphics.off()
# library(rgl)
# 
# target = G1
# plot3d(x=target[,1], y=target[,2], z=target[,3], col=colvec, size=10)
# target = G2$embedding
# plot3d(x=target[,1], y=target[,2], z=target[,3], col=colvec, size=10)
# target = G3$embedding
# plot3d(x=target[,1], y=target[,2], z=target[,3], col=colvec, size=5)
# plot3d(x=0, y=0, z=1, add=TRUE, col="cyan", size=20)
# 
# nn = base::nrow(G2$embedding)
# check_2 <- rep(0, nn)
# check_3 <- rep(0, nn)
# lorentz <- function(x){
#   return(sum(x[1:(length(x)-1)]^2)-(x[length(x)]^2))
# }
# for (i in 1:nn){
#   check_2[i] = sqrt(sum(G2$embedding[i,]^2))
#   check_3[i] = lorentz(G3$embedding[i,])
# }


# @param ... extra parameters including \describe{




# #  