#' Test Log-Exp Relationship
#' 
#' 
#' @keywords internal
#' @noRd
test_logexp <- function(d, sd=0.1, geometry=c("spherical","hyperbolic")){
  
  
  par_geom = match.arg(geometry)
  if (aux_equalstr(par_geom, "spherical")){
    x = rnorm(d, sd=sd)
    y = rnorm(d, sd=sd)
    
    x = x/sqrt(sum(x^2))
    y = y/sqrt(sum(y^2))
    
    yy = spherical_exp(x, spherical_log(x,y), 1)
  } else {
    x = hy_randvec(d, sd)
    y = hy_randvec(d, sd)
    
    yy = hyperbolic_exp(x, hyperbolic_log(x,y), 1)
  }
  
  output = list()
  output$true = y
  output$recon = yy
  return(output)
  
}

#' @keywords internal
#' @noRd
hy_randvec <- function(d, sd=0.1){
  x = rnorm(d, sd=sd)
  y = sqrt(sum(x^2)+1)
  return(c(x,y))
}

