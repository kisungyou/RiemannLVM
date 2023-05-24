# 
# # reference is manopt : https://github.com/NicolasBoumal/manopt/blob/master/manopt/manifolds/hyperbolic/hyperbolicfactory.m
# 
# rm(list=ls())
# 
# # personal understanding --------------------------------------------------
# hy_checker <- function(x){
#   return(sum(x[1:(length(x)-1)]^2)-x[length(x)]^2)
# }
# hy_inner <- function(x, y){
#   n = length(x)
#   return(sum(x[1:(n-1)]*y[1:(n-1)])-(x[n]*y[n]))
# }
# hy_dist1 <- function(u,v){
#   return(acosh(-hy_inner(u,v)))
# }
# hy_dist2 <- function(u,v){
#   z = u-v
#   return(2*asinh(0.5*sqrt(max(0, hy_inner(z,z)))))
# }
# hy_randvec <- function(d, sd=0.1){
#   x = rnorm(d, sd=sd)
#   y = sqrt(sum(x^2)+1)
#   return(c(x,y))
# }
# hy_normalize <- function(x){
#   return(x/sqrt(abs(hy_inner(x,x))))
# }
# 
# 
# # manopt operations -------------------------------------------------------
# M_inner <- function(x, u, v){
#   return(hy_inner(u,v))
# }
# M_proj <- function(x,u){
#   inners = hy_inner(x,u)
#   pu = u + x*inners
#   return(pu)
# }
# M_log <- function(x,y){
#   d = hy_dist2(x,y)
#   a = d/sinh(d)
#   a[is.na(a)]=1
#   u = M_proj(x, y*a)
#   return(u)
# }
# M_norm <- function(x,u){
#   return(sqrt(max(0, M_inner(x,u,u))))
# }
# M_exp <- function(x, u, t=1.0){
#   tu = t*u
#   mink_inner = hy_inner(tu,tu)
#   mink_norm  = sqrt(max(0, mink_inner))
#   
#   a = cosh(mink_norm)
#   b = sinh(mink_norm)/mink_norm
#   b[is.na(b)] = 1
#   
#   y = x*a + tu*b
#   return(y)
# }
# M_log2 <- function(mu, z){
#   alpha = -hy_inner(mu,z)
#   return((acosh(alpha)/sqrt((alpha^2) - 1))*(z-alpha*mu))
# }
# 
# 
# 
# 
# # example 1 : distance ----------------------------------------------------
# # random vectors
# x = hy_randvec(5)
# y = hy_randvec(5)
# 
# print(paste0("distance method 1: ", round(hy_dist1(x,y),3)))
# print(paste0("distance method 2: ", round(hy_dist2(x,y),3)))
# 
# 
# # example 2 : normalize ---------------------------------------------------
# z = 2*hy_randvec(5, sd=abs(rnorm(1, sd=0.05)))
# zz = hy_normalize(z)
# hy_checker(zz)
# 
# 
# # example 3 : norm of log vectors -----------------------------------------
# # random vectors
# x = hy_randvec(5, sd=0.01)
# y = hy_randvec(5, sd=0.01)
# 
# # log
# M_log(x,y)
# M_log2(x,y)
# 
# # log-exp : fine!
# M_exp(x, M_log(x,y))
# y
# 
# # gradient of the distance function
# u = hy_randvec(5)
# v = hy_randvec(5)
# -M_log(u,v)/hy_dist1(u,v)
# 
# gra-1/sqrt((hy_inner(u,v)^2)-1)



# hy_randvec <- function(d, sd=0.1){
#   x = rnorm(d, sd=sd)
#   y = sqrt(sum(x^2)+1)
#   return(c(x,y))
# }
# x1 = hy_randvec(5)
# x2 = hy_randvec(5)
# 
# hyperbolic_exp(x1, hyperbolic_log(x1,x2), 1)
