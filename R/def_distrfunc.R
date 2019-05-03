require(statmod)

dgumbel <- function(x, a, b) 1/b * exp((a - x)/b) * exp(-exp((a -+     x)/b)) 
pgumbel <- function(q, a, b) exp(-exp((a - q)/b)) 
qgumbel <- function(p, a, b) a - b * log(-log(p))

#In Anders et al. (2016) https://doi.org/10.1037/met0000066, the shifted Wald distribution
#has three parameters: theta, alpha, and gamma  (Fig. 1, rightmost panel); alternative
#parametization is: mu, lambda, and tau
#verification with the R & Stan implementation of Mr. Unadon
# https://mrunadon.github.io/Shifted-Wald-distribution-for-response-time-data-using-R-and-Stan/
#
# The seqmodels package implement the shifted Wald as IG with 3 parameters: kappa, xi, and tau
# All the above mentioned parameters of shifted Wald have the following relationship
# mu = alpha/gamma = kappa/xi
# tau = theta: the shift parameter
# lambda = alpha^2



#wrapper function revise inverse gaussian in statmod to 3 parameter
dIGstod <<- function(x, kappa, xi, tau, dispersion=1,
                    log =F) {
  
  mean <- kappa/xi
  shape <- kappa^2
  x.shift <- x - tau
  dval <- statmod::dinvgauss( x.shift, 
                              mean = mean, 
                              shape = shape,
                              dispersion = dispersion,
                              log=log )
  
  return(dval)

}
pIGstod <<- function(q, kappa, xi, tau, dispersion=1,
                    log.p =F, lower.tail=T) {
  
  mean <- kappa/xi
  shape <- kappa^2
  q.shift <- q - tau
  pval <- statmod::pinvgauss( q.shift, 
                              mean = mean, 
                              shape = shape,
                              dispersion = dispersion,
                              log.p = log.p,
                              lower.tail = lower.tail
                              )
  
  return(pval)
  
}

qIGstod <<- function(p, kappa, xi, tau, dispersion=1, 
                    lower.tail=TRUE, log.p=FALSE,
                    maxit=200L, tol=1e-14, trace=FALSE) {
  
  mean <- kappa/xi
  shape <- kappa^2
  
  qval <- statmod::qinvgauss( p, 
                              mean = mean, 
                              shape = shape,
                              dispersion = dispersion,
                              lower.tail = lower.tail,
                              log.p = log.p,
                              maxit = maxit,
                              tol = tol,
                              trace = trace)
  qval <- qval + tau
  
  return(qval)
  
}