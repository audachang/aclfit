require(gamlss.dist)
require(seqmodels)
require(statmod)

dexGAUS <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE)
{

  if (any(sigma <= 0))
    #stop(paste("sigma must be greater than 0 ", "\n", ""))
    return(NaN)
  if (any(tau <= 0))
    #stop(paste("tau must be greater than 0 ", "\n", ""))
    return(NaN)
  ly <- length(x)
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)
  tau <- rep(tau, length = ly)
  z <- x - mu - ((sigma^2)/tau)
  logfy <- ifelse(tau > 0.05 * sigma,
                  -log(tau) - (z + (sigma^2/(2 *tau)))/tau +
                    log(pnorm(z/sigma)),
                  dnorm(x, mean = mu, sd = sigma,log = TRUE))
  if (log == FALSE)
    fy <- exp(logfy)
  else fy <- logfy
  fy
}

pexGAUS <- function (q, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(sigma <= 0))
    #stop(paste("sigma must be greater than 0 ", "\n", ""))
    return(NaN)
  if (any(tau <= 0))
    #stop(paste("tau must be greater than 0 ", "\n", ""))
    return(NaN)
  ly <- length(q)
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)
  tau <- rep(tau, length = ly)
  index <- seq(along = q)
  z <- q - mu - ((sigma^2)/tau)
  cdf <- ifelse(tau > 0.05 * sigma,
                pnorm((q - mu)/sigma) -
                  pnorm(z/sigma) *
                  exp(((mu + (sigma^2/tau))^2 - (mu^2) -
                         2 * q * ((sigma^2)/tau))/(2 * sigma^2)),
                pnorm(q, mean = mu, sd = sigma))

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

qexGAUS <- function (p, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)
{
  h1 <- function(q) {
    pexGAUS(q, mu = mu[i], sigma = sigma[i], tau = tau[i]) -
      p[i]
  }
  h <- function(q) {
    pexGAUS(q, mu = mu[i], sigma = sigma[i], tau = tau[i])
  }
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(tau <= 0))
    stop(paste("tau must be positive", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  lp <- max(length(p), length(mu), length(sigma), length(tau))
  p <- rep(p, length = lp)
  sigma <- rep(sigma, length = lp)
  mu <- rep(mu, length = lp)
  tau <- rep(tau, length = lp)
  q <- rep(0, lp)
  for (i in seq(along = p)) {
    if (h(mu[i]) < p[i]) {
      interval <- c(mu[i], mu[i] + sigma[i])
      j <- 2
      while (h(interval[2]) < p[i]) {
        interval[2] <- mu[i] + j * sigma[i]
        j <- j + 1
      }
    }
    else {
      interval <- c(mu[i] - sigma[i], mu[i])
      j <- 2
      while (h(interval[1]) > p[i]) {
        interval[1] <- mu[i] - j * sigma[i]
        j <- j + 1
      }
    }
    q[i] <- uniroot(h1, interval)$root
  }
  q
}

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



dinvgauss.seq <- function (x, kappa, xi, tau = as.numeric(c(0)), sigma = as.numeric(c(1)),
                           ln = FALSE)
{
  tryCatch({
    .Call("seqmodels_dinvgauss", PACKAGE = "seqmodels", x, kappa,
          xi, tau, sigma, ln)},
    error = function(err){
      return()
    })
}


pinvgauss.seq <- function (q, kappa, xi, tau = as.numeric(c(0)), sigma = as.numeric(c(1)),
                           bounds = 3, em_stop = 20, err = 1e-08)
{
  tryCatch({
    .Call("seqmodels_pinvgauss", PACKAGE = "seqmodels", q, kappa,
          xi, tau, sigma, bounds, em_stop, err)},
    error =function(err){
      return()
    }
  )
}

qinvgauss.seq <- function (p, kappa, xi, tau = as.numeric(c(0)), sigma = as.numeric(c(1)),
                           ln = FALSE, lower_tail = TRUE, ni = FALSE)
{

  tryCatch({
    .Call("seqmodels_qinvgauss", PACKAGE = "seqmodels", p, kappa,
          xi, tau, sigma, bounds, em_stop,err)},
    error = function(err){
      return()
    })
}



dexwald.seq <- function (x, kappa, xi, tau, sigma = as.numeric(c(1)), ln = FALSE,
                         ni = FALSE)
{
  tryCatch({
    #message(sprintf("kappa=%.3f,xi=%.3f,tau=%.3f",kappa,xi,tau))

    .Call("seqmodels_dexwald", PACKAGE = "seqmodels", x, kappa,
          xi, tau, sigma, ln, ni)},
    error = function(err){
      return()
    })
}


pexwald.seq <- function (q, kappa, xi, tau = as.numeric(c(0)), sigma = as.numeric(c(1)),
                         ln = FALSE, lower_tail = TRUE, ni = FALSE)
{

  tryCatch({
    .Call("seqmodels_pexwald", PACKAGE = "seqmodels", q, kappa,
          xi, tau, sigma, ln, lower_tail, ni)},
    error = function(err){
      return()
    })
}

qexwald.seq <- function (p, kappa, xi, tau = as.numeric(c(0)))
{
  sigma = as.numeric(c(1))

  bounds=.9
  em_stop=1e5
  err=9
  tryCatch({
    .Call("seqmodels_qexwald", PACKAGE = "seqmodels", p, kappa,
          xi, tau, sigma, bounds, em_stop,err)},
    error = function(err){
      return()
    })
}



dexwald.ln.seq <- function (x, kappa, xi, tau)
{
  #message(sprintf("kappa=%.3f,xi=%.3f,tau=%.3f"))
  sigma = as.numeric(c(1))
  ln = T
  ni = FALSE
  tryCatch({
    .Call("seqmodels_dexwald", PACKAGE = "seqmodels", x, kappa,
          xi, tau, sigma, ln, ni)},
    error = function(err){
      return()
    })
}


pexwald.ln.seq <- function (q, kappa, xi, tau = as.numeric(c(0)))
{
  sigma = as.numeric(c(1))
  ln = T
  lower_tail = TRUE
  ni = FALSE
  tryCatch({
    .Call("seqmodels_pexwald", PACKAGE = "seqmodels", q, kappa,
          xi, tau, sigma, ln, lower_tail, ni)},
    error = function(err){
      return()
    })
}

qexwald.ln.seq <- function (p, kappa, xi, tau = as.numeric(c(0)))
{
  sigma = as.numeric(c(1))

  bounds=.9
  em_stop=1e5
  err=9
  tryCatch({
    .Call("seqmodels_qexwald", PACKAGE = "seqmodels", p, kappa,
          xi, tau, sigma, bounds, em_stop,err)},
    error = function(err){
      return()
    })
}



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
