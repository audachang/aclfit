
distrn <- c("IGstod","exGAUS","norm","weibull","lnorm","logis.stats","gumbel.evd")

startflag <- c(T,T,F,F,F,T,T) #need starting flag or not

mlist <- c(1, 2, 3, 4, 5, 6, 7) #model index to be fitted

svals <- list(
    weibull = c(),
    exGAUS = list(mu = 500 , sigma=100, tau = 100),
    norm =c(),
    lnorm = c(),
    logis.stats = list(location=100, scale=300),
    invgauss.seq=list(kappa=156.62, xi=0.21, tau=152), #namely Wald distr.
    #exwald.ln.seq=list(kappa=156.62, xi=0.21, tau=152),
    IGstod = list(kappa=156.62, xi=0.21, tau=152), #revised statmod
    exwald.seq=list(kappa=156.62, xi=0.21, tau=152),
    gumbel.evd = list(loc=1500,scale=300)

)
lb <- list(
    weibull = c(-Inf, -Inf),
    exGAUS = c(-Inf,-Inf, -Inf),
    lnorm = c(-Inf, -Inf),
    logis = c(-Inf, -Inf),
    IGstod=c(-Inf,-Inf,-Inf),
    #exwald.ln.seq=c(-Inf,-Inf, -Inf),
    exwald.seq=c(-Inf,-Inf, -Inf)
)

ub <- list(
    weibull = c(Inf,Inf),
    exGAUS = c(Inf,Inf,Inf),
    lnorm = c(Inf,Inf),
    logis = c(Inf,Inf),
    IGstod=c(Inf,Inf,Inf), #invGauss/statmod pakckage.
    #exwald.ln.seq=c(Inf,Inf,Inf),
    exwald.seq=c(Inf,Inf,Inf)
)

lnflag <- c(F,F,F,F,F,F)

#mth <- c("mle","mle","mle","mle","mle","mle")
mth <- c("mle","mle","mle","mle","mle","mle","mle")

# optim.method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#            "Brent"),
oM <- list(
  weilbull = "Nelder-Mead",
  exGAUS = "SANN",
  lnorm = "Nelder-Mead",
  norm =  "Nelder-Mead",
  logis.stats = "Nelder-Mead",
  invgauss.seq="Nelder-Mead", #namely Wald distr.
  #exwald.ln.seq=list(kappa=156.62, xi=0.21, tau=152),
  IGstod = "Nelder-Mead", #revised statmod
  exwald.seq="SANN",
  gumbel.evd = "Nelder-Mead"  
  
)

qmeset <- which(mth =="qme")

pbs<- list(
    weibull = c(0.25,0.75),
    exGAUS = c(0.25,0.5,0.75),
    lnorm = c(),
    logis.stats = c(),
    exwald.seq = c(0.25,0.5,0.75),
    invgauss.seq = c(0.25,0.5,0.75)
)


# plot.legend <- c("weibull","exgaussian","lognormal",
#                  "logllogistic", "exwald.seq")



# plot.legend <- c("weibull","exgaussian","lognormal",
#                  "logllogistic", "shifted-wald","exwald")

ppflag <- T
qqflag <- T
cdflag <- T

#In Anders et al. (2016) https://doi.org/10.1037/met0000066, the shifted Wald distribution
#has three parameters: theta, alpha, and gamma 
#(Fig. 1, rightmost panel); 
#alternative parametization is: mu, lambda, and tau
#verification with the R & Stan implementation of Mr. Unadon
# https://mrunadon.github.io/Shifted-Wald-distribution-for-response-time-data-using-R-and-Stan/
#
# The seqmodels package implement the shifted Wald as IG with 3 parameters: kappa, xi, and tau
# All the above mentioned parameters of shifted Wald have the following relationship
# mu = alpha/gamma = kappa/xi (mu is mean in statmod)
# tau = theta: the shift parameter
# lambda = alpha^2 = 1/phi = kappa^2
# (phi in statmod; dispersion parameter; lambda is shape)
#
# invgauss in Seqmodels
#Random generation, density, distribution, and 
#quantile functions for the shifted inverse gaussian 
# (or Wald) distribution, parameterized for Brownian motion. 
# kappa refers to the threshold (alpha), 
# xi refers to the rate of evidence accumulation 
# towards this threshold (gamma) 
# tau is the shift in response times (theta) and
# sigma refers to the within-trial variability for 
# the rate of evidence accumulation 
# (the coefficient of drift, typically fixed to 1). 
# 
# The following are equal
# seqmodels::dinvgauss( 0.1, kappa = 0.3, xi = 2.0, tau = 0 )
# statmod::dinvgauss( 0.1, mean = 0.15, shape =  .09)
# Where kappa^2 = shape