# aclfit
A R package wrapping fitting routines in ACL@NCU

fitdistplus is an excellent package for fitting distributions. Some typically used functions for reaction time analysis are defined soundly in other packages (seqmodels, statmod, or gamlss.dist), including error handling. However, fitdistplus does not like it when a error message is returned instead of NA. Therefore, wrapper functions are defined here for more smoothly using the functions such as Shifted Wald or exGaussian distributions. 
