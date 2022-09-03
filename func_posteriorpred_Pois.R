sample_posterior_Pois <- function(chain, 
                                  theta, 
                                  inits, 
                                  times, 
                                  model, 
                                  ndraw, 
                                  nburn, 
                                  progress="text") {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(chain, parametersOnly = TRUE, thin=1, numSamples=ndraw, start=nburn)
  
  fit <- plyr::adply(.data=sample, .margins=1, .progress=progress, .parallel=F, .fun=function(x) {
    
    #Define the theta as the estimated parameters from the current draw and the fixed parameters
    theta_sample <- c(x, theta[!is.na(theta)])
    
    #Simulate trajectory for selected theta
    foo <- match.fun(model)(times, inits, theta_sample)
    foo$simobs <- sapply(foo$inc, function(y) rpois(n=1, lambda=y))
    
    return(foo)
  })
  
  colnames(fit)[1] <- "replicate"
  
  quantiles <- plyr::ddply(.data=fit, 
                           .variables="time", 
                           function(x) quantile(x[, which(colnames(fit)=="simobs")], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
  colnames(quantiles) <- c("time", "low95PPI", "median", "up95PPI")
  
  return(quantiles)
  
} 
