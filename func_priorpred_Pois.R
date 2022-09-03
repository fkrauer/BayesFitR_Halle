sample_prior_Pois <- function(prior, model, times, inits, 
                              theta, index, estpars, nsims) {
  
  traj <- vector("list", nsims)
  for (i in 1:nsims) {
    
    theta_est <- prior$sampler() # sample from the prior
    names(theta_est) <- estpars
    theta_updated <- c(theta[-index], theta_est)
    inc <- match.fun(model)(times, inits, theta_updated)[,c("inc")]
    inc_obs <- sapply(inc, function(x) rpois(1, x))
    traj[[i]] <- inc_obs
  }
  
  trajdf <- t(plyr::ldply(traj, rbind))
  quantiles <- data.frame(t(apply(trajdf, 1, quantile, prob = c(0.0, 1.0))))
  quantiles$time <- times
  colnames(quantiles)[1:2] <- c("min", "max")
  
  return(quantiles)
}
