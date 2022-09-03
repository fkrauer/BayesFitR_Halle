ll_Pois <- function(model, theta, inits, times, data) {
  
  traj <- match.fun(model)(times, inits, theta)
  datapoint <- data$obs
  modelpoint <- traj$inc
  
  if (any(is.na(modelpoint))) {
    ll <- -Inf
  } 
  
  
  else { # Minimize the log likelihood
    ll <- sum(dpois(x=datapoint, 
                    lambda=modelpoint,
                    log=TRUE), na.rm=TRUE)
  }
  
  return(ll)
}

