# Script to fit a SEIR ODE  model with BayesianTools
# Fabienne Krauer, 15.08.2022
# fabienne.krauer@lshtm.ac.uk



# Household stuff  -------------------------------

library(deSolve)
library(tidyverse)
library(BayesianTools)
library(plyr)
set.seed(42)


# ODE Model -------------------------------

model <- function(times, inits, theta) {
  
  SEIR <- function(times, inits, theta) {
      
      S = inits[["S"]]
      E = inits[["E"]]
      I = inits[["I"]]
      R = inits[["R"]]
      N = S + E + I + R
    
      beta <- theta[["beta"]]
      sigma <- theta[["sigma"]]
      gamma <- theta[["gamma"]]
      
      dS <- -beta*S*I/N
      dE <- beta*S*I/N - sigma*E
      dI <- sigma*E - gamma*I
      dR <- gamma*I
      dC <- beta*S*I/N
      list(c(dS, dE, dI, dR, dC))
    
  }  
  
  traj <- data.frame(lsoda(inits, times, SEIR, theta))
  # Calculate the incidence per time step from the cumulative state:
  traj$inc <- c(inits["I"], diff(traj$C))
  return(traj)
  
}

# Parameters -------------------------------

beta <- 0.6
sigma <- 1/5 # =1/duration of latency
gamma <- 1/5 # = 1/duration of infectiousness

theta <- c(beta=beta, sigma=sigma, gamma=gamma)
Npop <- 10000
inits <- c("S"=Npop-1, "E"=0, "I"=1, "R"=0, "C"=1) 

# Time steps modelled ---------------------------------------------

# This vector defines the time points at which the ODE solution is saved.
# CAVE: the unit of the time steps MUST BE the same as the unit of your parameters
# i.e. if gamma = 1/5 days, then one time step (e.g. 1 to 2) is 1 day.
# Note: integration is computationally costly, save only the time steps you need. 

times <- seq(1:100)

# Run model once
traj <- model(times, inits, theta)

# Sanity check: Sum of all states must be constant and equivalent to Npop
summary(rowSums(traj[,c(2:5)]))


# Data -------------------------------

# Generate data from the trajectory assuming a Poisson process observation model
data <- data.frame(obs=sapply(traj$inc, function(x) rpois(1, x)),
                   time = times)

ggplot() + 
  geom_line(data=traj, aes(x=time, y=inc)) +
  geom_point(data=data, aes(x=time, y=obs))


# Prior -------------------------------

# Estimated params
estpars <- c("beta", "sigma") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params

# Priors
lower = c("beta"=1e-6, "sigma"=1e-6, "gamma"=1e-6)
upper = c("beta"=3, "sigma"=1, "gamma"=1)


# Create uniform prior
prior_unif <- createUniformPrior(lower=lower[estpars], 
                                 upper=upper[estpars])


# Likelihood function -------------------------------

# This is the loglik density function for the Poisson distribution
ll_Pois <- function(theta, inits, times, data) {
  
  traj <- model(times, inits, theta)
  datapoint <- data$obs
  modelpoint <- traj$inc

  if (any(is.na(modelpoint)) || any(modelpoint==0.0)) {
    ll <- -Inf
  } 
  
  
  else { # Minimize the log likelihood
    ll <- sum(dpois(x=datapoint, 
                    lambda=modelpoint,
                    log=TRUE), na.rm=TRUE)
  }
  
  return(ll)
}


# test
ll_Pois(theta, inits, times, data)


# Make a wrapper around the loglik function so that it is compatible how BT returns the estimated params
ll_Pois_wrapper <- function(par) {
  
  parX = theta
  parX[index] = par

  return(ll_Pois(theta=parX, 
                inits=inits, 
                times=times, 
                data=data))
}    

# Test
ll_Pois_wrapper(theta[index])



## Run DEMCMC ---------------------------------------------------

# MCMC settings
nchains <- 2
iter_per_chain <- 30000
mcmc_settings <- list(iterations = 3*iter_per_chain, 
                      nrChains = nchains)
sampler_algo <- "DEzs"

bayesianSetup_unif <- createBayesianSetup(prior = prior_unif,
                                     likelihood = ll_Pois_wrapper,
                                     names = names(theta[index]),
                                     parallel = FALSE)

system.time({chain_unif <- runMCMC(bayesianSetup = bayesianSetup_unif, 
                              sampler = sampler_algo, 
                              settings = mcmc_settings)})



# Diagnostics ------------------------------

plot(chain_unif)

# plot posterior after burn-in
nburn <- 10000
plot(chain_unif, parametersOnly = TRUE, start =nburn)

# check convergence
# PSRF gives the scale reduction factor for each parameter. 
# A PSRF of 1 means that variance between and within the chains are equal, which is the goal. 
# As a rule of thumb, a if PSRF < 1.1 or 1.05, convergence is very likely. But also check the traceplots!
gelmanDiagnostics(chain_unif, plot=TRUE, start=nburn)

# check parameter correlation
correlationPlot(chain_unif, start=nburn)

# check summary stats with burn-in
summary(chain_unif, start = nburn)




# Postprocessing --------------------------


# This function draws n samples from the posterior to calculate the posterior predictive uncertainty (95%)
trajsim_Pois <- function(chain, theta, inits, times, model, ndraw, nburn) {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(chain, parametersOnly = TRUE, thin=1, numSamples=ndraw, start=nburn)
  
  fit <- adply(.data=sample, .margins=1, .progress="text", .parallel=F, .fun=function(x) {
    
    #Define the theta as the estimated parameters from the current draw and the fixed parameters
    theta_sample <- c(x, theta[!is.na(theta)])
    
    #Simulate trajectory for selected theta
    foo <- match.fun(model)(times, inits, theta_sample)
    foo$simobs <- sapply(foo$inc, function(y) rpois(n=1, lambda=y))
    #foo$simobs <- sapply(foo$inc, function(y) rpois(n=1, lambda=y*theta_sample[["rho"]]))
    
    return(foo)
  })
  
  colnames(fit)[1] <- "replicate"
  return(fit)
  
} 


# Simulate with posteriors
fit <- trajsim_Pois(chain_unif, theta, inits, times, model, ndraw=999, nburn=nburn)
fit$replicate <- as.numeric(fit$replicate)

# Calculate quantiles of the trajectory
fit_quantiles <- plyr::ddply(.data=fit, 
                             .variables="time", 
                             function(x) quantile(x[, which(colnames(fit)=="simobs")], prob = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)) 
colnames(fit_quantiles) <- c("time", "low95", "low50", "median", "up50", "up95")

fit_quantiles <- merge(fit_quantiles, data[,c("time", "obs")], by="time", all=T)

# Plot fit
ggplot(fit_quantiles) +
  geom_point(aes(x=time, y=obs)) +
  geom_line(aes(x=time, y=median), color="red") +
  geom_ribbon(aes(x=time, ymin=low95, ymax=up95), alpha=0.2, fill="red")


# Example 2: Custom priors -----------------------

# Vectors that define the parameters of the prior distributions for each ODE model parameter
par1 <- c("beta"=0.5, "sigma"=2)
par2 <- c("beta"=0.1, "sigma"=5)
lower = c("beta"=1e-6, "sigma"=1e-6, "gamma"=1e-6)
upper = c("beta"=3, "sigma"=1, "gamma"=1)

# create a custom prior: beta as a NT() Prior, sigma as a Beta() prior
# Cave: ensure that the order of the parameters in the density function corresponds to the order in theta[index]
density <- function(par) {
  
  return(
    msm::dtnorm(par[1], # Beta 
                mean = par1[["beta"]], 
                sd = par2[["beta"]], 
                lower = lower[["beta"]], 
                upper = upper[["beta"]], 
                log = TRUE) + 
      dbeta(par[2],  # Sigma
            shape1 = par1[["sigma"]], 
            shape2 = par2[["sigma"]], 
            log = TRUE)
    
  )
}

sampler <-  function(n=1){
  
  return(cbind(
    
    msm::rtnorm(n, # Beta 
                mean = par1[["beta"]], 
                sd = par2[["beta"]], 
                lower = lower[["beta"]], 
                upper = upper[["beta"]]),
    rbeta(n, # Sigma
          shape1 =  par1[["sigma"]], 
          shape2 = par2[["sigma"]])
  ))
  
}

prior_custom <- createPrior(density=density, 
                            sampler=sampler, 
                            lower=lower, 
                            upper=upper)

# Fit again

bayesianSetup_custom <- createBayesianSetup(prior = prior_custom,
                                     likelihood = ll_Pois_wrapper,
                                     names = names(theta[index]),
                                     parallel = FALSE)

system.time({chain_custom <- runMCMC(bayesianSetup = bayesianSetup_custom, 
                              sampler = sampler_algo, 
                              settings = mcmc_settings)})

# Diagnostics
plot(chain_custom)

# plot posterior after burn-in
nburn <- 20000
plot(chain_custom, parametersOnly = TRUE, start =nburn)

# check convergence
gelmanDiagnostics(chain_custom, plot=TRUE, start=nburn)

# check parameter correlation
correlationPlot(chain_custom, start=nburn)

# check summary stats with burn-in
summary(chain_custom, start = nburn)

# Simulate with posteriors
fit <- trajsim_Pois(chain_custom, theta, inits, times, 
                    model, ndraw=999, nburn=nburn)
fit$replicate <- as.numeric(fit$replicate)

# Calculate quantiles of the trajectory
fit_quantiles <- plyr::ddply(.data=fit, 
                             .variables="time", 
                             function(x) quantile(x[, which(colnames(fit)=="simobs")], prob = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)) 
colnames(fit_quantiles) <- c("time", "low95", "low50", "median", "up50", "up95")

fit_quantiles <- merge(fit_quantiles, data[,c("time", "obs")], by="time", all=T)

# Plot fit
ggplot(fit_quantiles) +
  geom_point(aes(x=time, y=obs)) +
  geom_line(aes(x=time, y=median), color="red") +
  geom_ribbon(aes(x=time, ymin=low95, ymax=up95), alpha=0.2, fill="red")


########################################################################
# EXERCISES

# For all exercises, append the number of the exercise to your modified functions
# e.g. "prior1 <- ... " for a prior created in exercise 1 and "chain1 <- ..." for the chain


# Exercise 1: [~15 mins] ------------------------------------

# Modify the custom prior so that the 'beta' parameter follows a Gamma distribution with 
# a mean of 0.75 and a variance of 0.375 and fit the model again

# 

# Exercise 2: [~15 mins] ------------------------------------

# Modify the code to incorporate a constant underreporting. 
# Hint: Underreporting is often modelled as a fraction x of all cases being observed (with x [0,1])
# You can model it straightforward by scaling the simulated incidence with rho at each time step.

# a) Generate new data from the same ODE model but including a reported fraction "rho" = 0.4 (i.e. 40% reported)

# b) Modify the prior. # Use the prior function from the previous exercise. 
# Assume a Beta() prior distribution for rho with shape1=1 and shape2=1
# Note: you also need to modify the vectors that define the lower/upper bounds
# and the vectors with the parameters ("par1") for the prior distributions

# c) Modify the ll function to account for the underreporting

# d) Fit the model with the new ll function. 

# e) Inspect the chain and the model fit to the new data



# Exercise 3: [~30 mins]------------------------------------


# Fit the model with a Negative Binomial observation model instead of a Poisson.
# Use an additional parameter "k" to model the overdispersion.

# a) Read up how the Negative Binomial distribution is modelled in R.
# Decide which parametrization is most useful for count data. 
# Use the following vector "foo" and the sampling function and 
# explore different values for k to understand how it influences the variance
foo <- data.frame("data" = c(1,8,13,25,46,66,78,130,245,355,455,321,294,239,178,156,76,44,28,18,7,2),
                  "time" = 1:22)

sample_k <- function(k) {
  
  foo_obs <- sapply(foo$data, function(x) rnbinom(100, 
                                                  mu = x, 
                                                  size=k))
  
  quantiles <- t(apply(foo_obs, 2, quantile, probs=c(0.025, 0.975), names=TRUE))
  colnames(quantiles) <- c("low_95CI", "up_95CI")
  foo <- cbind(foo, quantiles)
  
  kplot <- ggplot(foo) + 
    geom_point(aes(x=time, y=data)) +
    geom_ribbon(aes(x=time, ymin=low_95CI, ymax=up_95CI), alpha=0.5, fill="red")
  
  return(kplot)
  
}

# What is the relationship between k and the variance in this parametrization?
# What happens if we use a slighly different parametrization of "size = 1/k"
# and why would this be preferable over "size = k"?


# b) Generate new data from the same model, but this time with a Negative Binomial process
# using the parametrization you think is best (see a)

# c) Modify the prior. Use the prior function from the previous exercise. 
# Assume a Beta() prior distribution for k with shape1=1 and shape2=1


# d) Adapt the ll function and the wrapper function. 
# Keep the underreporting from exercise 2.

# d) Fit the model with the new ll function

# e) Inspect the chain and the model fit to the new data






# Model comparison ----------------------------------

# Calculate the Bayes Factor of two models
M1 <- marginalLikelihood(chains_model1, start=nburn)
M2 <- marginalLikelihood(chains_model2, start=nburn)

exp(M1$ln.ML - M2$ln.ML)
# BF > 1 means the evidence is in favor of M1


# Parallel ----------------------------------------

cpus <- 2 # this may not be faster!


# Options for parallel computing
opts <- list(packages=list("BayesianTools", "deSolve"), 
             variables=list("theta", "data", "index", "times", "model", "inits","loglik"), 
             dlls=NULL)

bayesianSetup <- createBayesianSetup(prior = prior,
                                     likelihood = ll_wrapper,
                                     names = names(theta[index]),
                                     parallel = cpus,
                                     parallelOptions = opts)

system.time({chain <- runMCMC(bayesianSetup = bayesianSetup, 
                              sampler = sampler, 
                              settings = mcmc_settings)})

stopParallel(bayesianSetup)






# Extension --------------------

# The chain can be extended if it is not converged yet. CAVE: this does not work when parallelised
chain_ext <- runMCMC(bayesianSetup = chain, sampler = sampler_algo, settings = mcmc.settings)


