# Script to fit a SEIR ODE  model with BayesianTools
# Author: Fabienne Krauer
# last updated: 15.08.2022
# contact: fabienne.krauer@lshtm.ac.uk


# For all exercises, append the number of the exercise to your modified functions or objects
# e.g. "prior1 <- ... " for a prior created in exercise 1 and "chain1 <- ..." for the chain

# NEVER NEVER adjust the window panes while you are fitting. 
# It could crash the R session ("the bomb") and you need to start from scratch.


# Household stuff  -------------------------------

library(deSolve)
library(tidyverse)
library(BayesianTools)
library(plyr)
set.seed(42)

# Prep   -------------------------------

# load models
source("func_models.R")

# load functions for prior and posterior predictive and loglik function
source("func_priorpred_Pois.R")
source("func_posteriorpred_Pois.R")
source("func_ll_Pois.R")



# Exercise 1:  -------------------------------

# Use the SEIR model and the parameters from the tutorial for this exercise: 
beta <- 0.6
sigma <- 1/5 # =1/duration of latency
gamma <- 1/5 # = 1/duration of infectiousness
theta <- c(beta=beta, sigma=sigma, gamma=gamma)
inits <- c("S"=10000-1, "E"=0, "I"=1, "R"=0, "C"=1) 

times <- seq(1:100)

traj <- model_SEIR(times, inits, theta)

data1 <- readRDS("data_ex1.rds") # load the synthetic data

# Define the parameters to be estimated
estpars <- c("beta", "sigma") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params

# Define the boundaries of the Priors (we don't need it for gamma, but do it anyway)
lower = c("beta"=1e-6, "sigma"=1e-6)
upper = c("beta"=5, "sigma"=1)


# Write a custom prior such that the 
# 'beta' parameter follows a Gamma distribution with a mean = 0.75 and a variance = 0.375. 
# For the 'sigma' parameter, use the prior from the tutorial
# HINT: check the wikipedia page to understand the properties of the Gamma distribution

density1 <- ...

sampler1 <- ...

# Setup the wrapper for the Loglik function:
ll_Pois_wrapper <- ...


# Use the same MCMC settings as in the tutorial:
nchains <- 2
iter_per_chain <- 20000
mcmc_settings <- list(iterations = iter_per_chain, 
                      nrChains = nchains)
sampler_algo <- "Metropolis"

## Run MCMC
bayesianSetup1 <- createBayesianSetup(prior = prior1,
                                      likelihood = ll_Pois_wrapper,
                                      names = names(theta[index]),
                                      parallel = FALSE)

system.time({chain1 <- runMCMC(bayesianSetup = bayesianSetup1, 
                               sampler = sampler_algo, 
                               settings = mcmc_settings)})

# If the chains are converged, check the posterior predictive (the fit)


# Exercise 2:  ------------------------------------

# Modify the code to incorporate a constant underreporting. 

# a) Generate new data ("data2") from the same ODE model but only 40% of all cases are reported
# Hint: Underreporting is often modelled as a fraction x of all cases being observed (with x [0,1])
# You can model it straightforward by scaling the simulated incidence with rho at each time step.

data2 <- ...

# b) Modify the prior. --> Use the prior function from the previous exercise as a template 
# Assume a Beta() prior distribution for rho with shape1=1 and shape2=1
# Note: you also need to modify the vectors that define the lower/upper bounds
# and the vector that defines the parameters to be estimated (estpar)

lower2 <- ...
upper2 <- ...

sampler2 <- ...
density2 <- ...

# c) Modify the ll function to account for the underreporting

ll_Pois2 <- 
ll_Pois2_wrapper <- 

# d) Fit the model with the new ll function. 
# Use the following settings for the fitting:
nchains2 <- 2
iter_per_chain2 <- 60000
sampler_algo2 <- "DEzs"

mcmc_settings2 <- list(iterations = iter_per_chain2, 
                      nrChains = nchains2)

bayesianSetup2 <- createBayesianSetup(prior = prior2,
                                      likelihood = ll_Pois2_wrapper,
                                      names = names(theta2[index2]),
                                      parallel = FALSE)

system.time({chain2 <- runMCMC(bayesianSetup = bayesianSetup2, 
                               sampler = sampler_algo2, 
                               settings = mcmc_settings2)})


# e) Inspect the chain and the model fit to the new data.
# You will need to adjust the function sample_posterior_Pois to account for the underreporting

sample_posterior_Pois2 <- ...

fit_quantiles2 <- ...

# Exercise 3: ------------------------------------


# Fit the model with a Negative Binomial observation model instead of a Poisson.
# Use an additional parameter "k" to model the overdispersion.

# a) Read up how the Negative Binomial distribution is modelled in R
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
# Plugin different values for k and observe what happens to the dispersion. 

# What happens if we use a slighly different parametrization of "size = 1/k"
# and why would this be preferable over "size = k"?


# b) Generate new data from the same SEIR model, but this time with a Negative Binomial process
# using the parametrization with k = 0.1 (which parametrization makes more sense for k = 0.1)?

# c) Modify the prior. Use the prior function from the previous exercise as a template 
# Assume a Beta() prior distribution for k with shape1=1 and shape2=1

# d) Adapt the loglik function and the wrapper function. 
# Keep the underreporting from exercise 2.

# d) Fit the model with the new loglik function
# Use the following MCMC settings
mcmc_settings3 <- list(iterations = 150000, 
                      nrChains = 2)

bayesianSetup3 <- createBayesianSetup(prior = prior3,
                                      likelihood = ll_NB3_wrapper,
                                      names = names(theta3[index3]),
                                      parallel = FALSE)

# This will take about 10 minutes to fit. You can start with Exercise 4 in the meanwhile
system.time({chain3 <- runMCMC(bayesianSetup = bayesianSetup3, 
                               sampler = sampler_algo2, 
                               settings = mcmc_settings3)})


# e) Inspect the chain and the model fit to the new data
# You will need to adjust the posterior predictive to account for the new loglik function

sample_posterior_NB3 <- ...

fit_quantiles3 <- ...

# f) Can you spot an issue with the chains and with some of the marginal parameters? 
# How can this potentially be improved?

par1 = c("beta"=1.5, "sigma"=15.0, "rho"=1.0, "k"=1.0) 
par2 = c("beta"=2.0, "sigma"=50.0, "rho"=1.0, "k"=1.0)

density3 <- function(par) {
  
  return(
    dgamma(par[1], # Beta 
           shape = par1[["beta"]], 
           rate =  par2[["beta"]], 
           log = TRUE) + 
      dbeta(par[2],  # Sigma
            shape1 = par1[["sigma"]], 
            shape2 = par2[["sigma"]], 
            log = TRUE) +
      dbeta(par[3], # rho
            shape1 = par1[["rho"]], 
            shape2 = par2[["rho"]], 
            log = TRUE) +
      dbeta(par[4], # rho
            shape1 = par1[["k"]], 
            shape2 = par2[["k"]], 
            log = TRUE) 
    
  )
}

sampler3 <-  function(n=1){
  
  return(cbind(
    
    rgamma(n,
           shape = par1[["beta"]], 
           rate = par2[["beta"]]), 
    rbeta(n, 
          shape1 =  par1[["sigma"]], 
          shape2 = par2[["sigma"]]),
    rbeta(n, 
          shape1 =  par1[["rho"]], 
          shape2 = par2[["rho"]]),
    rbeta(n, 
          shape1 =  par1[["k"]], 
          shape2 = par2[["k"]])
  ))
  
}

prior3 <- createPrior(density=density3, 
                      sampler=sampler3, 
                      lower=lower3, 
                      upper=upper3)








# Model comparison ----------------------------------

# Calculate the Bayes Factor of two models
M1 <- marginalLikelihood(chains_model1, start=nburn)
M2 <- marginalLikelihood(chains_model2, start=nburn)

exp(M1$ln.ML - M2$ln.ML)
# BF > 1 means the evidence is in favor of M1

