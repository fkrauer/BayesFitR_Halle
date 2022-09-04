# Halle Summer School
# Exercises to fit an ODE model with BayesianTools
# Author: Fabienne Krauer
# last updated: 04.09.2022
# contact: fabienne.krauer@lshtm.ac.uk


# Some notes:

# For all exercises, append the number of the exercise to your modified functions or objects
# e.g. "prior1 <- ... " for a prior created in exercise 1 and "chain1 <- ..." for the chain

# You should try to go through the exercises in ascending order.
# If you cannot figure one out, check and run the solutions before moving on to the next exercise.

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
source("func_model.R")

# load functions for prior and posterior predictive and loglik function
source("func_priorpred_Pois.R")
source("func_posteriorpred_Pois.R")
source("func_ll_Pois.R")
source("func_params.R")


# Exercise 1:  -------------------------------

# a) load and plot the  synthetic data
data1 <- readRDS("data_ex1.rds") 

# b) Adjust the prior:
# We will fit only to parameters for now:
estpars <- c("beta", "sigma")
index <- which(names(theta) %in% estpars) # index of estimated params
lower = c("beta"=1e-6, "sigma"=1e-6)
upper = c("beta"=5, "sigma"=1)

# Write a custom prior such that the 
# 'beta' parameter follows a Gamma distribution with a mean = 0.75 and a variance = 0.375. 
# For the 'sigma' parameter, use the prior from the tutorial
# HINT: check the wikipedia page to understand the properties of the Gamma distribution

density1 <- ...

sampler1 <- ...

# c) Setup the wrapper for the Loglik function:
ll_Pois_wrapper <- ...


# d) Use the same MCMC settings as in the tutorial:
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
# Assume that you don't know the underreporting and have to fit it. 
# For this, assume a Beta() distribution for rho prior with shape1=1 and shape2=1
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
# using the parametrization with k = 0.1 (which means that we use size =1/k in the NegBi function). 

# c) Modify the prior. Use the prior function from the previous exercise as a template 
# Assume that you don't know the dispersion parameter k and have to fit it.
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

# If you have some time left at the end of the exercises, 
# you can re-fit this model with an approach you think improves the posteriors



# Exercise 4: ------------------------------------

# For this exercise, you will use simulated data from model_SEIRS, but you
# will not know the true parameter values or how the data were generated.
# Use your knowledge from exercises 1-3 to fit a model of your choice to the data.

# To make sure you can fit a model in useful time, 
# you can make the following assumptions:

# 1. The disease is introduced in a completely susceptible population of size = 100000
# by 1 infected individual (thus S=100000-1)
# 2. The disease has a latent stage (E)
# 3. Individuals become temporarily immune for 180 days, after which they become
# suceptible again
# 4. Not all cases are reported, studies suggest that the reported fraction may be between 5 and 20%.
# 5. The latent period lasts for 6 days
# 6. The duration of infectiousness is estimated between 4 and 8 days
# 7. The data consist of daily new reported cases
# 8. The population size is constant, no deaths, births or migrations
# 9. The disease appears less transmissible than the one from exercise 1-3 

# We will add temporary (waning) immunity to the model:

# a) Read and plot the data
data4 <- readRDS("data_ex4.rds")

# b) Use model_SEIR and adjust the code to reflect the waning of immunity
model_SEIRS <- ...

# c) Set up a vector of inits and a vector of theta
inits4 <- ...
theta4 <- ...

# d) Set up an appropriate prior
estpars4 <- ...
index4 <- ...

lower4 <- ...
upper4 <- ...

density4 <- ...
sampler4 <- ...

# e) Choose a function for the observation process and set up an appropriate loglikelihood function and a wrapper for BT
# you can also use a function from any of the previous exercises, if you
# think it is the appropriate observation process

ll_4 <- ...
ll_4_wrapper <- ...

# f) fit the model. Here are some suggestions for settings:
# These settings will take about 20 minutes to fit
mcmc_settings4 <- list(iterations = 90000, 
                       nrChains = 2)

bayesianSetup4 <- createBayesianSetup(prior = prior4,
                                      likelihood = ll_NB4_wrapper,
                                      names = names(theta4[index4]),
                                      parallel = FALSE)

system.time({chain4 <- runMCMC(bayesianSetup = bayesianSetup4, 
                               sampler = "DEzs", 
                               settings = mcmc_settings4)})

