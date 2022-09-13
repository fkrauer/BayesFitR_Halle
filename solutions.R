# Halle Summer School
# Solutions to the exercises
# Author: Fabienne Krauer
# last updated: 04.09.2022
# contact: fabienne.krauer@lshtm.ac.uk


# Exercise 1 -------------------------------------------------------------------

#b)
# The Gamma distribution in R is modelled with dgamma() and rgamma(). 
# The function allows two different parametrizations: 
# Gamma(shape, rate) or Gamma(shape, scale). 

# Assuming a mean (E(X)) of 0.75 and a variance (Var(X)) of 0.375, we first calculate the
# shape/rate OR shape/scale parameters of the Gamma distribution. It doesn't matter which you 
# choose, just make sure you use the same parametrization for the prior.
# The shape/rate parametrization is more common in Bayesian statistics (where shape=alpha, rate=beta)

# Since mean = alpha/beta and var = alpha/beta^2, we can plugin the requested values and solve
# this system of two equations:

# alpha/beta = 0.75 
# -> alpha = 0.75*beta 
# -> (0.75*beta)/beta^2 = 0.375 
# -> 0.75/beta = 0.375 
# -> beta = 2 
# -> alpha/2 = 0.75 
# -> alpha = 1.5
# Result: alpha=shape=1.5 and beta=rate=2

# The adapted prior becomes:
par1 = c("beta"=1.5, "sigma"=2) 
par2 = c("beta"=2.0, "sigma"=5)


density1 <- function(par) {
  
  return(
    dgamma(par[1], # Beta 
                shape = par1[["beta"]], 
                rate =  par2[["beta"]], 
                log = TRUE) + 
      dbeta(par[2],  # Sigma
            shape1 = par1[["sigma"]], 
            shape2 = par2[["sigma"]], 
            log = TRUE) 
    
  )
}

sampler1 <-  function(n=1){
  
  return(cbind(
    
    rgamma(n,
           shape = par1[["beta"]], 
           rate = par2[["beta"]]), 
    rbeta(n, 
          shape1 =  par1[["sigma"]], 
          shape2 = par2[["sigma"]])
  ))
  
}

prior1 <- createPrior(density=density1, 
                             sampler=sampler1, 
                             lower=lower, 
                             upper=upper)

# C) 
# Setup the wrapper
ll_Pois_wrapper <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(ll_Pois(model=model_SEIR,
                 theta=parX, 
                 inits=inits, 
                 times=times, 
                 data=data1))
}    

# Test (should return the same as ll_pois())
ll_Pois_wrapper(theta[index])

# d) Check the fit
fit_quantiles1 <- sample_posterior_Pois(chain1, theta, inits, times, model_SEIR, ndraw=500, nburn=10000, progress="text")

# Plot fit
ggplot() +
  geom_point(data=data1, aes(x=time, y=obs)) +
  geom_line(data=fit_quantiles1, aes(x=time, y=median), color="red") +
  geom_ribbon(data=fit_quantiles1, aes(x=time, ymin=low95PPI, ymax=up95PPI), alpha=0.2, fill="red")


# Exercise 2 -------------------------------------------------------------------

# a) Generate new data

# Update theta
rho <- 0.4
theta2 <- c(beta=beta, sigma=sigma, gamma=gamma, rho=rho)

data2 <- data.frame(obs=sapply(traj1$inc, function(x) rpois(1, x * theta2[["rho"]])),
                    time = times)

ggplot() + 
  geom_line(data=traj1, aes(x=time, y=inc)) +
  geom_point(data=data2, aes(x=time, y=obs))


# b) Modify prior

lower2 = c("beta"=0, "sigma"=0, "rho"=0)
upper2 = c("beta"=Inf, "sigma"=1.0, "rho"=1.0)
par1 = c("beta"=1.5, "sigma"=2.0, "rho"=1.0) 
par2 = c("beta"=2.0, "sigma"=5.0, "rho"=1.0)

estpars2 <- c("beta", "sigma", "rho") # parameters to estimate, can be modified
index2 <- which(names(theta2) %in% estpars2) # index of estimated params
theta2[index2]

density2 <- function(par) {
  
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
            log = TRUE) 
    
  )
}

sampler2 <-  function(n=1){
  
  return(cbind(
    
    rgamma(n,
           shape = par1[["beta"]], 
           rate = par2[["beta"]]), 
    rbeta(n, 
          shape1 =  par1[["sigma"]], 
          shape2 = par2[["sigma"]]),
    rbeta(n, 
          shape1 =  par1[["rho"]], 
          shape2 = par2[["rho"]])
  ))
  
}

prior2 <- createPrior(density=density2, 
                      sampler=sampler2, 
                      lower=lower2, 
                      upper=upper2)


# c) Modify the ll function

ll_Pois2 <- function(model, theta, inits, times, data) {
  
  traj <- match.fun(model)(times, inits, theta)
  datapoint <- data$obs
  modelpoint <- traj$inc
  
  if (any(is.na(modelpoint))) {
    ll <- -Inf
  } 
  
  
  else { # Minimize the log likelihood
    ll <- sum(dpois(x=datapoint, 
                    lambda=modelpoint * theta[["rho"]],
                    log=TRUE), na.rm=TRUE)
  }
  
  return(ll)
}


# test
ll_Pois2(model_SEIR, theta2, inits, times, data2)


# Make a wrapper around the loglik function so that it is compatible how BT returns the estimated params
ll_Pois2_wrapper <- function(par) {
  
  parX = theta2
  parX[index2] = par
  
  return(ll_Pois2(model=model_SEIR,
                  theta=parX, 
                 inits=inits, 
                 times=times, 
                 data=data2))
}    

# Test
ll_Pois2_wrapper(theta2[index2])

# e) Inspect the chain and the model fit to the new data
plot(chain2)


# Update the trajsim function to incorporate the underreporting
sample_posterior_Pois2 <- function(chain, 
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
    foo$simobs <- sapply(foo$inc, function(y) rpois(n=1, lambda=y*theta_sample[["rho"]]))
    
    
    return(foo)
  })
  
  quantiles <- plyr::ddply(.data=fit, 
                           .variables="time", 
                           function(x) quantile(x[, which(colnames(fit)=="simobs")], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
  colnames(quantiles) <- c("time", "low95PPI", "median", "up95PPI")
  
  return(quantiles)
  
} 


fit_quantiles2 <- sample_posterior_Pois2(chain2, theta2, inits, times, model_SEIR, ndraw=500, nburn=10000, progress="none")

# Plot fit
ggplot() +
  geom_point(data=data2, aes(x=time, y=obs)) +
  geom_line(data=fit_quantiles2, aes(x=time, y=median), color="red") +
  geom_ribbon(data=fit_quantiles2, aes(x=time, ymin=low95PPI, ymax=up95PPI), alpha=0.2, fill="red")


# Exercise 3 -------------------------------------------------------------------

# The Negative Binomial in R is modelled with dnbinom() / rnbinom().
# It has two different parametrizations:
# (size/prob) for the classic binomial trial/successes parametrization and
# (size/mu) for when the observations are positive integers, 
# with size being the overdispersion parameter. It must be strictly positive 
# We commonly use the second parametrization for modelling count data. 


# a) Sample from foo and explore how k changes the variance
sample_k(0.1)
sample_k(1)
sample_k(10)
sample_k(100)

# --> the larger k, the smaller the variance

sample_k_inv <- function(k) {
  
  foo_obs <- sapply(foo$data, function(x) rnbinom(100, 
                                                  mu = x, 
                                                  size=1/k))
  
  quantiles <- t(apply(foo_obs, 2, quantile, probs=c(0.025, 0.975), names=TRUE))
  colnames(quantiles) <- c("low_95CI", "up_95CI")
  foo <- cbind(foo, quantiles)
  
  kplot <- ggplot(foo) + 
    geom_point(aes(x=time, y=data)) +
    geom_ribbon(aes(x=time, ymin=low_95CI, ymax=up_95CI), alpha=0.5, fill="red")
  
  return(kplot)
  
}

sample_k_inv(0.01)
sample_k_inv(0.1)
sample_k_inv(1)
sample_k_inv(10)
sample_k_inv(100)

# --> the smaller k, the smaller the variance. Values above 1 make little difference in the variance
# --> this property is useful because we can define a finite uniform prior [0,1] rather than [0, Inf]


# b) Generate new data
k = 0.1
theta3 <- c(beta=beta, sigma=sigma, gamma=gamma, rho=rho, k=k)

traj3 <- model(times, inits, theta)
data3 <- data.frame(obs=sapply(traj3$inc, function(x) rnbinom(1, 
                                                             mu = x*theta3[["rho"]], 
                                                             size=1/theta3[["k"]])),
                    time = times)

ggplot() + 
  geom_line(data=traj3, aes(x=time, y=inc)) +
  geom_point(data=data3, aes(x=time, y=obs))


# c) Modify the prior

lower3 = c("beta"=0, "sigma"=0, "rho"=0, "k"=0.0)
upper3 = c("beta"=5.0, "sigma"=1.0, "rho"=1.0, "k"=1.0)
par1 = c("beta"=1.5, "sigma"=2.0, "rho"=1.0, "k"=1.0) 
par2 = c("beta"=2.0, "sigma"=5.0, "rho"=1.0, "k"=1.0)

estpars3 <- c("beta", "sigma", "rho", "k") # parameters to estimate, can be modified
index3 <- which(names(theta3) %in% estpars3) # index of estimated params
theta3[index3]

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


# d) Adapt the ll function and the wrapper function

ll_NB3 <- function(model, theta, inits, times, data) {
  
  traj <- match.fun(model)(times, inits, theta)
  datapoint <- data$obs
  modelpoint <- traj$inc
  
  if (any(is.na(modelpoint))) {
    ll <- -Inf
  } 
  
  
  else { # Minimize the log likelihood
    ll <- sum(dnbinom(x=datapoint,
                      mu=modelpoint * theta[["rho"]],
                      size=1/theta[["k"]],
                      log=TRUE), na.rm=TRUE)
  }
  
  return(ll)
}

# test
ll_NB3(model_SEIR, theta3, inits, times, data3)

# Wrapper
ll_NB3_wrapper <- function(par) {
  
  parX = theta3
  parX[index3] = par
  
  return(ll_NB3(model=model_SEIR,
                theta=parX, 
                  inits=inits, 
                  times=times, 
                  data=data3))
}    

# Test
ll_NB3_wrapper(theta3[index3])

# e) Fit

# f) Assess diagnostics and fit
plot(chain3)

nburn = 20000

plot(chain3, parametersOnly = TRUE, start =nburn)

gelmanDiagnostics(chain3, plot=TRUE, start=nburn)

correlationPlot(chain3, start=nburn)

# Modify trajsim function to accommodate NegBin reporting
sample_posterior_NB3 <- function(chain, 
                                 theta, 
                                 inits, 
                                 times, 
                                 model, 
                                 ndraw, 
                                 nburn, 
                                 progress="text") {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(chain, parametersOnly = TRUE, thin=1, numSamples=ndraw, start=nburn)
  
  fit <- adply(.data=sample, .margins=1, .progress=progress, .parallel=F, .fun=function(x) {
    
    #Define the theta as the estimated parameters from the current draw and the fixed parameters
    theta_sample <- c(x, theta[!is.na(theta)])
    
    #Simulate trajectory for selected theta
    foo <- match.fun(model)(times, inits, theta_sample)
    foo$simobs <- sapply(foo$inc, function(y) rnbinom(1, 
                                                      mu = y*theta[["rho"]], 
                                                      size=1/theta[["k"]]))

    return(foo)
  })
  
  quantiles <- plyr::ddply(.data=fit, 
                           .variables="time", 
                           function(x) quantile(x[, which(colnames(fit)=="simobs")], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
  colnames(quantiles) <- c("time", "low95PPI", "median", "up95PPI")
  
  return(quantiles)
  
  
} 

# Assess fit
fit_quantiles3 <- sample_posterior_NB3(chain3, theta3, inits, times, model_SEIR, ndraw=500, nburn=10000, progress="none")

# Plot fit
ggplot() +
  geom_point(data=data3, aes(x=time, y=obs)) +
  geom_line(data=fit_quantiles3, aes(x=time, y=median), color="red") +
  geom_ribbon(data=fit_quantiles3, aes(x=time, ymin=low95PPI, ymax=up95PPI), alpha=0.2, fill="red")


# f) beta and sigma are strongly correlated (as before), but now 
# beta and sigma are also very wide and not very informative. 
# This could potentially be improved by tightening the priors for beta and sigma, e.g.

# beta ~ Gamma(4.5, 6.0)
# sigma ~ Beta(15, 50)

# If you have some time left at the end of the exercises, 
# you can re-fit this model with tighter priors


# Exercise 4: ------------------------------------


model_SIR_age <- function(times, inits, parms) {
  
  SIR_age <- function(times, state, parms) {
    
    N <- (state[(0*parms[["n.age"]]+1):(1*parms[["n.age"]])] + state[(1*parms[["n.age"]]+1):(2*parms[["n.age"]])] + state[(2*parms[["n.age"]]+1):(3*parms[["n.age"]])])
    
    flow_I <- (state[(0*parms[["n.age"]]+1):(1*parms[["n.age"]])] * parms[["beta"]]) * (parms[["contacts"]] %*% ((state[(1*parms[["n.age"]]+1):(2*parms[["n.age"]])]))/N)
    
    dS_N <- -flow_I + parms[["sigma"]] * state[(2*parms[["n.age"]]+1):(3*parms[["n.age"]])]
    dI_N <-  flow_I                                                                           - parms[["gamma"]] * state[(1*parms[["n.age"]]+1):(2*parms[["n.age"]])]
    dR_N <-           - parms[["sigma"]] * state[(2*parms[["n.age"]]+1):(3*parms[["n.age"]])] + parms[["gamma"]] * state[(1*parms[["n.age"]]+1):(2*parms[["n.age"]])] 
    
    list(c(dS_N, dI_N, dR_N, flow_I))
    
  }  
  
  traj <- data.frame(lsoda(y = inits, times = times, func = SIR_age, parms = params))
  # Calculate the incidence per time step from the cumulative state:
  # This is very hacky, can be improved
  traj$inc1 <- c(inits["flow_I1"], diff(traj$flow_I1))
  traj$inc2 <- c(inits["flow_I2"], diff(traj$flow_I2))
  traj$inc3 <- c(inits["flow_I3"], diff(traj$flow_I3))
  traj$inc4 <- c(inits["flow_I4"], diff(traj$flow_I4))
  
  # reshape to long to match data
  traj <- traj %>% select(time, inc1, inc2, inc3, inc4) %>% 
              pivot_longer(traj, cols=c(2:5),
                                names_to = "compartment",
                                values_to = "cases") %>% 
              arrange(time, compartment)
  
  return(traj)
  
}


traj4 <- model_SIR_age(times4, inits4, parms)
ggplot(traj4) +
  geom_line(aes(x=time, y=cases, colour=as.factor(compartment))) +
  geom_point(data=data4, aes(x=DAY, y=CASES, colour=as.factor(AGE)))

# LL
ll_Pois4 <- function(model, theta, parms, inits, times, data) {
  
  parms$gamma <- theta
  
  traj <- match.fun(model)(times, inits, parms)
  datapoint <- data$CASES
  modelpoint <- traj$cases
  
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

ll_Pois4(model_SIR_age, params$theta, params, inits, times, data4)

# LL wrapper
ll_Pois4_wrapper <- function(par) {
  
  parX = theta4
  parX[index4] = par
  
  return(ll_Pois4(model=model_SIR_age,
                  parms=params,
                  theta=parX, 
                  inits=inits4, 
                  times=times4, 
                  data=data4))
}    

# Test
ll_Pois4_wrapper(theta4[index4])


# Priors
lower4 <- c(0.0, 0.0, 0.0, 0.0)
upper4 <- c(0.6, 0.6, 0.6, 0.6)
names(lower4) <- names(upper4) <- names(theta4)

prior4 <- createUniformPrior(lower=lower4, 
                                 upper=upper4)

nchains4 <- 2
iter_per_chain4 <- 60000
sampler_algo4 <- "DEzs"

mcmc_settings4 <- list(iterations = iter_per_chain4, 
                       nrChains = nchains4)

bayesianSetup4 <- createBayesianSetup(prior = prior4,
                                      likelihood = ll_Pois4_wrapper,
                                      names = names(theta4[index4]),
                                      parallel = FALSE)

system.time({chain4 <- runMCMC(bayesianSetup = bayesianSetup4, 
                               sampler = sampler_algo4, 
                               settings = mcmc_settings4)})



# Exercise 5: ------------------------------------


model_SEIRS <- function(times, inits, theta) {
  
  SEIR <- function(times, inits, theta) {
    
    S = inits[["S"]]
    E = inits[["E"]]
    I = inits[["I"]]
    R = inits[["R"]]
    N = S + E + I + R
    
    beta <- theta[["beta"]]
    sigma <- theta[["sigma"]]
    gamma <- theta[["gamma"]]
    omega <- theta[["omega"]]
    
    dS <- -beta*S*I/N + omega*R
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- gamma*I - omega*R
    dC <- beta*S*I/N
    list(c(dS, dE, dI, dR, dC))
    
  }  
  
  traj <- data.frame(lsoda(inits, times, SEIR, theta))
  # Calculate the incidence per time step from the cumulative state:
  traj$inc <- c(inits["I"], diff(traj$C))
  return(traj)
  
}


# The data were generated as follows:
beta <- 0.2
omega <- 1/180
rho <- 0.12
sigma <- 1/6
gamma <- 1/7
k <- 0.05

theta4 <- c(beta=beta, sigma=sigma, gamma=gamma, omega=omega, rho=rho, k=k)
inits4 <- c("S"=100000-1, "E"=0, "I"=1, "R"=0, "C"=1) 

times4 <- seq(1:1000)

traj4 <- model_SEIRS(times4, inits4, theta4)

ggplot(traj4) +
  geom_line(aes(x=time, y=inc))

data4 <- data.frame(obs=sapply(traj4$inc, function(x) rnbinom(1, 
                                                              mu = x * rho, 
                                                              size=1/k)),
                   time = times4)
ggplot() + 
  geom_line(data=traj4, aes(x=time, y=inc)) +
  geom_point(data=data4, aes(x=time, y=obs))

saveRDS(data4, "data_ex4.rds")


# Set up the parameter vector and priors
estpars4 <- c("beta", "gamma", "rho", "k") # parameters to estimate, can be modified
index4 <- which(names(theta4) %in% estpars4) # index of estimated params
theta4[index4]

# We will fit a uniform prior for rho, and beta priors for the rest
lower4 = c("beta"=0,  "gamma"=0, "rho"=0.05, "k"=0.0)
upper4 = c("beta"=1.0,  "gamma" = 1.0, "rho"=0.2, "k"=1.0)
par1 = c("beta"=1.5,  "gamma"=2, "k"=1.0) 
par2 = c("beta"=4.0,  "gamma"=9, "k"=50)


density4 <- function(par) {
  
  return(
    dbeta(par[1], # Beta 
          shape1 = par1[["beta"]], 
          shape2 =  par2[["beta"]], 
           log = TRUE) +
      
      dbeta(par[2], # gamma
            shape1 = par1[["gamma"]], 
            shape2 = par2[["gamma"]], 
            log = TRUE) +
      
      dunif(par[3], # rho
            min = lower4[["rho"]], 
            max = upper4[["rho"]], 
            log = TRUE) +
  
      dbeta(par[4], # k
            shape1 = par1[["k"]], 
            shape2 = par2[["k"]], 
            log = TRUE) 
    
  )
}

sampler4 <-  function(n=1){
  
  return(cbind(
    
    rbeta(n,
          shape1 = par1[["beta"]], 
          shape2 = par2[["beta"]]), 
    
    rbeta(n, 
          shape1 =  par1[["gamma"]], 
          shape2 = par2[["gamma"]]),
    
    runif(n, 
          min = lower4[["rho"]], 
          max = upper4[["rho"]]),
    
    rbeta(n, 
          shape1 =  par1[["k"]], 
          shape2 = par2[["k"]])
  ))
  
}

prior4 <- createPrior(density=density4, 
                      sampler=sampler4, 
                      lower=lower4, 
                      upper=upper4)


# We will use the loglik NB function from exercise 3:
ll_NB3(model_SEIRS, theta4, inits4, times4, data4)

# Update the wrapper
ll_NB4_wrapper <- function(par) {
  
  parX = theta4
  parX[index4] = par
  
  return(ll_NB3(model=model_SEIRS,
                theta=parX, 
                inits=inits4, 
                times=times4, 
                data=data4))
}    

# Test
ll_NB4_wrapper(theta4[index4])

# Setup the MCMC and run
mcmc_settings4 <- list(iterations = 90000, 
                       nrChains = 2)

bayesianSetup4 <- createBayesianSetup(prior = prior4,
                                          likelihood = ll_NB4_wrapper,
                                          names = names(theta4[index4]),
                                          parallel = FALSE)

system.time({chain4 <- runMCMC(bayesianSetup = bayesianSetup4, 
                                    sampler = "DEzs", 
                                    settings = mcmc_settings4)})



