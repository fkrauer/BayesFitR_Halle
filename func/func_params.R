# Params
beta <- 0.55
sigma <- 1/5 # =1/duration of latency
gamma <- 1/5 # = 1/duration of infectiousness
theta <- c(beta=beta, sigma=sigma, gamma=gamma)
inits <- c("S"=10000-1, "E"=0, "I"=1, "R"=0, "C"=1) 

times <- seq(1:100)