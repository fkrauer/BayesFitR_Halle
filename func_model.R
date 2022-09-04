model_SEIR <- function(times, inits, theta) {
  
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
