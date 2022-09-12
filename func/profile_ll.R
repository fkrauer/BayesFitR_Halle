library(deSolve)
library(tidyverse)
library(plotly)
set.seed(42)

theme_set(theme_classic())

# Prep   -------------------------------------


source("func_model.R")
source("func_ll_Pois.R")
source("func_params.R")
data <- readRDS("data_ex1.rds")


# plot the likelihood profile
calc_ll <- function(x, y, z){
    theta[["beta"]] <- x
    theta[["sigma"]] <- y
    theta["gamma"] <- z
    return(ll_Pois(model_SEIR, theta, inits, times, data))
}


# Simulate   -------------------------------------

# Simulate all parameters in a given range
range <- seq(0.05, 1, by=.05)

profile_LL <- matrix(NA, nrow=length(range)^length(theta), ncol=length(theta)+1)
count <- 1
for (i in range) {
  for (j in range) {
    for (k in range) {
      profile_LL[count,1] <- i
      profile_LL[count,2] <- j
      profile_LL[count,3] <- k
      profile_LL[count,4] <- calc_ll(i,j,k)
      count <- count+1
    }
  }
}

profile_LL <- data.frame(profile_LL)
colnames(profile_LL) <- c("beta", "sigma", "gamma", "ll")
profile_LL <- profile_LL %>% filter(!is.infinite(ll))


# 1D profile likelihood --------------------------------

profile_LL %>% 
  filter(sigma==theta[["sigma"]] & gamma==theta[["gamma"]]) %>% 
ggplot(.) +
  geom_line(aes(x=beta, y=ll)) +
  geom_vline(aes(xintercept=theta[["beta"]]), colour="red") +
  ylab("log likelihood")

profile_LL %>% 
  filter(beta==theta[["beta"]] & gamma==theta[["gamma"]]) %>% 
  ggplot(.) +
  geom_line(aes(x=sigma, y=ll)) +
  geom_vline(aes(xintercept=theta[["sigma"]]), colour="red") +
  ylab("log likelihood")

profile_LL %>% filter(sigma==theta[["sigma"]] & beta==theta[["beta"]]) %>% 
  ggplot(.) +
  geom_line(aes(x=gamma, y=ll)) +
  geom_vline(aes(xintercept=theta[["gamma"]]), colour="red") +
  ylab("log likelihood")


# 2D profile likelihood --------------------------------

# Beta vs. Sigma

#2D plot
ggplot(profile_LL[profile_LL$gamma==0.2,]) +
  geom_contour_filled(aes(x=beta, y=sigma, z=ll), bins=50) +
  theme(legend.position = "none") +
  geom_point(aes(x=theta[["beta"]], y=theta[["sigma"]]), colour="red")


# 3D plot
z <- profile_LL %>% filter(gamma==0.2) %>% 
                select(beta, sigma, ll) %>% 
              pivot_wider(id_cols = beta, names_from=sigma, values_from = ll) %>% 
              select(-beta) %>% as.matrix() 

landscape3D <- plot_ly(z = ~ z, 
                       x = ~ range, 
                       y = ~ range,
                       type="surface") %>% 
              hide_colorbar() %>% 
              add_trace(x=sigma,
                        y=beta,
                        z=max(z),
                        mode='markers',
                        type = "scatter3d",
                        marker= list(color="red")) %>% 
              layout(scene = list(xaxis=list(title="sigma"),
                                  yaxis=list(title="beta"), 
                                  zaxis=list(title="log likelihood"))) 

landscape3D


# Beta vs. Gamma

#2D plot
ggplot(profile_LL[profile_LL$sigma==0.2,]) +
  geom_contour_filled(aes(x=beta, y=gamma, z=ll), bins=50) +
  theme(legend.position = "none") +
  geom_point(aes(x=theta[["beta"]], y=theta[["gamma"]]), colour="red")


# 3D plot
z <- profile_LL %>% filter(sigma==0.2) %>% 
  select(beta, gamma, ll) %>% 
  pivot_wider(id_cols = beta, names_from=gamma, values_from = ll) %>% 
  select(-beta) %>% as.matrix()

landscape3D <- plot_ly(z = ~ z, 
                       x = ~ range, 
                       y = ~ range,
                       type="surface") %>% 
  hide_colorbar() %>% 
  add_trace(x=gamma,
            y=beta,
            z=max(z),
            mode='markers',
            type = "scatter3d",
            marker= list(color="red")) %>% 
  layout(scene = list(xaxis=list(title="gamma"),
                      yaxis=list(title="beta"), 
                      zaxis=list(title="log likelihood"))) 

landscape3D

# Sigma vs. Gamma

#2D plot
ggplot(profile_LL[profile_LL$beta==0.55,]) +
  geom_contour_filled(aes(x=sigma, y=gamma, z=ll), bins=50) +
  theme(legend.position = "none") +
  geom_point(aes(x=theta[["sigma"]], y=theta[["gamma"]]), colour="red")


# 3D plot
z <- profile_LL %>% filter(beta==0.55) %>% 
  select(sigma, gamma, ll) %>% 
  pivot_wider(id_cols = sigma, names_from=gamma, values_from = ll) %>% 
  select(-sigma) %>% as.matrix()

landscape3D <- plot_ly(z = ~ z, 
                       x = ~ range, 
                       y = ~ range,
                       type="surface") %>% 
  hide_colorbar() %>% 
  add_trace(x=gamma,
            y=sigma,
            z=max(z),
            mode='markers',
            type = "scatter3d",
            marker= list(color="red")) %>% 
  layout(scene = list(xaxis=list(title="gamma"),
                      yaxis=list(title="sigma"), 
                      zaxis=list(title="log likelihood"))) 

landscape3D



