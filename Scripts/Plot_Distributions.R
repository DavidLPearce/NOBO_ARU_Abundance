# --------------------------------------------------------------
# 
#                         Prior Check
#
# --------------------------------------------------------------

library(tidyverse)
library(ggfortify)

# -----------------------------
# Normal
# -----------------------------

# Parameters

mean <- 0

# sigma <- 0.75
# tau <- 1 / sigma^2

tau <- 0.2
sigma <- 1 / sqrt(tau) # exp(sigma)

# precision <- 2
# standdev <- 1 / sqrt(precision) 

# standdev = 3
# precision <- 1 / standdev^2
# precision
# variance <- standdev^2 # 1/precision
# variance

 
# # 95% of values fall in this range
# qnorm(0.025, mean = mean, sd = sigma)   
# qnorm(0.975, mean = mean, sd = sigma)  
# 
# # Which corresponds to delta values:
# exp(qnorm(0.025, mean = mean, sd = sigma))  
# exp(qnorm(0.975, mean = mean, sd = sigma))    

# Plot
ggdistribution(
  dnorm,
  x = seq(-10 * sigma, 10 * sigma, length.out = 10000),
  mean = mean,
  sd = sigma,
  xlab = "Value",
  ylab = "Density",
  fill = "black"
) +
  ggtitle(paste0("Normal Distribution\n", 
                 "mean = ",  round(mean, 3), "\n",
                 "tau = ", tau, "\n",
                 "sigma = ", round(sigma, 3)
                 )) +
  theme_classic()



# -----------------------------
# Log-Normal
# -----------------------------

# Parameters
mean = 6
logmean = log(mean)
# sigma = 1.5
# tau = 1/sigma^2

tau <- 0.2
sigma <- 1 / sqrt(tau) 

# Plot
ggdistribution(
  dlnorm,
  x = seq(0, 20, length.out = 1000),
  meanlog = logmean,
  sdlog = sigma,
  xlab = "Value",
  ylab = "Density",
  fill = "black"
) +
  ggtitle(paste0("Log-Normal Distribution\n",
                 "mean = ",  round(mean, 3), "\n",
                 "log scale mean = ",  round(logmean, 3), "\n",
                 "tau = ", tau, "\n",
                 "sigma = ", round(sigma, 3)
  )) +
  theme_classic()

# -----------------------------
# zero-truncated Normal
# -----------------------------

# ztNorm function
dtruncnorm <- function(x, mean, sd) {
  ifelse(x > 0,
         dnorm(x, mean, sd) /
           (1 - pnorm(0, mean, sd)),
         0)
}

# Parameters

mean <- log(6)

sigma <- 0.5
tau <- 1 / sigma^2

# tau <- 4
# sigma <- 1 / sqrt(tau)

# Plot
ggdistribution(
  dtruncnorm,
  x = seq(0, 5 * sigma + mean, length.out = 2000),
  mean = mean,
  sd = sigma,
  xlab = "Value",
  ylab = "Density",
  fill = "black"
) +
  ggtitle(paste0(
    "Zero-Truncated Normal Distribution\n",
    "mean = ", mean, "\n",
    "sigma = ", round(sigma, 3), "\n",
    "tau = ", round(tau, 3)
  )) +
  theme_classic()



# -----------------------------
# Gamma
# -----------------------------

# Parameters
shape <- 1
rate  <- 0.1

# Plot
ggdistribution(
  dgamma,
  x = seq(0, qgamma(0.999, shape, rate), length.out = 1000),
  shape = shape,
  rate = rate,
  xlab = "Value",
  ylab = "Density",
  fill = "black"
) +
  ggtitle(paste0("Gamma Distribution\n",
                 "shape = ", shape, "\n",
                 "rate = ", rate
                 )) +
  theme_classic()
