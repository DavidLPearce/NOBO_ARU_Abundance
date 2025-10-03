
# Load library
library(tidyverse)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")


# ---------------------------------------
# Conspecific Density Simulation 
# ---------------------------------------

# Constrain kappa1 to be between  0.6 +/- 0.4
kappa1 <- rnorm(50, mean = 0.6, sd = 0.4)   #rnorm(50, mean = 2, sd = 1)   # 

# Constrain kappa2 to be between -0.06 +/- 0.04
kappa2 <- rnorm(50, mean = -0.06, sd = 0.02)  

# Intercept: log(6) so that exp(gamma0) = 6 calls per 30 min for one individual, 2 calls per 10 mins
gamma0 <- log(2)  
exp(gamma0)

# Variance of random effect: precision
tau <- rgamma(1, shape = 2, rate = 2) # rgamma(1, shape = 50, rate = 25) 

# Convert tau to standard deviation (precision interpretation)
sigma <- sqrt(1 / tau)

# Number of surveys
num_surveys <- 14 

# Generate survey random effects
# raneff <- rnorm(num_surveys, mean = 0, sd = sigma)
#raneff <- pmax(pmin(rnorm(num_surveys, mean = 0, sd = sigma), 1), -1)  
raneff <- rt(num_surveys, df = 10) * sigma  # df = 10 gives a moderate tail


# Define a sequence of N values to evaluate the function
N <- 0:10

# Create an empty data frame to store the results
results <- data.frame()

# Simulation
for (i in 1:length(kappa1)) {
  
  # Assign a survey  
  survey <- sample(1:num_surveys, 1, replace = TRUE)
  
  # Call rate model
  #delta <-  exp(gamma0+ raneff[survey] ) + ((kappa1[i]*N) + (kappa2[i]*(N)^2)) # 
  delta <- pmax(0, exp(gamma0 + raneff[survey]) + (kappa1[i] * N) + (kappa2[i] * (N)^2))
  
  # Temp holding df
  tmp <- data.frame(N = N, delta = delta, kappa1 = kappa1[i], kappa2 = kappa2[i], 
                    combination = paste0("k1=", round(kappa1[i], 2), ", k2=", round(kappa2[i], 2)))
  
  # Combine Results
  results <- rbind(results, tmp)
}

# Take a look
head(results, 50)


# Color scheme
colors <- scico::scico(100, palette = "berlin")

# Plot 
ggplot(results, aes(x = N, y = delta, color = combination)) +
  geom_line(linewidth = 1) +  
  labs(x = "Calling Males (N)", y = "Call Rate (δ)",
       title = "Conspecifics on Call Rate Simulation") +
  scale_color_manual(values = colors) +   
  theme_minimal() +   
  theme(legend.title = element_blank())   

