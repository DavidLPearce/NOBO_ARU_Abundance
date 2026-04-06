


library(tidyverse)


# Load in capture data
pc_dat <- read.csv("./Data/Point_Count_Data/NOBO_PC_Summer2024data.csv")

# Load in site covariates
site_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")
colnames(site_covs)[1] <- "PointNum"# Rename SiteID to PointNum for matching



 
# Subset to surveys 2 and 3  
 
pc_CMR <- pc_dat %>%
  filter(Survey %in% c(2, 3))

 
# Check the first few rows
head(pc_CMR, 10)


# Function to get first detection interval
get_first_detection <- function(int1, int2, int3, int4) {
  # Handle NA rows (no detections during survey)
  if(all(is.na(c(int1, int2, int3, int4)))) {
    return(NA)  
  }
  
  # Find first detection
  if(int1 == 1) return(1)
  if(int2 == 1) return(2)
  if(int3 == 1) return(3)
  if(int4 == 1) return(4)
  return(5)  
}

# Add first detection variable
pc_CMR <- pc_CMR %>%
  rowwise() %>%
  mutate(first_detection = get_first_detection(int1, int2, int3, int4)) %>%
  ungroup()

# Remove rows with NA 
pc_CMR <- pc_CMR %>%
  filter(!is.na(first_detection))

 
# Distribution of first detections
print(table(pc_CMR$first_detection))

# Add a Unique identifier  
pc_CMR <- pc_CMR %>%
     mutate(UniqueID = paste(PointNum, Survey, NOBOnum, sep = "_"))

# Recode Survey to Visit Survey 2 > Visit 1, Survey 3 > Visit 2
pc_CMR <- pc_CMR %>%
  mutate(Visit = ifelse(Survey == 2, 1, 2))

# Counts by visit
print(table(pc_CMR$Survey, pc_CMR$Visit))
 
# Inspect
pc_CMR %>%
  select(PointNum, 
         Visit, 
         NOBOnum, 
         int1, int2, int3, int4, 
         first_detection) %>%
  head(15) %>%
  print()


# Define constants
S <- 10  # Number of sites
K <- 2   # Number of survey occasions 
J <- 4   # Number of time intervals


# Create ordered site list
site_list <- sort(unique(pc_CMR$PointNum))
 

# Check visits at each site 
print(table(pc_CMR$PointNum, pc_CMR$Visit))

# missing a few
 
# All 10 sites were surveyed twice
all_sites <- 1:S
all_visits <- 1:K

complete_surveys <- expand.grid(
  PointNum = all_sites,
  Visit = all_visits
) %>%
  arrange(PointNum, Visit)
 
print(complete_surveys)
 
# ----------------------
# Covariates 
# ----------------------

# Convert Wind to numeric
pc_CMR <- pc_CMR %>%
  mutate(Wind = as.numeric(Wind.Beau.Code))

# Number of levels
Wind_Lvls <- length(unique(pc_CMR$Wind))

# Get wind cov for every visit 
survey_covariates <- pc_dat %>%
  filter(Survey %in% c(2, 3)) %>%
  mutate(Visit = ifelse(Survey == 2, 1, 2)) %>%
  group_by(PointNum, Visit) %>%
  summarise(
    Wind = first(Wind.Beau.Code),
    Survey = first(Survey),
    Date = first(Date),
    .groups = "drop"
  ) %>%
  arrange(PointNum, Visit)
 
print(survey_covariates)
 
# Create Wind matrix: site x visit 
Wind_matrix <- matrix(NA, nrow = S, ncol = K)

for(i in 1:nrow(survey_covariates)) {
  site_num <- survey_covariates$PointNum[i]
  visit <- survey_covariates$Visit[i]
  Wind_matrix[site_num, visit] <- survey_covariates$Wind[i]
}
print(Wind_matrix)

# Detection covariates  
VegDens <- as.vector(scale(site_covs$vegDens50m))

# Abundance covariates 
Herb_COH <- as.vector(scale(site_covs$herb_COH))
Woody_SPLIT <- as.vector(scale(site_covs$woody_SPLIT))
 







# ------------------------------------------
# Number of individuals detected
# ------------------------------------------
nind <- nrow(pc_CMR)
cat("=== Individual Detection Data ===\n")
cat("Number of individuals detected:", nind, "\n\n")

# ------------------------------------------
# Create detection variable (y)
# ------------------------------------------
y <- pc_CMR$first_detection

cat("Distribution of first detections:\n")
print(table(y))
cat("\n")

# ------------------------------------------
# Create site assignment for each individual
# ------------------------------------------
# Site number (1-10) for each detected individual
site <- pc_CMR$PointNum

cat("Individuals per site:\n")
print(table(site))
cat("\n")

# ------------------------------------------
# Create visit assignment for each individual
# ------------------------------------------
visit <- pc_CMR$Visit

cat("Individuals per visit:\n")
print(table(visit))
cat("\n")

# ------------------------------------------
# Verify data alignment
# ------------------------------------------
cat("Sample of individual data (first 15):\n")
individual_sample <- data.frame(
  Individual = 1:min(15, nind),
  UniqueID = pc_CMR$UniqueID[1:min(15, nind)],
  PointNum = pc_CMR$PointNum[1:min(15, nind)],
  site = site[1:min(15, nind)],
  Survey = pc_CMR$Survey[1:min(15, nind)],
  visit = visit[1:min(15, nind)],
  first_det = y[1:min(15, nind)]
)
print(individual_sample)

# ------------------------------------------
# Check for any issues
# ------------------------------------------
cat("\nData validation:\n")
cat("Range of site:", range(site), "should be [1, 10]\n")
cat("Range of visit:", range(visit), "should be [1, 2]\n")
cat("Range of first_detection:", range(y), "should be [1, 4]\n")







# ------------------------------------------
# Data Augmentation
# ------------------------------------------
# M = augmented population size (detected + potential undetected)
M <- nind + 100  # Add 100 potential individuals

cat("=== Data Augmentation ===\n")
cat("Detected individuals (nind):", nind, "\n")
cat("Augmented population (M):", M, "\n")
cat("Augmented individuals:", M - nind, "\n\n")

# ------------------------------------------
# Augment y (detection data)
# ------------------------------------------
# Detected individuals: 1-4 (their first detection interval)
# Augmented individuals: 5 (never detected)
y_aug <- c(y, rep(5, M - nind))

cat("Augmented y length:", length(y_aug), "\n")
cat("Detected (1-4):", sum(y_aug <= 4), "\n")
cat("Never detected (5):", sum(y_aug == 5), "\n\n")

# ------------------------------------------
# Augment site assignment
# ------------------------------------------
# Detected individuals: known sites (2-10)
# Augmented individuals: NA (unknown, to be estimated)
site_aug <- c(site, rep(NA, M - nind))

cat("Augmented site length:", length(site_aug), "\n")
cat("Known sites:", sum(!is.na(site_aug)), "\n")
cat("Unknown sites (NA):", sum(is.na(site_aug)), "\n\n")

# ------------------------------------------
# Augment visit assignment
# ------------------------------------------
# Detected individuals: known visits (1 or 2)
# Augmented individuals: NA (unknown, to be estimated)
visit_aug <- c(visit, rep(NA, M - nind))

cat("Augmented visit length:", length(visit_aug), "\n")
cat("Known visits:", sum(!is.na(visit_aug)), "\n")
cat("Unknown visits (NA):", sum(is.na(visit_aug)), "\n\n")

# ------------------------------------------
# Summary
# ------------------------------------------
cat("=== Augmented Data Summary ===\n")
cat("Total individuals in model (M):", M, "\n")
cat("  - Detected:", nind, "\n")
cat("  - Augmented (potential):", M - nind, "\n")
cat("\nAugmented data structures:\n")
cat("  - y_aug: length", length(y_aug), "\n")
cat("  - site_aug: length", length(site_aug), "\n")
cat("  - visit_aug: length", length(visit_aug), "\n")







# ------------------------------------------
# Set up prediction data (optional - can skip for now)
# ------------------------------------------
U <- 0  # No unsampled sites for prediction
Herb_COH_pred <- numeric(0)
Woody_SPLIT_pred <- numeric(0)
rho <- numeric(0)

# ------------------------------------------
# Bundle all data for JAGS
# ------------------------------------------
data <- list(
  # Sample size
  M = M,
  S = S,
  K = K,
  nind = nind,
  
  # Individual-level data
  y = y_aug,
  site = site_aug,
  visit = visit_aug,
  
  # Survey-level covariates
  Wind = Wind_matrix,        # [S × K] matrix
  Wind_Lvls = Wind_Lvls,
  
  # Site-level covariates
  VegDens = VegDens,         # [S] vector - for detection
  Herb_COH = Herb_COH,       # [S] vector - for abundance
  Woody_SPLIT = Woody_SPLIT, # [S] vector - for abundance
  
  # Prediction (empty for now)
  U = U,
  Herb_COH_pred = Herb_COH_pred,
  Woody_SPLIT_pred = Woody_SPLIT_pred,
  rho = rho
)

# ------------------------------------------
# Check the bundled data structure
# ------------------------------------------
cat("=== JAGS Data Bundle ===\n\n")
str(data)

cat("\n=== Dimension Checks ===\n")
cat("M (augmented pop):", data$M, "\n")
cat("S (sites):", data$S, "\n")
cat("K (visits):", data$K, "\n")
cat("nind (detected):", data$nind, "\n")
cat("Wind_Lvls:", data$Wind_Lvls, "\n\n")

cat("Individual-level data:\n")
cat("  y: length", length(data$y), "should equal M:", data$M, "\n")
cat("  site: length", length(data$site), "should equal M:", data$M, "\n")
cat("  visit: length", length(data$visit), "should equal M:", data$M, "\n\n")

cat("Survey-level covariates:\n")
cat("  Wind: dimensions", dim(data$Wind), "should be [S, K]:", c(data$S, data$K), "\n\n")

cat("Site-level covariates:\n")
cat("  VegDens: length", length(data$VegDens), "should equal S:", data$S, "\n")
cat("  Herb_COH: length", length(data$Herb_COH), "should equal S:", data$S, "\n")
cat("  Woody_SPLIT: length", length(data$Woody_SPLIT), "should equal S:", data$S, "\n")

cat("\n data bundle complete and ready for JAGS!\n")




#################


cat("
model {

  # ---------------------------------
  # Priors: Abundance Model
  # ---------------------------------
  
  beta0 ~ dnorm(0, 0.1)          # Intercept
  beta1 ~ dnorm(0, 0.1)          # Herbaceous cohesion effect
  beta2 ~ dnorm(0, 0.1)          # Woody splitting effect
  
  # ---------------------------------
  # Priors: Detection Model
  # ---------------------------------
  
  # Baseline detection probability (per interval)
  p0 ~ dunif(0, 1)
  alpha0 <- log(p0 / (1 - p0))   # Logit scale
  
  # Wind effect (categorical - 3 levels: 0, 1, 2)
  for(w in 1:Wind_Lvls) {
    alpha1[w] ~ dnorm(0, 0.1)
  }
  
  # Vegetation density effect
  alpha2 ~ dnorm(0, 0.1)
  
  # Site random effect (accounts for pseudoreplication across visits)
  tau_site ~ dgamma(0.1, 0.1)
  sd_site <- 1 / sqrt(tau_site)
  
  for(s in 1:S) {
    site_re[s] ~ dnorm(0, tau_site)
  }
  
  # ---------------------------------
  # Abundance Model (Spatial Process)
  # ---------------------------------
  
  for(s in 1:S) {
    # Expected abundance at site s
    log(lambda[s]) <- beta0 + beta1 * Herb_COH[s] + beta2 * Woody_SPLIT[s]
    
    # Site probability (for site membership)
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance (sum of individuals assigned to site s)
    N[s] <- sum(z[] * equals(group[], s))
  }
  
  # Overall inclusion probability
  psi <- sum(lambda[]) / M
  
  # ---------------------------------
  # Individual-Level Model
  # ---------------------------------
  
  for(i in 1:M) {
    
    # Presence/absence (data augmentation variable)
    z[i] ~ dbern(psi)
    
    # Site membership (which site is individual i at?)
    group[i] ~ dcat(probs[])
    
    # Visit assignment (which visit was individual i detected in?)
    # Uniform across K visits
    visit_assignment[i] ~ dcat(visit_probs[])
    
    # Detection probability for this individual
    # Depends on: baseline + Wind at their site-visit + VegDens at their site + site random effect
    logit(p[i]) <- alpha0 + 
                   alpha1[Wind[group[i], visit_assignment[i]]] + 
                   alpha2 * VegDens[group[i]] + 
                   site_re[group[i]]
    
    # ---------------------------------
    # Time-of-Detection (Removal) Model
    # ---------------------------------
    
    # Multinomial cell probabilities
    # pi[i,1] = detected in interval 1
    # pi[i,2] = detected in interval 2 (not detected in 1)
    # pi[i,3] = detected in interval 3 (not detected in 1 or 2)
    # pi[i,4] = detected in interval 4 (not detected in 1, 2, or 3)
    # pi[i,5] = never detected
    
    pi[i,1] <- p[i] * z[i]
    pi[i,2] <- (1 - p[i]) * p[i] * z[i]
    pi[i,3] <- pow(1 - p[i], 2) * p[i] * z[i]
    pi[i,4] <- pow(1 - p[i], 3) * p[i] * z[i]
    pi[i,5] <- 1 - sum(pi[i, 1:4])
    
    # Likelihood: categorical distribution
    y[i] ~ dcat(pi[i, ])
  }
  
  # Uniform visit probabilities
  for(k in 1:K) {
    visit_probs[k] <- 1 / K
  }
  
  # ---------------------------------
  # Derived Quantities
  # ---------------------------------
  
  # Total abundance at sampled sites
  N_tot <- sum(z[])
  
  # Mean abundance per site
  N_mean <- N_tot / S
  
  # Expected lambda (mean and total)
  lambda_mean <- mean(lambda[])
  lambda_tot <- sum(lambda[])
  
}
", fill = TRUE, file = "./jags_models/Model_TOD_MCR.txt")

cat("✅ JAGS model written to: ./jags_models/Model_TOD_MCR.txt\n")





# ------------------------------------------
# Initial values function
# ------------------------------------------
inits <- function() {
  list(
    # Detection parameters
    p0 = runif(1, 0.3, 0.7),
    alpha1 = rnorm(Wind_Lvls, 0, 0.5),
    alpha2 = rnorm(1, 0, 0.5),
    
    # Abundance parameters
    beta0 = rnorm(1, 0, 0.5),
    beta1 = rnorm(1, 0, 0.5),
    beta2 = rnorm(1, 0, 0.5),
    
    # Data augmentation variable
    z = c(rep(1, nind), rbinom(M - nind, 1, 0.5)),
    
    # Site random effects
    site_re = rnorm(S, 0, 0.1)
  )
}

# Test the initial values function
cat("Testing initial values function...\n")
test_inits <- inits()
str(test_inits)
cat("✅ Initial values function working!\n\n")

# ------------------------------------------
# Parameters to monitor
# ------------------------------------------
params <- c(
  # Detection parameters
  "p0", "alpha0", "alpha1", "alpha2",
  
  # Abundance parameters
  "beta0", "beta1", "beta2",
  
  # Random effects
  "sd_site", "site_re",
  
  # Abundance estimates
  "N", "N_tot", "N_mean",
  
  # Lambda (expected abundance)
  "lambda", "lambda_mean", "lambda_tot",
  
  # Inclusion probability
  "psi"
)

cat("Parameters to monitor:\n")
print(params)
cat("\n✅ Parameters list ready!\n")











# ------------------------------------------
# Recode Wind matrix to be 1-indexed for JAGS
# ------------------------------------------
cat("Original Wind values:", sort(unique(as.vector(Wind_matrix))), "\n")

# Add 1 to all Wind values (0→1, 1→2, 2→3)
Wind_matrix_jags <- Wind_matrix + 1

cat("Recoded Wind values:", sort(unique(as.vector(Wind_matrix_jags))), "\n")

# Update Wind_Lvls
Wind_Lvls <- length(unique(as.vector(Wind_matrix_jags)))
cat("Wind_Lvls:", Wind_Lvls, "\n\n")

# ------------------------------------------
# Update the data bundle
# ------------------------------------------
data$Wind <- Wind_matrix_jags
data$Wind_Lvls <- Wind_Lvls

cat("Updated data bundle:\n")
cat("Wind matrix:\n")
print(data$Wind)
cat("\nWind_Lvls:", data$Wind_Lvls, "\n")

# ------------------------------------------
# Verify the fix
# ------------------------------------------
cat("\nVerification:\n")
cat("Min Wind value:", min(data$Wind), "should be >= 1\n")
cat("Max Wind value:", max(data$Wind), "should be <=", data$Wind_Lvls, "\n")





##############


# ------------------------------------------
# Load jagsUI package
# ------------------------------------------
library(jagsUI)

# ------------------------------------------
# MCMC settings
# ------------------------------------------
n_adapt <- 50000      # Adaptation phase
n_burnin <- 0     # Burn-in iterations
n_iter <- 30000      # Total iterations (including burn-in)
n_thin <- 5          # Thinning rate
n_chains <- 3        # Number of chains

cat("=== MCMC Settings ===\n")
cat("Adaptation:", n_adapt, "\n")
cat("Burn-in:", n_burnin, "\n")
cat("Total iterations:", n_iter, "\n")
cat("Thinning:", n_thin, "\n")
cat("Chains:", n_chains, "\n")
cat("Samples per chain:", (n_iter - n_burnin) / n_thin, "\n\n")

# ------------------------------------------
# Run the model!
# ------------------------------------------
cat("Starting JAGS model...\n")
cat("This may take a few minutes...\n\n")

start_time <- Sys.time()

fm1 <- jags(
  data = data,
  parameters.to.save = params,
  inits = inits,
  model.file = "./jags_models/Model_TOD_MCR.txt",
  n.iter = n_iter,
  n.burnin = n_burnin,
  n.chains = n_chains,
  n.thin = n_thin,
  n.adapt = n_adapt,
  parallel = TRUE,
  n.cores = 3,
  verbose = TRUE,
  DIC = FALSE
)

end_time <- Sys.time()
run_time <- end_time - start_time

cat("\n✅ Model completed!\n")
cat("Run time:", round(run_time, 2), attr(run_time, "units"), "\n\n")

# ------------------------------------------
# Print model summary
# ------------------------------------------
print(fm1)

# ------------------------------------------
# Check convergence (Rhat values)
# ------------------------------------------
cat("\n=== Convergence Check ===\n")
rhat_values <- fm1$summary[, "Rhat"]
cat("Range of Rhat values:", round(min(rhat_values, na.rm = TRUE), 3), 
    "to", round(max(rhat_values, na.rm = TRUE), 3), "\n")
cat("Parameters with Rhat > 1.1:", sum(rhat_values > 1.1, na.rm = TRUE), "\n")

if(sum(rhat_values > 1.1, na.rm = TRUE) == 0) {
  cat("✅ All parameters converged (Rhat < 1.1)!\n")
} else {
  cat("⚠️ Some parameters may not have converged. Consider longer chains.\n")
}



# Trace plots
MCMCvis::MCMCtrace(fm1, 
                   params = c(  
                     # Detection parameters
                     "p0", "alpha0", "alpha1", "alpha2",
                     
                     # Abundance parameters
                     "beta0", "beta1", "beta2",
                     
                     # Random effects
                     "sd_site", "site_re",
                     
                     # Abundance estimates
                     "N", "N_tot", "N_mean",
                     
                     # Lambda (expected abundance)
                     "lambda", "lambda_mean", "lambda_tot",
                     
                     # Inclusion probability
                     "psi"
                              
                   ),
                   pdf = T,
                   filename = "TracePlots_PC_TOD.pdf",
                   wd = "./Figures"
)
