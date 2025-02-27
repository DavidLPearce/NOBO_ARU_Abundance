
# install.packages('plotly')
# install.packages('reshape2')

# Load necessary packages
library(ggplot2)
library(plotly)
library(reshape2)

# -----------------------------
# Woody proportion
# -----------------------------

# Set number of samples
n_samples <- 10000  # Number of draws from priors
n_points <- 1000    # Number of WoodyPrp values per sample

# Define priors (U-shaped effect)
beta1_samples <- rnorm(n_samples, -1, 1)  # Slightly negative to slope downward initially
beta2_samples <- rnorm(n_samples, 1, 0.5)  # Positive for U-shape

# Generate WoodyPrp values
WoodyPrp <- seq(0, 1, length.out = n_points)

# Compute the effect
effect_samples <- matrix(NA, nrow = n_samples, ncol = length(WoodyPrp))

for (i in 1:n_samples) {
  effect_samples[i, ] <- beta1_samples[i] * WoodyPrp + beta2_samples[i] * (WoodyPrp - 0.6)^2
}

# Convert to a dataframe
effect_df <- data.frame(
  WoodyPrp = rep(WoodyPrp, n_samples),
  Effect = as.vector(effect_samples),
  Sample = rep(1:n_samples, each = length(WoodyPrp))
)

# Plot
ggplot(effect_df, aes(x = WoodyPrp, y = Effect, group = Sample)) +
  geom_line(alpha = 0.1, color = "blue") +
  labs(title = "Prior Predictive Check: U-Shaped Effect of WoodyPrp",
       subtitle = "Effect is high at low/high WoodyPrp, lowest at ~0.6",
       x = "WoodyProportion (WoodyPrp)",
       y = "Predicted Effect") +
  theme_minimal()



# -----------------------------
# Mean Patch area : Number of patches
# -----------------------------


set.seed(42)  # For reproducibility

# Generate Mean Patch Area (MPA) values between 0 and 1
MPA <- runif(500, min = 0, max = 1)  

# Define inverse exponential relationship parameters
a <- 1000  # Maximum NP when MPA is near 0
b <- 5     # Controls steepness of relationship

# Simulate Number of Patches (NP) with some noise
NP <- a * exp(-b * MPA) + rnorm(500, mean = 0, sd = 50)

# Ensure NP is within realistic bounds (no negative values)
NP <- pmax(NP, 1)

# Create a dataframe
simulated_data <- data.frame(MPA = MPA, NP = NP)

# plot
ggplot(simulated_data, aes(x = MPA, y = NP)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", col = "forestgreen") +
  labs(title = "Simulated Relationship: Mean Patch Area vs. Number of Patches",
       x = "Mean Patch Area (MPA)",
       y = "Number of Patches (NP)") +
  theme_minimal()


# Standardize predictors (MPA & NP) and ensure they are numeric
simulated_data$MPA_scaled <- as.numeric(scale(simulated_data$MPA))
simulated_data$NP_scaled <- as.numeric(scale(simulated_data$NP))

# Fit an interaction model where NP is log-transformed
interaction_model <- lm(log(NP) ~ MPA_scaled * NP_scaled, data = simulated_data)

# Create prediction grid
MPA_seq <- seq(min(simulated_data$MPA_scaled), max(simulated_data$MPA_scaled), length.out = 50)
NP_seq <- seq(min(simulated_data$NP_scaled), max(simulated_data$NP_scaled), length.out = 50)
prediction_grid <- expand.grid(MPA_scaled = MPA_seq, NP_scaled = NP_seq)

# Predict interaction effect
prediction_grid$predicted_effect <- predict(interaction_model, newdata = prediction_grid)

# Reshape for plotting
prediction_matrix <- acast(prediction_grid, MPA_scaled ~ NP_scaled, value.var = "predicted_effect")

# 3D Surface Plot
plot_ly(x = unique(prediction_grid$MPA_scaled), 
        y = unique(prediction_grid$NP_scaled), 
        z = prediction_matrix, 
        type = "surface",
        colors = c("red", "yellow", "green")) %>%
  layout(title = "3D Surface Plot: MPA & NP Interaction",
         scene = list(xaxis = list(title = "Mean Patch Area (Scaled)"),
                      yaxis = list(title = "Number of Patches (Scaled)"),
                      zaxis = list(title = "Predicted Effect"),
                      camera = list(eye = list(x = -1.5, y = 1.5, z = 1.2))))  # Adjust viewing angle

# Heatmap Representation
ggplot(prediction_grid, aes(x = MPA_scaled, y = NP_scaled, fill = predicted_effect)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +  # Color scheme for better visualization
  labs(title = "Heatmap of MPA & NP Interaction Effect",
       x = "Mean Patch Area (Scaled)",
       y = "Number of Patches (Scaled)",
       fill = "Predicted Effect") +
  theme_minimal()


# Contour Plot
ggplot(prediction_grid, aes(x = MPA_scaled, y = NP_scaled, z = predicted_effect)) +
  geom_contour_filled(bins = 15) +  # Filled contour with 15 levels
  scale_fill_viridis_d(option = "magma") +  # Use "d" for discrete values
  labs(title = "Contour Plot: MPA & NP Interaction Effect",
       x = "Mean Patch Area (Scaled)",
       y = "Number of Patches (Scaled)",
       fill = "Predicted Effect") +
  theme_minimal()


