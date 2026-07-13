# =============================================================================
# Description:
# This script performs a Monte Carlo simulation to estimate uncertainty in
# global carbon biomass from hyphal networks. It uses a global mean and
# standard deviation of hyphal length (in km), simulates variation in hyphal
# radius, and calculates carbon mass (in megatons, MT). The results include
# a summary data frame of mean and SD estimates and a plot visualizing how
# hyphal radius variation affects global carbon biomass.
# =============================================================================

##### Load necessary libraries #####
library(dplyr)
library(purrr)
library(ggplot2)

set.seed(123) # For reproducibility

##### Define constants for Monte Carlo simulation #####
pi <- 3.1415
dry_mass_ratio <- 0.21
cell_density <- 1.1
carbon_ratio <- 0.5
mean_radius <- 2.75e-9    # in km (since hyphal length is in km)
sd_radius <- 2.1e-10     # in km
simulation_size <- 1e6   # Number of simulations

# Global hyphal length parameters (in km)
global_mean_length <- 1.1e17
global_sd_length <- global_mean_length*.118
RMSE <- 1.5e9

##### Monte Carlo Simulation Function with RMSE propagation #####
perform_monte_carlo_global <- function(mean_length, sd_length, RMSE, n) {
  
  # Random hyphal lengths (km), including RMSE error term
  lengths <- rnorm(n, mean_length, sd_length) + rnorm(n, 0, RMSE)
  
  # Random hyphal radii (km)
  radii <- rnorm(n, mean_radius, sd_radius)
  
  # Calculate carbon biomass
  volume_hyp <- lengths * (pi * radii^2)   # Volume in km^3
  
  carbon <- volume_hyp * (dry_mass_ratio * carbon_ratio * cell_density)
  carbon_mt <- carbon * 1000  # convert to Mt
  
  tibble(
    Simulation_Size = n,
    Mean_Carbon_MT = mean(carbon_mt, na.rm = TRUE),
    SD_Carbon_MT = sd(carbon_mt, na.rm = TRUE)
  )
}

##### Run the updated simulation #####
results_df <- perform_monte_carlo_global(
  mean_length = global_mean_length,
  sd_length = global_sd_length,
  RMSE = RMSE,
  n = simulation_size
)

print(results_df)


##### Histogram simulation with RMSE #####
carbon_samples <- {
  lengths <- rnorm(simulation_size, global_mean_length, global_sd_length) +
    rnorm(simulation_size, 0, RMSE)
  
  radii <- rnorm(simulation_size, mean_radius, sd_radius)
  
  volume_hyp <- lengths * (pi * radii^2)
  
  carbon <- volume_hyp * (dry_mass_ratio * carbon_ratio * cell_density)
  carbon * 1000
}

ggplot(data.frame(Carbon_MT = carbon_samples), aes(x = Carbon_MT)) +
  geom_histogram(bins = 60, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Monte Carlo Simulation of Global Hyphal Carbon Biomass with RMSE Propagation",
    x = "Carbon Biomass (MT)",
    y = "Frequency"
  )
