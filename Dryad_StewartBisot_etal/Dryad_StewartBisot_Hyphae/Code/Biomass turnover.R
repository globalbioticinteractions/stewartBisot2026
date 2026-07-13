# =============================================================================
# Description:
# This script estimates seasonal turnover of AM fungal spore biomass using
# observed hyphal density data. The workflow loads and cleans the dataset,
# fits mixed effects models with season as a fixed effect and study identity
# as a random intercept, and compares seasonal and null models. Seasonal
# marginal means are scaled to a global biomass estimate to produce seasonal
# stock fluctuations. These fluctuations are used to calculate
# and a global turnover rate, with uncertainty propagated through Monte Carlo sampling.
# A second calculation estimates turnover based on an assumed global carbon
# input to hyphae from Hawkins 2023.
# =============================================================================


# --- Load libraries ---
library(dplyr)      
library(lme4)       
library(emmeans)    

# ================================================================
# 1. Load and clean data
# ================================================================
spore <- read.csv("hyphal_density_full_metadata") 

spore_clean <- spore %>%
  mutate(
    season = factor(season,
                    levels = c("Spring", "Summer", "Autumn", "Winter"), # ensure Fall included
                    ordered = TRUE),
    hyphal_density_m_cm3 = as.numeric(HyphalLength_m_cm3_final)
  ) %>%
  filter(!is.na(season), !is.na(hyphal_density_m_cm3))

# ================================================================
# 2. Fit mixed-effects models (frequentist)
# ================================================================
# Null model
null_mod <- lmer(
  log(hyphal_density_m_cm3 + 1) ~ 1 + (1 | DOI),
  data = spore_clean
)

# Seasonal model
season_mod <- lmer(
  (hyphal_density_m_cm3) ~ season + (1 | DOI),
  data = spore_clean
)
unique(spore_clean$season)
# ================================================================
# 3. Model comparison (AIC)
# ================================================================
print(AIC(null_mod, season_mod))

# ================================================================
# 4. Seasonal marginal means scaled to global stock
# ================================================================
B_total <- 297  # Mt
B_sd    <- 42   # Mt

emm <- emmeans(season_mod, ~ season, type = "response")
emm_df <- as.data.frame(emm)

# Relative seasonal fractions
rel <- with(emm_df, response / mean(response))

# Scale to global biomass
B_season <- B_total * rel

# ================================================================
# 5. Compute residence time τ
# ================================================================
peak_idx <- which.max(B_season)
peak_B   <- B_season[peak_idx]
trough_B <- B_season[(peak_idx %% length(B_season)) + 1]

# Quarter year step
delta_t <- 0.25

F_out <- (peak_B - trough_B) / delta_t
tau   <- B_total / F_out

# ================================================================
# 6. Uncertainty propagation (Monte Carlo)
# ================================================================
set.seed(1)
n <- 40000
tau_draw <- rnorm(n, B_total, B_sd) / F_out
tau_CI   <- quantile(tau_draw, c(0.025, 0.5, 0.975))

cat("Residence time τ (years):\n")
print(tau_CI)

# ================================================================
# 7. Compute turnover rate k (per year)
# ================================================================
k_draw   <- 1 / tau_draw
k_median <- median(k_draw)
k_sd     <- sd(k_draw)  # standard deviation as error term
k_CI     <- quantile(k_draw, c(0.025, 0.975))

cat("\nTurnover rate k (per year):\n")
cat(sprintf("Median: %.3f ± %.3f yr⁻¹\n", k_median, k_sd))
cat(sprintf("95%% CI: [%.3f, %.3f] yr⁻¹\n", k_CI[1], k_CI[2]))

# ================================================================
# Turnover calculation assuming 1000 Mt C annual input
# ================================================================

# Known values
B_total <- 300   # standing stock (Mt C)
B_sd    <- 50    # uncertainty (Mt C)
Input   <- 1000  # annual plant C input to hyphae (Mt C)

# Turnover rate formula: T = Input / B_total ± (Input * B_sd) / B_total^2
T_mean <- Input / B_total
T_sd   <- (Input * B_sd) / (B_total^2)

# Report result
cat(sprintf("Turnover rate T = %.2f ± %.2f times per year\n", T_mean, T_sd))

