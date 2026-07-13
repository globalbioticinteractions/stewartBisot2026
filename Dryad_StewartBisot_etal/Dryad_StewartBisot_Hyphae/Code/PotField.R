# =============================================================================
# Description:
# This script analyzes observed hyphal density to test for differences between
# pot based and field based studies, and to evaluate effects of artificial
# versus natural light within pot studies. The workflow loads and cleans the
# dataset, filters observations with complete metadata, and fits mixed effects
# models with study DOI as a random intercept. Separate models compare pot
# versus field conditions and artificial versus natural light conditions.
# The script produces marginal effect plots and raw data visualizations for
# both comparisons and exports the resulting figures.
# =============================================================================


# === Load libraries ===
library(lme4)
library(lmerTest)
library(dplyr)
library(sjPlot)
library(tidyverse)
library(ggplot2)
library(patchwork)

# === Load and clean data ===
data <- read.csv("hyphal_density_full_metadata")

data$Growth_Form_Final[data$Growth_Form_Final == ""] <- NA
data$pot_field <- tolower(trimws(data$pot_field))
data$pot_field[data$pot_field == ""] <- NA
data$pot_artificial_natural_light <- tolower(trimws(data$pot_artificial_natural_light))
data$pot_artificial_natural_light[data$pot_artificial_natural_light == ""] <- NA

# === === 1. POT vs FIELD COMPARISON === ===

# Filter to complete entries
fielddata <- data %>%
  filter(!is.na(HyphalLength_m_cm3_final), !is.na(pot_field), !is.na(DOI))

# Fit null and full models
null_pf <- lmer(HyphalLength_m_cm3_final ~ 1 + (1 | DOI), data = fielddata, REML = FALSE)
mod_pf <- lmer(HyphalLength_m_cm3_final ~ pot_field + (1 | DOI), data = fielddata, REML = FALSE)

# Marginal effect plot
p1 <- plot_model(mod_pf, type = "eff", terms = "pot_field", show.data = TRUE) +
  ggtitle("Marginal Effect: Pot vs Field")

# Raw data plot
p2 <- ggplot(fielddata, aes(x = pot_field, y = log(HyphalLength_m_cm3_final), fill = pot_field)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  scale_fill_manual(values = c("field" = "forestgreen", "pot" = "steelblue")) +
  theme_minimal(base_size = 14) +
  labs(title = "Observed Hyphal Density by Study Type", x = "Study Type", y = "log(Hyphal Density m/cm³)") +
  theme(legend.position = "none")

# === === 2. ARTIFICIAL vs NATURAL LIGHT (POT ONLY) === ===

# Filter pot-only
potdata <- data %>%
  filter(
    pot_field == "pot",
    !is.na(HyphalLength_m_cm3_final),
    !is.na(pot_artificial_natural_light),
    !is.na(DOI)
  )

# Fit models
null_light <- lmer(HyphalLength_m_cm3_final ~ 1 + (1 | DOI), data = potdata, REML = FALSE)
mod_light <- lmer(HyphalLength_m_cm3_final ~ pot_artificial_natural_light + (1 | DOI), data = potdata, REML = FALSE)
anova(null_light,mod_light)


# Marginal effect plot
p3 <- plot_model(mod_light, type = "eff", terms = "pot_artificial_natural_light", show.data = TRUE) +
  ggtitle("Marginal Effect: Artificial vs Natural Light")

# Raw data plot
p4 <- ggplot(potdata, aes(x = pot_artificial_natural_light, y = log(HyphalLength_m_cm3_final), fill = pot_artificial_natural_light)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  scale_fill_manual(values = c("natural" = "goldenrod", "artificial" = "dodgerblue")) +
  theme_minimal(base_size = 14) +
  labs(title = "Observed Hyphal Density by Light Type", x = "Light Source", y = "log(Hyphal Density m/cm³)") +
  theme(legend.position = "none")


# === Final plot: Pot vs Field ===
plot_pf <- ggplot(fielddata, aes(x = pot_field, y = log(HyphalLength_m_cm3_final), fill = pot_field)) +
  geom_jitter(aes(color = pot_field), width = 0.2, size = 1.2, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("field" = "forestgreen", "pot" = "steelblue")) +
  scale_color_manual(values = c("field" = "forestgreen", "pot" = "steelblue")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hyphal Density by Study Type", x = "Study Type", y = "log(Hyphal Density m/cm³)") +
  theme(legend.position = "none")

ggsave(
  "pot_field.pdf",device = "pdf",
  plot = plot_pf, width = 6, height = 5, dpi = 300
)

# === Final plot: Artificial vs Natural Light ===
plot_light <- ggplot(potdata, aes(x = pot_artificial_natural_light, y = log(HyphalLength_m_cm3_final), fill = pot_artificial_natural_light)) +
  geom_jitter(aes(color = pot_artificial_natural_light), width = 0.2, size = 1.2, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("natural" = "goldenrod", "artificial" = "dodgerblue")) +
  scale_color_manual(values = c("natural" = "goldenrod", "artificial" = "dodgerblue")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hyphal Density by Light Type (Pot Studies)", x = "Light Source", y = "log(Hyphal Density m/cm³)") +
  theme(legend.position = "none")


ggsave(
  "/pot_light.pdf",device = "pdf",
  plot = plot_light, width = 6, height = 5, dpi = 300
)
