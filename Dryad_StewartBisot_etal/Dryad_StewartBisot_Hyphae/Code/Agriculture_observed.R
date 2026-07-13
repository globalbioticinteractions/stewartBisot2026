# =============================================================================
# Description:
# This script analyzes observed hyphal density measurements to evaluate how
# values differ between cropland and non cropland soils across global biomes.
# The workflow loads the observed dataset, assigns each case to its biome,
# filters to biomes with sufficient representation of both land use types, and
# fits a mixed effects model with crop status, biome, and their interaction as
# fixed effects and study as a random intercept. Estimated marginal means are
# computed globally and per biome, backtransformed, and visualized with
# standard error intervals. The script also produces a figure showing sample
# sizes for each biome and crop category.
# =============================================================================

library(tidyverse)
library(sf)
library(lme4)
library(lmerTest)
library(emmeans)

# === Load data ===
df <- read.csv("hyphal_density_full_metadata")
ecoregions <- st_read("Ecoregions2017.shp") #https://ecoregions.appspot.com/ 
biome_lookup <- ecoregions %>%
  st_drop_geometry() %>%
  distinct(BIOME_NUM, BIOME_NAME) %>%
  mutate(RESOLVE = paste0("Biome ", BIOME_NUM))

df <- df %>%
  mutate(RESOLVE = paste0("Biome ", biomes)) %>%
  left_join(biome_lookup, by = "RESOLVE") %>%
  mutate(
    crop = factor(crops_ESA, levels = c(0, 1), labels = c("Non-Crop", "Crop")),
    hyphae = HyphalLength_m_cm3_final,
    study = as.factor(DOI),
    biome = as.factor(BIOME_NAME)
  ) %>%
  filter(biome != "Rock and Ice") %>%
  drop_na(hyphae, crop, phosphorus, CHELSA_BIO_Annual_Precipitation,
          nitrogen, isric.soil.proportion.of.sand, MODIS_NPP, study, biome)

# === Filter to biomes with ≥5 samples for both crop types ===
valid_biomes <- df %>%
  count(biome, crop) %>%
  pivot_wider(names_from = crop, values_from = n, values_fill = 0) %>%
  filter(Crop >= 5, `Non-Crop` >= 5) %>%
  pull(biome)
df <- df %>% filter(biome %in% valid_biomes)

# === Fit models on log-transformed hyphae ===
df <- df %>% mutate(hyphae = (hyphae)) # add small offset if needed

model <- lmer(hyphae ~ crop * biome + phosphorus +
                CHELSA_BIO_Annual_Precipitation + nitrogen +
                isric.soil.proportion.of.sand + EarthEnvTopoMed_AspectCosine +
                (1| study), data = df)
model_null <- lmer(hyphae ~ 1 +
                     (1 | study), data = df)
anova(model_null, model)

# === emmeans (log scale, backtransform) ===
emm_global <- emmeans(model, ~ crop, type = "response")
emm_biome <- emmeans(model, ~ crop | biome, type = "response")

# === Prepare for plot (with SEs) ===
global_df <- as.data.frame(emm_global) %>%
  rename(response = emmean, SE = SE) %>%
  mutate(biome = "Global",
         ymin = response - SE,
         ymax = response + SE)

biome_df <- as.data.frame(emm_biome) %>%
  rename(response = emmean, SE = SE) %>%
  mutate(ymin = response - SE,
         ymax = response + SE)

combined_df <- bind_rows(global_df, biome_df) %>%
  mutate(biome = factor(biome, levels = c("Global", sort(unique(biome[biome != "Global"])))))

# === Plot with SE error bars ===
p_density <- ggplot(combined_df, aes(x = biome, y = (response), fill = crop)) +
  geom_col(position = position_dodge(0.6), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(0.6), width = 0.2) +
  scale_fill_manual(values = c("#CCCCCC", "#2166AC")) +
  theme_classic() +
  labs(x = NULL, y = "Hyphal density (m/cm³, backtransformed)", fill = "Crop status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  annotate("text", x = 1, y = max(combined_df$ymax, na.rm = TRUE) * 1.05,
           label = "Global", fontface = "bold", hjust = 0.5) +
  annotate("text", x = length(unique(combined_df$biome))/2 + 1,
           y = max(combined_df$ymax, na.rm = TRUE) * 1.05,
           label = "Biomes", fontface = "bold", hjust = 0.5)

print(p_density)

# === Plot 2: Sample sizes
p_samples <- df %>%
  count(biome, crop) %>%
  ggplot(aes(x = biome, y = n, fill = crop)) +
  geom_col(position = position_dodge(0.6), width = 0.6, color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#CCCCCC", "#2166AC")) +
  labs(x = NULL, y = "Number of samples", fill = "Crop status") +
  ggtitle("Sample sizes by biome and crop") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_samples)

fig_dir <- "/"
ggsave("HyphalDensity_Global_Biome_SE_log.pdf", p_density, path = fig_dir, width = 12, height = 8, scale = 0.7)
ggsave("SampleCounts_by_Biome_and_Crop_log.pdf", p_samples, path = fig_dir, width = 10, height = 6, scale = 0.7)
