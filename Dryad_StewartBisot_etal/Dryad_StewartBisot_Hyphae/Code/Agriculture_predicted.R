# =============================================================================
# Description:
# This script quantifies differences in predicted hyphal network density between
# cropland and non cropland areas across global biomes. The workflow loads mean
# and uncertainty rasters, applies cropland masks, computes zonal statistics by
# biome, and calculates the difference in hyphal density between land use types
# together with propagated standard deviations. The script produces global and
# per biome summaries, generates visualizations of absolute and percent
# differences, and estimates the total amount of hyphal network that could be
# restored under a non cropland baseline.
# =============================================================================

library(terra)
library(dplyr)
library(sf)
library(ggplot2)

# Load rasters
hyphae <- rast("hyphal_density_m_cm3_Classified_mean.tif")
hyphae_sd_rast <- rast("/hyphal_density_m_cm3_Classified_sd.tif")
biome_raster <- rast("ResolveBiome.tif")
esa <- rast("/ESA_Worldcover.tif")

# Convert to km/km² (same as before)
km <- 100000; depth <- 15; volume_length <- km * km * depth
hyphae_km <- hyphae * volume_length * 0.001
hyphae_sd_km <- hyphae_sd_rast * volume_length * 0.001

# Resample
esa_resampled <- resample(esa, hyphae, method = "near")
biome_resampled <- resample(biome_raster, hyphae, method = "near")

# Crop/non-crop masks
crop_mask <- esa_resampled; crop_mask[crop_mask != 40] <- NA; crop_mask[crop_mask == 40] <- 1
noncrop_mask <- esa_resampled; noncrop_mask[noncrop_mask == 40] <- NA; noncrop_mask[!is.na(noncrop_mask)] <- 1

# Masked mean and SD rasters
hyphae_crop <- mask(hyphae_km, crop_mask, maskvalues = NA)
hyphae_noncrop <- mask(hyphae_km, noncrop_mask, maskvalues = NA)
hyphae_sd_crop <- mask(hyphae_sd_km, crop_mask, maskvalues = NA)
hyphae_sd_noncrop <- mask(hyphae_sd_km, noncrop_mask, maskvalues = NA)

# Zonal mean and SD per biome
crop_mean    <- as.data.frame(zonal(hyphae_crop, biome_resampled, fun = "mean", na.rm = TRUE))
noncrop_mean <- as.data.frame(zonal(hyphae_noncrop, biome_resampled, fun = "mean", na.rm = TRUE))
crop_sd      <- as.data.frame(zonal(hyphae_sd_crop, biome_resampled, fun = "mean", na.rm = TRUE))
noncrop_sd   <- as.data.frame(zonal(hyphae_sd_noncrop, biome_resampled, fun = "mean", na.rm = TRUE))
names(crop_mean)    <- c("biome_id", "hyphae_mean_crop")
names(noncrop_mean) <- c("biome_id", "hyphae_mean_noncrop")
names(crop_sd)      <- c("biome_id", "hyphae_sd_crop_km")
names(noncrop_sd)   <- c("biome_id", "hyphae_sd_noncrop_km")

# Ecoregion/biome lookup
ecoregions <- st_read("/Ecoregions2017.shp") # https://ecoregions.appspot.com/ 
biome_lookup <- ecoregions %>% st_drop_geometry() %>% distinct(BIOME_NUM, BIOME_NAME)

# Combine all and calculate differences
results <- crop_mean %>%
  left_join(crop_sd, by = "biome_id") %>%
  left_join(noncrop_mean, by = "biome_id") %>%
  left_join(noncrop_sd, by = "biome_id") %>%
  left_join(biome_lookup, by = c("biome_id" = "BIOME_NUM")) %>%
  mutate(
    mean_diff = hyphae_mean_noncrop - hyphae_mean_crop,
    sd_diff = sqrt(hyphae_sd_crop_km^2 + hyphae_sd_noncrop_km^2)
  )

# Filter to valid biomes
results_plot <- results %>%
  filter(!is.na(BIOME_NAME), BIOME_NAME != "N/A", !grepl("Tundra", BIOME_NAME, ignore.case=TRUE)) %>%
  mutate(
    effect = case_when(
      mean_diff > 0 ~ "Lost (restorable)",
      mean_diff < 0 ~ "Gain (under cropland)"
    )
  )

# Global mean and SD bar
global_row <- data.frame(
  BIOME_NAME = "Global",
  mean_diff = mean(results_plot$mean_diff, na.rm = TRUE),
  sd_diff = sqrt(sum(results_plot$sd_diff^2, na.rm = TRUE)) / nrow(results_plot),
  effect = "Lost (restorable)"
)

plot_df <- bind_rows(global_row, results_plot)

# Order for plotting: global first, then biomes descending by mean_diff
biome_levels <- plot_df %>%
  filter(BIOME_NAME != "Global") %>%
  arrange(desc(mean_diff)) %>%
  pull(BIOME_NAME)
plot_df$BIOME_NAME <- factor(plot_df$BIOME_NAME, levels = c("Global", biome_levels))

# Convert to million km for axis
plot_df <- plot_df %>%
  mutate(
    mean_diff = mean_diff / 1e6,
    sd_diff = sd_diff / 1e6
  )

max_y <- max(plot_df$mean_diff + plot_df$sd_diff, na.rm = TRUE)

# Plot using the per-biome SD as error bars
p <- ggplot(plot_df, aes(x = BIOME_NAME, y = mean_diff, fill = effect)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2, colour = "black") +
  scale_fill_manual(values = c("Lost (restorable)" = "steelblue", "Gain (under cropland)" = "tomato")) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean Difference (million km hyphae/km²)\n(Non-crop - Crop, error bars = SD)",
    fill = "Effect"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 1, color = "black") +
  annotate("text", x = 1, y = max_y * 1.05,
           label = "Global", fontface = "bold", hjust = 0.5, size = 5) +
  annotate("text", x = (length(levels(plot_df$BIOME_NAME)) + 1)/2 + 0.5,
           y = max_y * 1.05,
           label = "Biomes", fontface = "bold", hjust = 0.5, size = 5)

fig_dir <- ""
ggsave("HyphalDensity_Difference_Global_and_Biome_PerBiomeSD.pdf", p, path = fig_dir, width = 12, height = 7, scale = 0.7)
print(p)


# Total change in hyphae length based on agriculture (global vaue)
# Calculate area_crop_km2 if not already in your results (run only if missing)
pixel_area_km2 <- cellSize(hyphae, unit = "km")
crop_area_km2 <- mask(pixel_area_km2, crop_mask, maskvalues = NA)
crop_area_df <- as.data.frame(zonal(crop_area_km2, biome_resampled, fun = "sum", na.rm = TRUE))
names(crop_area_df) <- c("biome_id", "area_crop_km2")
results_plot <- results_plot %>% left_join(crop_area_df, by = "biome_id")

# Multiply mean_diff and sd_diff by area_crop_km2 to get total for each biome, mean_diff and sd_diff are in km/km², area_crop_km2 in km², so product is in km)
biome_totals <- results_plot %>%
  mutate(
    total_restore = mean_diff * area_crop_km2,         # km hyphae (may be negative for gain)
    total_restore_sd = sd_diff * area_crop_km2         # SD in km hyphae
  )

# Sum across biomes for the global estimate (ignore NA or Tundra, etc.)
global_total_restore <- sum(biome_totals$total_restore, na.rm = TRUE)
global_total_restore_sd <- sqrt(sum(biome_totals$total_restore_sd^2, na.rm = TRUE)) # propagate SDs as independent

cat("Global restoration potential: ", signif(global_total_restore, 3), " km hyphae\n")
cat("Global uncertainty (SD): ", signif(global_total_restore_sd, 3), " km hyphae\n")
# Save table if desired
write.csv(biome_totals, file.path(fig_dir, "HyphalDensity_PerBiome.csv"), row.names = FALSE)


library(tidyverse)

# === Load data table exported above ===
biome_totals <- read.csv("/HyphalDensity_PerBiome.csv")

# === Calculate percent difference and SD ===
plot_df <- biome_totals %>%
  mutate(
    percent_diff = 100 * (hyphae_mean_noncrop - hyphae_mean_crop) / hyphae_mean_crop,
    percent_sd = 100 * sqrt(hyphae_sd_crop_km^2 + hyphae_sd_noncrop_km^2) / hyphae_mean_crop,
    effect = ifelse(percent_diff >= 0, "Lost (restorable)", "Gain (under cropland)")
  )

# === Add Global row (mean across all biomes) ===
global_row <- tibble(
  BIOME_NAME = "Global",
  percent_diff = mean(plot_df$percent_diff, na.rm = TRUE),
  percent_sd = mean(plot_df$percent_sd, na.rm = TRUE),
  effect = "Lost (restorable)"
)

plot_df <- bind_rows(global_row, plot_df)

# === Reorder BIOME_NAME: Global first, then Lost (descending), then Gain (ascending) ===
plot_df <- plot_df %>%
  arrange(factor(BIOME_NAME == "Global", levels = c(TRUE, FALSE)),
          desc(effect == "Lost (restorable)"),
          if_else(effect == "Lost (restorable)", -percent_diff, percent_diff)) %>%
  mutate(BIOME_NAME = factor(BIOME_NAME, levels = BIOME_NAME))

max_y <- max(plot_df$percent_diff + plot_df$percent_sd, na.rm = TRUE)

# === Plot ===
p <- ggplot(plot_df, aes(x = BIOME_NAME, y = percent_diff, fill = effect)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = percent_diff - percent_sd, ymax = percent_diff + percent_sd),
                width = 0.2, colour = "black") +
  scale_fill_manual(values = c("Lost (restorable)" = "steelblue", "Gain (under cropland)" = "tomato")) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Percent Difference in Hyphal Density\n(Non-crop vs Crop, error bars = SD)",
    fill = "Effect"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 1, color = "black") +
  annotate("text", x = 1, y = max_y * 1.05,
           label = "Global", fontface = "bold", hjust = 0.5, size = 5) +
  annotate("text", x = (length(levels(plot_df$BIOME_NAME)) + 1)/2 + 0.5,
           y = max_y * 1.05,
           label = "Biomes", fontface = "bold", hjust = 0.5, size = 5)

print(p)

fig_dir <- ""
ggsave("HyphalDensity_Difference_Global_and_Biome_PerBiomeSD_Percent.pdf", p, path = fig_dir, width = 12, height = 7, scale = 0.7)


# === Plot ===
p <- ggplot(plot_df, aes(x = BIOME_NAME, y = percent_diff, fill = effect)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = percent_diff - percent_sd, ymax = percent_diff + percent_sd),
                width = 0.2, colour = "black") +
  scale_fill_manual(values = c("Lost (restorable)" = "steelblue", "Gain (under cropland)" = "tomato")) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Percent Difference in Hyphal Density\n(Non-crop vs Crop, error bars = SD)",
    fill = "Effect"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 1, color = "black") +
  annotate("text", x = 1, y = max_y * 1.05,
           label = "Global", fontface = "bold", hjust = 0.5, size = 5) +
  annotate("text", x = (length(levels(plot_df$BIOME_NAME)) + 1)/2 + 0.5,
           y = max_y * 1.05,
           label = "Biomes", fontface = "bold", hjust = 0.5, size = 5)

print(p)


