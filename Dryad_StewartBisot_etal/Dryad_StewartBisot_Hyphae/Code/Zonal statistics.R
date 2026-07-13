# =============================================================================
# Description:
# This script computes biome level, ecoregion level, and global summaries of
# predicted hyphal network density from raster outputs. The workflow aligns
# biome, ecoregion, and hyphal density rasters, converts predicted hyphal
# length from m cm⁻³ to km km⁻² to a 15 cm soil depth calculates zonal
# means, standard deviations, totals, and coefficients of variation for each
# biome and ecoregion. Global values are derived by aggregating across all
# biomes. The script produces summary tables for biomes, ecoregions, and the
# global estimate and exports them for downstream analysis.
# =============================================================================


##### Libraries & Setup #####
library(terra)
library(tidyverse)
library(sf)

# === Paths ===
biome_raster <- rast("ResolveBiome.tif")
ecoregion_raster <- rast("Ecoregion_GEE.tif")
length_raster <- rast("hyphal_density_m_cm3_Classified_mean.tif")
cov_raster <- rast("hyphal_density_m_cm3_Classified_cov.tif")
sd_raster <- rast("hyphal_density_m_cm3_Classified_sd.tif")
ecoregion_shp <- st_read("Ecoregions2017.shp")

# === Align rasters ===
biome_raster <- resample(biome_raster, length_raster, method = "near")
ecoregion_raster <- resample(ecoregion_raster, length_raster, method = "near")
sd_raster <- resample(sd_raster, length_raster, method = "bilinear")

# === Conversion constants ===
km <- 100000
depth <- 15
volume_length <- km * km * depth

##### A. Per-Biome Stats #####
length_image <- length_raster * volume_length
length_km_km2 <- length_image * 0.001
sd_total_km <- sd_raster * volume_length * 0.001

biome_mean_km_km2 <- zonal(length_km_km2, biome_raster, fun = mean, na.rm = TRUE)
biome_sd_km_km2 <- zonal(length_km_km2, biome_raster, fun = sd, na.rm = TRUE)
biome_sum_km <- zonal(length_km_km2, biome_raster, fun = sum, na.rm = TRUE)
biome_mean_m_cm3 <- zonal(length_raster, biome_raster, fun = mean, na.rm = TRUE)
biome_sd_m_cm3 <- zonal(length_raster, biome_raster, fun = sd, na.rm = TRUE)
biome_cov <- zonal(cov_raster, biome_raster, fun = mean, na.rm = TRUE)
biome_total_sd_km <- zonal(sd_total_km, biome_raster, fun = sum, na.rm = TRUE)

colnames(biome_mean_km_km2) <- c("BiomeID", "mean_km_km2")
colnames(biome_sd_km_km2) <- c("BiomeID", "sd_km_km2")
colnames(biome_sum_km) <- c("BiomeID", "total_km")
colnames(biome_mean_m_cm3) <- c("BiomeID", "mean_m_cm3")
colnames(biome_sd_m_cm3) <- c("BiomeID", "sd_m_cm3")
colnames(biome_cov) <- c("BiomeID", "cov")
colnames(biome_total_sd_km) <- c("BiomeID", "total_sd_km")

biome_labels <- ecoregion_shp %>%
  st_drop_geometry() %>%
  select(BIOME_NUM, BIOME_NAME) %>%
  distinct() %>%
  rename(BiomeID = BIOME_NUM)

biome_summary <- reduce(
  list(biome_mean_km_km2, biome_sd_km_km2, biome_sum_km, biome_mean_m_cm3, biome_sd_m_cm3, biome_cov, biome_total_sd_km),
  left_join,
  by = "BiomeID"
) %>%
  left_join(biome_labels, by = "BiomeID") %>%
  select(BIOME_NAME, mean_m_cm3, sd_m_cm3, cov, mean_km_km2, sd_km_km2, total_km, total_sd_km)

##### B. Global Stats #####
global_sum_km <- sum(biome_summary$total_km, na.rm = TRUE)
global_mean_km_km2 <- mean(biome_summary$mean_km_km2, na.rm = TRUE)
global_sd_km_km2 <- mean(biome_summary$sd_km_km2, na.rm = TRUE)
global_mean_m_cm3 <- mean(biome_summary$mean_m_cm3, na.rm = TRUE)
global_sd_m_cm3 <- mean(biome_summary$sd_m_cm3, na.rm = TRUE)
global_cov <- mean(biome_summary$cov, na.rm = TRUE)
global_total_sd_km <- sum(biome_summary$total_sd_km, na.rm = TRUE)

cat("\n🌍 Global Stats (based on biome-level aggregation):\n") #Ignoring RMSE 
cat("Total hyphal length (km, 15cm):", global_sum_km, "\n")
cat("Mean hyphal length (km/km², 15cm):", global_mean_km_km2, "\n")
cat("SD hyphal length (km/km², 15cm):", global_sd_km_km2, "\n")
cat("Mean hyphal density (m/cm³):", global_mean_m_cm3, "\n")
cat("SD hyphal density (m/cm³):", global_sd_m_cm3, "\n")
cat("Mean CoV:", global_cov, "\n")
cat("Global total hyphal length uncertainty (km):", global_total_sd_km, "\n")

global_summary <- tibble(
  total_km_15cm = global_sum_km,
  mean_km_km2_15cm = global_mean_km_km2,
  sd_km_km2_15cm = global_sd_km_km2,
  mean_m_cm3 = global_mean_m_cm3,
  sd_m_cm3 = global_sd_m_cm3,
  cov = global_cov,
  global_total_sd_km = global_total_sd_km
)

##### C. Per-Ecoregion Stats #####
eco_mean_km_km2 <- zonal(length_km_km2, ecoregion_raster, fun = mean, na.rm = TRUE)
eco_sd_km_km2 <- zonal(length_km_km2, ecoregion_raster, fun = sd, na.rm = TRUE)
eco_sum_km <- zonal(length_km_km2, ecoregion_raster, fun = sum, na.rm = TRUE)
eco_mean_m_cm3 <- zonal(length_raster, ecoregion_raster, fun = mean, na.rm = TRUE)
eco_sd_m_cm3 <- zonal(length_raster, ecoregion_raster, fun = sd, na.rm = TRUE)
eco_cov <- zonal(cov_raster, ecoregion_raster, fun = mean, na.rm = TRUE)
eco_total_sd_km <- zonal(sd_total_km, ecoregion_raster, fun = sum, na.rm = TRUE)

colnames(eco_mean_km_km2) <- c("ECO_ID", "mean_km_km2")
colnames(eco_sd_km_km2) <- c("ECO_ID", "sd_km_km2")
colnames(eco_sum_km) <- c("ECO_ID", "total_km")
colnames(eco_mean_m_cm3) <- c("ECO_ID", "mean_m_cm3")
colnames(eco_sd_m_cm3) <- c("ECO_ID", "sd_m_cm3")
colnames(eco_cov) <- c("ECO_ID", "cov")
colnames(eco_total_sd_km) <- c("ECO_ID", "total_sd_km")

ecoregion_df <- st_drop_geometry(ecoregion_shp) %>%
  select(ECO_ID, BIOME_NAME) %>%
  distinct()

ecoregion_summary <- reduce(
  list(eco_mean_km_km2, eco_sd_km_km2, eco_sum_km, eco_mean_m_cm3, eco_sd_m_cm3, eco_cov, eco_total_sd_km),
  left_join,
  by = "ECO_ID"
) %>%
  left_join(ecoregion_df, by = "ECO_ID") %>%
  select(ECO_ID, BIOME_NAME, mean_m_cm3, sd_m_cm3, cov, mean_km_km2, sd_km_km2, total_km, total_sd_km)

##### 📤 Export Area #####
write_csv(biome_summary, "/Biome_Hyphal_Summary_update.csv")
write_csv(ecoregion_summary, "coregion_Hyphal_Summary.csv")
write_csv(global_summary, "/Global_Hyphal_Summary.csv")

