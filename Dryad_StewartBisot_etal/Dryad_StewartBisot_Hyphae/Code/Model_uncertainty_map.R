#Description: 
#This script compares the model standard deviation to the RMSE by computing their ratio for each #pixel in a global hyphal density map. The rasters are coarsened to reduce processing time, #reclassified into two uncertainty classes, converted to a data frame, and mapped globally. The #output highlights regions where uncertainty is lower or higher than the model RMSE.


##### Libraries & Setup #####
library(terra)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(progress)

# === Load rasters ===
mean_raster <- rast("/hyphal_density_m_cm3_Classified_mean.tif")
sd_relative <- rast("hyphal_density_m_cm3_Classified_sd.tif")
rmse_val <- 0.31  # m/cm³


# === Coarsen map to speed up 
fact <- 60  # 0.0083 * 60 ~ 0.5 degree
mean_coarse <- aggregate(mean_raster, fact = fact, fun = mean, na.rm = TRUE)
sd_rel_coarse <- aggregate(sd_relative, fact = fact, fun = mean, na.rm = TRUE)
gc()

# === Compute SD/RMSE binary classification ===
sd_abs <- mean_coarse * sd_rel_coarse
sd_ratio <- sd_abs / rmse_val
sd_ratio[sd_ratio <= 0] <- NA
sd_binary <- classify(sd_ratio, rcl = matrix(c(-Inf, 1, 1,
                                               1,  Inf, 2), 
                                             ncol = 3, byrow = TRUE), 
                      include.lowest = TRUE)


# === Convert to data frame for fast plot ===
df <- as.data.frame(sd_binary, xy = TRUE, na.rm = TRUE)
colnames(df)[3] <- "class"
df$class <- factor(df$class, levels = c(1, 2), labels = c("SD / RMSE < 1", "SD / RMSE ≥ 1"))

# === Load country shapefile ===
world <- ne_countries(scale = "small", returnclass = "sf")

# === Global Plot ===
ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = class)) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.1) +
  scale_fill_manual(values = c("green", "purple")) +
  coord_sf(expand = FALSE) +
  theme_classic() +
  labs(title = "SD / RMSE Ratio (Global)", fill = NULL) +
  theme(legend.position = "top")
ggsave("SD_RMSE_binary_global_ggplot.pdf", width = 11, height = 6, dpi = 300,device = "pdf")
gc()
