# =============================================================================
# Description:
# This script trains a Random Forest model to predict observed hyphal density
# and computes SHAP values to quantify the influence of environmental and
# study level variables. The workflow loads and filters the dataset, prepares
# predictor variables, fits a Random Forest model, and computes SHAP values
# using fastshap. It generates custom dependence plots, beeswarm plots, and
# variable importance summaries, and produces additional visualizations for
# seasonal effects. The script also exports SHAP based plots for selected
# predictors and saves importance and beeswarm figures for further analysis.
# =============================================================================


library(tidyverse)
library(ranger)
library(fastshap)
library(shapviz)

# === Load and filter data ===
data <- read.csv("hyphal_density_m_cm3.csv") 
# Log-transform requested covariates
data <- data %>% #Match model
  mutate(
    nitrogen = log(nitrogen),
    phosphorus = log(phosphorus)) 


# Define predictors
core_covs <- c("CGIAR_PET", "CHELSA_BIO_Annual_Mean_Temperature", "CHELSA_BIO_Annual_Precipitation",
               "CHELSA_BIO_Max_Temperature_of_Warmest_Month", "CHELSA_BIO_Precipitation_Seasonality",
               "EarthEnvTexture_CoOfVar_EVI", "EarthEnvTexture_Correlation_EVI", "EarthEnvTexture_Homogeneity_EVI",
               "EarthEnvTopoMed_AspectCosine", "EarthEnvTopoMed_AspectSine", "EarthEnvTopoMed_Elevation",
               "EarthEnvTopoMed_Slope", "EarthEnvTopoMed_TopoPositionIndex", "EsaCci_BurntAreasProbability",
               "GHS_Population_Density", "MODIS_NPP", "SG_Depth_to_bedrock", "SG_SOC_Content_005cm",
               "SG_Soil_pH_H2O_005cm", "crops_ESA", "cultivated_grassland_2010_2022", "earthEnvLandcover_class_barren",
               "earthEnvLandcover_class_cultivated_managed_vegetation", "earthEnvLandcover_class_evergreen_broadleaf_trees",
               "earthEnvLandcover_class_evergreen_decid_needleleaf_trees", "earthEnvLandcover_class_herbaceous_vegetation",
               "earthEnvLandcover_class_mixed_other_trees", "earthEnvLandcover_class_regularly_flooded_vegetation",
               "earthEnvLandcover_class_shrubs", "harmonized_aboveground_biomass", "harmonized_belowground_biomass",
               "isric_soil_proportion_of_sand", "myco_veg_cover", "nitrogen", "phosphorus",
               "sampling_intensity", "soil_depth_depth_cm", "sample_size", "soil_core_volume_cm3","Growth_Form_Grass","Growth_Form_Herb","Growth_Form_Shrub","Growth_Form_Tree", "cultivated1_wild0_Pironon2023"
)

season_vars <- grep("Season_NorthernHemHarmonized", colnames(data), value = TRUE)
predictors <- c(core_covs, season_vars)
data_model <- data %>%
  select(all_of(predictors), hyphal_density_m_cm3) %>%
  drop_na()

# === Train model ===
set.seed(123)
rf_model <- ranger(log(hyphal_density_m_cm3) ~ ., data = data_model, num.trees = 250)
print(rf_model)

# === SHAP computation ===
X <- data_model %>% select(-hyphal_density_m_cm3)

shap_matrix <- fastshap::explain(rf_model, X = X, pred_wrapper = function(object, newdata) {
  predict(object, data = newdata)$predictions
}, nsim = 1000)

shv <- shapviz(shap_matrix, X = X)
df_shaps <- data.frame(shv$S)
df_input <- data.frame(shv$X)

# === Custom SHAP plot ===
plot_shap_custom <- function(var, label = NULL) {
  df_plot <- data.frame(
    SHAP = df_shaps[[var]],
    Input = df_input[[var]]
  )
  
  ggplot(df_plot, aes(x = Input, y = SHAP)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(shape = 21, fill = "#98BFE2", color = "black", size = 3, alpha = 1) +
    stat_smooth(method = "loess", color = "black", lwd = 1.5, span = 2, se = TRUE) +
    theme_classic() +
    xlab(label %||% var) + ylab("SHAP") +
    ggtitle(paste("SHAP:", label %||% var))
}

vars_to_plot <- c("crops_ESA", season_vars, "soil_core_volume_cm3")
plot_list <- lapply(vars_to_plot, plot_shap_custom)
names(plot_list) <- vars_to_plot
for (p in plot_list) print(p)

# === Summary and beeswarm ===
sv_importance(shv, show_numbers = TRUE, max_display = Inf)
sv_importance(shv, kind = "beeswarm")
sv_dependence(shv, v = "phosphorus")

# === Season-specific SHAPs ===
valid_season_vars <- season_vars[!str_detect(season_vars, "NA")]

season_df <- df_shaps %>%
  select(all_of(valid_season_vars)) %>%
  mutate(obs = row_number()) %>%
  pivot_longer(cols = -obs, names_to = "Season", values_to = "SHAP") %>%
  left_join(
    df_input %>%
      select(all_of(valid_season_vars)) %>%
      mutate(obs = row_number()) %>%
      pivot_longer(cols = -obs, names_to = "Season", values_to = "value"),
    by = c("obs", "Season")
  ) %>%
  filter(value == 1) %>%
  mutate(Season = str_replace(Season, "Season_NorthernHemHarmonized_", ""))

ggplot(season_df, aes(x = reorder(Season, SHAP), y = SHAP, fill = Season)) +
  geom_jitter(width = 0.2, shape = 21, color = "black", alpha = 0.7, size = 3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  ylab("SHAP value") +
  xlab("Season") +
  theme(legend.position = "none") +
  ggtitle("Seasonality")

# === Create plot ===
season_plot <- ggplot(season_df, aes(x = reorder(Season, SHAP), y = SHAP)) +
  geom_jitter(width = 0.2, shape = 21, color = "black",fill="#98BFE2", alpha = 0.7, size = 3) +
  stat_summary(fun = median, geom = "point", size=5, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  ylab("SHAP value") +
  xlab("Season") +
  theme(legend.position = "none") +
  ggtitle("Seasonality")

# === Save plot ===
ggsave(
  filename = "Seasonality.pdf",
  plot = season_plot,device = "pdf",
  scale = 0.7,
  width = 6,
  height = 5)


# === Export SHAP plots ===
save_shap_plots <- function(var_list, shap_df, input_df, out_dir = ".", width = 6, height = 5, dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE)
  for (var in var_list) {
    df_plot <- data.frame(
      SHAP = shap_df[[var]],
      Input = input_df[[var]]
    )
    
    p <- ggplot(df_plot, aes(x = Input, y = SHAP)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(shape = 21, fill = "#98BFE2", color = "black", size = 3, alpha = 1) +
      stat_smooth(method = "loess", color = "black", lwd = 1.5, span = 2, se = TRUE) +
      theme_classic() +
      xlab(var) + ylab("SHAP") +
      ggtitle(paste("SHAP:", var)) + ylim(-0.3,0.5)
    
    ggsave(filename = file.path("", paste0("SHAP_", var, ".pdf")), plot = p,
           width = width, height = height, dpi = dpi, scale = 0.7)
  }
}

vars_to_export <- c("") #Input specific variable to export
save_shap_plots(vars_to_export, df_shaps, df_input, out_dir = "SHAP_Plots")


# === Export beeswarm and SHAP importance plots ===
library(ggplot2)
library(patchwork)

# Generate plots
p1 <- sv_importance(shv, kind = "beeswarm", show_numbers = FALSE, max_display = 100) +
  ggtitle("SHAP Beeswarm Plot")

p2 <- sv_importance(shv, kind = "bar", show_numbers = TRUE, max_display = 100) +
  ggtitle("Mean Absolute SHAP")

# Save plots
ggsave(
  filename = "shap_beeswarm_plot.png",
  plot = p1, dpi = 300, scale = 0.9
)

 ggsave(
  filename = "shap_importance_plot.png",
  plot = p2, dpi = 300, scale = 0.7
)

