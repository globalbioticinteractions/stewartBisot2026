#Description: 
This script trains a Random Forest model to predict density patterns using a large set of environmental, soil, vegetation, land cover, and study level variables. The workflow loads the training data, prepares covariates, computes pairwise geographic distances among samples, and constructs spatial distance thresholds for evaluating spatial structure. It fits a nonspatial Random Forest model, tests a spatial version of the model with Moran eigenvector predictors, and performs spatial cross validation to quantify model performance. The script outputs model objects and spatial CV results that can be used to assess predictive accuracy and spatial autocorrelation in model residuals.

library(spatialRF)
library(tidyverse)
library(tidyverse)
library(ranger)
library(ggplot2)

# --- Load and prep ---
training <- read.csv("hyphal_density_m_cm3.csv")

training <- training %>% #Log transform to match model
  mutate(
    nitrogen   = log(nitrogen),
    phosphorus = log(phosphorus)
  )

# predictors
core_covs <- c(
  "CGIAR_PET", "CHELSA_BIO_Annual_Mean_Temperature", "CHELSA_BIO_Annual_Precipitation",
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
  "sampling_intensity", "soil_depth_depth_cm", "sample_size", "soil_core_volume_cm3",
  "Growth_Form_Grass", "Growth_Form_Herb", "Growth_Form_Shrub", "Growth_Form_Tree", 
  "cultivated1_wild0_Pironon2023"
)
season_vars <- grep("Season_NorthernHemHarmonized", colnames(data), value = TRUE) 
predictors <- c(core_covs, season_vars)
data<-data.frame(training)
predictor.variable.names <- training %>% dplyr::select(covariateList)
dependent.variable.name <- training %>% dplyr::select(hyphal_density_m_cm3)

#coordinates of the cases
xy <- data[, c("Pixel_Long", "Pixel_Lat")] %>% 
  rename(
    x = Pixel_Long,
    y = Pixel_Lat)

coords <- SpatialPoints(cbind(xy$x, xy$y), proj4string = CRS("+proj=longlat"))
dists <- pointDistance(coords, lonlat = FALSE)
dists_lonlat <- pointDistance(coords, lonlat = TRUE) #specify for lonlat
distance_matrix <- na.omit(dist(dists_lonlat) )

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(11000,50000,100000,500000,1000000,2000000,5000000)
distance.matrix <- dists_lonlat
distance.matrix[is.na(distance.matrix)] <- 0 
min(distance.matrix) 
max(distance.matrix)

#random seed for reproducibility
random.seed <- 5789
data<-data.frame(training)

model.non.spatial_full <- spatialRF::rf(
  data = data,
  n.cores = 4,
  dependent.variable.name = "hyphal_density_m_cm3",
  predictor.variable.names = covariateList,
  distance.matrix = (distance.matrix),
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = TRUE
)

# Construct a spatial model, decide against spatial model
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial_full,
  method = "mem.moran.sequential", #default method
  verbose = TRUE,
  n.cores = 4,
  max.spatial.predictors = 100,
  seed = random.seed
)

model.non.spatial_full_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial_full,
  xy = xy,       
  n.cores = 3,           #data coordinates
  repetitions = 10,         #number of spatial folds
  training.fraction = 0.7, #training data fraction on each fold
  metrics = "r.squared",
  seed = random.seed,
  verbose = TRUE
)

spatialRF::plot_evaluation(model.non.spatial_full_cv)

