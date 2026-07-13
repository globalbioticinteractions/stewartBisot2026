library(dplyr)
library(ggplot2)
library(maps)
library(patchwork)
library(cowplot)

# Load data
data <- read.csv("hyphal_density_full_metadata.csv")

data$pot_field <- tolower(trimws(data$pot_field))
data$Latitude <- data$latitude
data$Longitude <- data$longitude

# World map
world <- map_data("world")

# Split data
pot_data <- data %>%
  filter(
    pot_field == "pot",
    !is.na(Latitude),
    !is.na(Longitude),
    !is.na(HyphalLength_m_cm3_final)
  )

field_data <- data %>%
  filter(
    pot_field == "field",
    !is.na(Latitude),
    !is.na(Longitude),
    !is.na(HyphalLength_m_cm3_final)
  )

# Maps
map_pot <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = pot_data,
    aes(x = Longitude, y = Latitude),
    color = "#8CC891",
    size = 1.4,
    alpha = 0.7
  ) +
  coord_fixed(1.3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Pot studies",
    x = "Longitude",
    y = "Latitude"
  )

map_field <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = field_data,
    aes(x = Longitude, y = Latitude),
    color = "#B41A21",
    size = 1.4,
    alpha = 0.7
  ) +
  coord_fixed(1.3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Field studies",
    x = "Longitude",
    y = "Latitude"
  )

# Histograms
hist_pot <- ggplot(pot_data, aes(x = log(HyphalLength_m_cm3_final))) +
  geom_histogram(bins = 30, fill = "#8CC891", color = "black", alpha = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Pot studies",
    x = "Hyphal density (log m/cm³)",
    y = "Count"
  )

hist_field <- ggplot(field_data, aes(x = log(HyphalLength_m_cm3_final))) +
  geom_histogram(bins = 30, fill = "#B41A21", color = "black", alpha = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Field studies",
    x = "Hyphal density (log m/cm³)",
    y = "Count"
  )

# Combine layout with panel labels A–D
final_plot <- (map_pot | map_field) / (hist_pot | hist_field) +
  plot_annotation(tag_levels = "A")

ggsave(
  "pot_field_maps_and_histograms.pdf",
  plot = final_plot,
  width = 13,
  height = 8,
  dpi = 300
)
