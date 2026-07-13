# =============================================================================
# Description:
This script loads a dataset of hyphal density measurements, fits three alternative depth decay models (linear, exponential, logarithmic) with biome and study as random effects, selects the best model by AIC, and calculates the half decay depth with confidence intervals using the delta method. It also produces a plot of the best fitting model with the estimated half decay depth highlighted.
# =============================================================================

# ===============================================================
# Libraries
# ===============================================================
library(tidyverse)
library(lme4)
library(sjPlot)

# ===============================================================
# Load and clean data
# ===============================================================
data <- read.csv(
  "hyphal_density_m_cm3.csv")

data_mod <- data %>% #Sanity check
  filter(hyphal_density_m_cm3 > 0,
         soil_depth_depth_cm > 0)

# ===============================================================
# Fit decay models
# ===============================================================
m_linear <- lmer(hyphal_density_m_cm3 ~ soil_depth_depth_cm +
                   (1|Resolve_Biome) + (1|Title), 
                 data = data_mod)

m_exp <- lmer(log(hyphal_density_m_cm3) ~ soil_depth_depth_cm +
                (1|Resolve_Biome) + (1|Title), 
              data = data_mod)

m_log <- lmer(hyphal_density_m_cm3 ~ log(soil_depth_depth_cm) +
                (1|Resolve_Biome) + (1|Title), 
              data = data_mod)

# ===============================================================
# Compare models
# ===============================================================
model_list <- list(Linear = m_linear, Exponential = m_exp, Logarithmic = m_log)
aic_table <- data.frame(Model = names(model_list),
                        AIC = sapply(model_list, AIC)) %>%
  arrange(AIC)
print(aic_table)

best_model_name <- aic_table$Model[1]
cat("Best model by AIC:", best_model_name, "\n")
best_model <- model_list[[best_model_name]]
anova(m_linear,m_exp,m_log)
# ===============================================================
# Half-decay depth + Delta method SE and CI
# ===============================================================
vc <- vcov(best_model)

if(best_model_name == "Exponential"){
  alpha <- fixef(best_model)[["(Intercept)"]]
  beta  <- fixef(best_model)[["soil_depth_depth_cm"]]
  half_depth <- log(0.5) / beta
  
  # Delta method
  d_dbeta <- -log(0.5) / (beta^2)
  var_beta <- vc["soil_depth_depth_cm","soil_depth_depth_cm"]
  se_half <- sqrt(d_dbeta^2 * var_beta)
  
} else if(best_model_name == "Linear"){
  alpha <- fixef(best_model)[["(Intercept)"]]
  beta  <- fixef(best_model)[["soil_depth_depth_cm"]]
  half_depth <- (alpha/2 - alpha) / beta
  
  # Gradient wrt alpha and beta
  d_dalpha <- (0.5 - 1)/beta
  d_dbeta  <- -(alpha/2 - alpha)/(beta^2)
  grad <- c(d_dalpha, d_dbeta)
  var_mat <- vc[c("(Intercept)","soil_depth_depth_cm"),
                c("(Intercept)","soil_depth_depth_cm")]
  se_half <- sqrt(t(grad) %*% var_mat %*% grad)
  
} else if(best_model_name == "Logarithmic"){
  alpha <- fixef(best_model)[["(Intercept)"]]
  beta  <- fixef(best_model)[["log(soil_depth_depth_cm)"]]
  y0 <- alpha  # at depth=1
  half_value <- y0/2
  half_depth <- exp((half_value - alpha)/beta)
  
  # Approximate SE with delta method
  # f(alpha,beta) = exp((y0/2 - alpha)/beta), y0 = alpha
  d_dalpha <- ( (0.5 - 1)/beta ) * half_depth
  d_dbeta  <- -((half_value - alpha)/(beta^2)) * half_depth
  grad <- c(d_dalpha, d_dbeta)
  var_mat <- vc[c("(Intercept)","log(soil_depth_depth_cm)"),
                c("(Intercept)","log(soil_depth_depth_cm)")]
  se_half <- sqrt(t(grad) %*% var_mat %*% grad)
}

ci95 <- c(half_depth - 1.96*se_half, half_depth + 1.96*se_half)

cat("Half-decay depth:", round(half_depth,2), "cm\n")
cat("SE:", round(se_half,2), "cm\n")
cat("95% CI:", round(ci95[1],2), "-", round(ci95[2],2), "cm\n")

# ===============================================================
# Plot with mean + CI for half-depth
# ===============================================================
p <- plot_model(best_model, type = "pred") +
  xlab("Soil depth (cm)") +
  ylab("Hyphal density (log m/cm³)") +
  ggtitle(paste("Best model:", best_model_name)) +
  geom_rect(aes(xmin = ci95[1], xmax = ci95[2], ymin = -Inf, ymax = Inf),
            fill = "blue", alpha = 0.01, inherit.aes = FALSE) +
  geom_vline(xintercept = half_depth, linetype = "dashed", color = "red") 

p

ggsave(
  filename = "Deepsoil_best_withCI.pdf",
  plot = p,
  width = 6,
  height = 5
)


