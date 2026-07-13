# =============================================================================
# Description:
# This script analyzes observed hyphal density data to test whether values differ
# across plant growth forms, cultivation status, and their interaction. The
# workflow loads and filters the harmonized dataset, standardizes taxonomy and
# growth form labels, and fits mixed effects models with study DOI as a random
# intercept. Separate models evaluate effects of cultivation, growth form, and
# a combined growth form by cultivation interaction. The script generates
# predicted marginal effects, estimated marginal means with standard errors,
# and raw data summaries, and exports all visualizations for downstream use. 
# Plant growth form and cultivation was classified using Kew's Plants of the World Database
# =============================================================================


library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(dplyr)
library(ggplot2)
library(forcats)

# === Load and clean data ===
df <- read.csv("/Users/justinstewart/Dropbox/Collaborations/HyphaeBiomass/AM_networkdensity_github/Revision/Dryad/Databases/Hyphae/hyphal_density_full_metadata.csv")

df$Hyphae_m_cm3 <- df$HyphalLength_m_cm3_final
df$Growth_Form_Final <- ifelse(df$Growth_Form_Final == "" | is.na(df$Growth_Form_Final), "Mixed community", df$Growth_Form_Final)

# === Filter and clip extreme values ===
df <- df %>%
  filter(!is.na(DOI), !is.na(Hyphae_m_cm3), !is.na(cultivated_wild_Pironon2023)) %>%
  filter(Hyphae_m_cm3 <= quantile(Hyphae_m_cm3, 0.95, na.rm = FALSE))

# === Format cultivation column ===
df$cultivated_wild_Pironon2023[df$Growth_Form_Final == "Mixed community"] <- "Wild"
df$cultivated_wild_Pironon2023 <- tools::toTitleCase(tolower(df$cultivated_wild_Pironon2023))

# === Create combined variable ===
df$GrowthCult <- interaction(df$Growth_Form_Final, df$cultivated_wild_Pironon2023, sep = "_")

# === MODELS ===

# Cultivation only
mod_cult <- lmer(Hyphae_m_cm3 ~ cultivated_wild_Pironon2023 + (1 | DOI), data = df)
mod_cult_null <- lmer(Hyphae_m_cm3 ~ 1 + (1 | DOI), data = df)

# Growth form only
mod_growth <- lmer(Hyphae_m_cm3 ~ Growth_Form_Final + (1 | DOI), data = df)
mod_growth_null <- lmer(Hyphae_m_cm3 ~ 1 + (1 | DOI), data = df)

# Combined growth form × cultivation
mod_combo <- lmer(Hyphae_m_cm3 ~ GrowthCult + (1 | DOI), data = df)
mod_combo_null <- lmer(Hyphae_m_cm3 ~ 1 + (1 | DOI), data = df)

# === Model comparisons ===
anova(mod_cult_null, mod_cult)
anova(mod_growth_null, mod_growth)
anova(mod_combo_null, mod_combo)

# === Save marginal effect plots ===
p1 <- plot_model(mod_cult, type = "pred", title = "Predicted hyphal density by cultivation")
ggsave("/Pred_Hyphae_Cultivated_vs_Wild.png", p1, width = 6, height = 4, dpi = 300)

p2 <- plot_model(mod_growth, type = "pred", terms = "Growth_Form_Final", title = "Predicted hyphal density by growth form")
ggsave("Pred_Hyphae_GrowthForm.png", p2, width = 7, height = 5, dpi = 300)

p3 <- plot_model(mod_combo, type = "pred", terms = "GrowthCult", title = "Predicted hyphal density by GrowthForm × Cultivation")
ggsave("Pred_Hyphae_GrowthForm_Cultivation.png", p3, width = 9, height = 5, dpi = 300)

library(emmeans)
library(ggplot2)
library(dplyr)
library(forcats)

# === Estimated marginal means for mod_combo ===
emm_df <- emmeans(mod_combo, ~ GrowthCult) %>%
  as.data.frame() %>%
  mutate(
    Growth_Form = sub("_.*", "", GrowthCult),
    Cultivation = sub(".*_", "", GrowthCult),
    Growth_Form = fct_reorder(Growth_Form, emmean)
  )

# Capitalize cultivation
emm_df$Cultivation <- tools::toTitleCase(tolower(emm_df$Cultivation))


# === Plot EMMs with SEs ===
p4 <- ggplot(emm_df, aes(x = Growth_Form, y = emmean, fill = Cultivation)) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    alpha = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(
    values = c(
      "Cultivated" = "#E69F00",  # orange
      "Wild" = "#0072B2"         # blue
    )
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  ylab("Predicted hyphal density (m/cm³)") +
  xlab("") +
  ggtitle("EMMs ± SE: Hyphal density by growth form and cultivation") +
  theme(legend.position = "top")

ggsave("Pred_Hyphae_GrowthForm_Cultivation.pdf", p4, device="pdf", dpi = 300,scale=0.5)

# === Summary of raw data by group (mean ± SD) ===
summary_raw <- df %>%
  group_by(Growth_Form_Final, cultivated_wild_Pironon2023) %>%
  summarise(
    n = n(),
    mean = mean(Hyphae_m_cm3, na.rm = TRUE),
    sd = sd(Hyphae_m_cm3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cultivation = tools::toTitleCase(tolower(cultivated_wild_Pironon2023)),
    Growth_Form = fct_reorder(Growth_Form_Final, mean)
  )

# === Plot: Raw data summary (mean ± SD) ===
p_raw <- ggplot(summary_raw, aes(x = Growth_Form, y = mean, fill = Cultivation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Cultivated" = "#D73027", "Wild" = "#4575B4")) +
  coord_flip() +
  theme_classic(base_size = 12) +
  ylab("Observed hyphal density (m/cm³)") +
  xlab("") +
  ggtitle("Raw mean ± SD: Hyphal density by growth form and cultivation") +
  theme(legend.position = "top")

# === Save plot ===
ggsave("Raw_Hyphae_GrowthForm_Cultivation.pdf", p_raw, dpi = 300, width = 9, height = 5)


