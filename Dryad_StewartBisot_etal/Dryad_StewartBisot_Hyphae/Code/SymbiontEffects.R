# =============================================================================
# Description:
# This script evaluates differences in hyphal density among arbuscular
# mycorrhizal symbiont genera. The workflow loads and filters a harmonized
# dataset, updates genus level taxonomy, and fits a linear mixed model with
# symbiont genus as the fixed effect and study DOI as a random intercept.
# A null model is used for comparison. The script generates predicted effects
# for each genus, computes estimated marginal means with standard errors, and
# exports figures summarizing model predictions.
# =============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(dplyr)
library(ggplot2)
library(forcats)

# === Load and clean data ===
df <- read.csv("hyphal_density_full_metadata")
df$Hyphae_m_cm3 <- df$HyphalLength_m_cm3_final
df$Symbiont_Genus <- ifelse(df$symbiont_genus == "" | is.na(df$symbiont_genus), "Unclassified", df$symbiont_genus)
df$Symbiont_Genus[df$symbiont_genus == "Claroideoglomus"] <- "Entrophosphora" #Update taxonomy


# === Filter and clip ===
df <- df %>%
  filter(!is.na(DOI), !is.na(Hyphae_m_cm3), !is.na(Symbiont_Genus)) %>%
  filter(Hyphae_m_cm3 <= quantile(Hyphae_m_cm3, 0.95, na.rm = TRUE))

# === Remove Diversispora, small sample size ===
df <- df %>%
  filter(Symbiont_Genus != "Diversispora")

# === Model: symbiont genus only ===
mod_symb <- lmer(Hyphae_m_cm3 ~ Symbiont_Genus + (1 | DOI), data = df)
mod_symb_null <- lmer(Hyphae_m_cm3 ~ 1 + (1 | DOI), data = df)

# === Model comparison ===
anova(mod_symb_null, mod_symb)

# === Plot marginal effects ===
p <- plot_model(mod_symb, type = "pred", terms = "Symbiont_Genus", title = "Predicted hyphal density by symbiont genus")

ggsave("Pred_Hyphae_SymbiontGenus.pdf", p, width = 7, height = 5, dpi = 300)

# === Estimated marginal means with SE ===
emm_df <- emmeans(mod_symb, ~ Symbiont_Genus) %>%
  as.data.frame() %>%
  mutate(Symbiont_Genus = fct_reorder(Symbiont_Genus, emmean))

# === Plot EMMs ===
p_emm <- ggplot(emm_df, aes(x = Symbiont_Genus, y = emmean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "grey70", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  coord_flip() +
  theme_classic(base_size = 12) +
  ylab("Predicted hyphal density (m/cm³)") +
  xlab("") +
  ggtitle("EMMs ± SE: Hyphal density by symbiont genus")

ggsave("EMM_Hyphae_SymbiontGenus.pdf", p_emm, dpi = 300, width = 7, height = 5)
