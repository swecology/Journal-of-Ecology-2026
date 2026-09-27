##########################################################################################################
###### Drought legacy intensifies negative conspecific feedback in native but not invasive plants ########
############################# Correspondence: Wei Song (wei.song@unibe.ch) ###############################


library(lmerTest)
library(lme4) # the order must be like this (later than LmerTest), or there will be no p-values.
library(performance)
library(emmeans)
library(tidyverse)

data <- read.csv("Feedback data.csv", header = TRUE)
str(data)
data_mod <- data %>%
  mutate(
    Pot_nr            = as.numeric(Pot_nr),
    Species_status    = factor(Species_status, levels = c("Native", "Invasive")),
    Soil_conditioning = factor(Soil_conditioning, levels = c("Conspecific soil", 
                                                             "Native heterospecific soil", 
                                                             "Invasive heterospecific soil")),
    Water_treatment   = factor(Water_treatment, levels = c("Control", "Drought", 
                                                           "Flood", "Drought + Flood")),
    Replicates        = as.numeric(Replicates),
    Soil_origin_site  = as.character(Soil_origin_site)
  )
str(data_mod) # make sure that every variable is formatted correctly.

## effects on total biomass ## Fig.2
TBM <-lmer(sqrt(TBM) ~ Species_status * Soil_conditioning * Water_treatment  +
             (1|Species) + (1|Soil_origin_site), data = data_mod)
anova(TBM)

pairwise_TBM <- emmeans(TBM, pairwise ~ Water_treatment | Species_status * Soil_conditioning, 
                        adjust = "bonferroni")
print(pairwise_TBM$contrasts)
# use emmeans package to extract model fitting data and plot it.


## test invasional meltdown ## Fig.S7
meltdown <- emmeans(TBM, pairwise ~ Soil_conditioning | Species_status * Water_treatment, 
                         adjust = "bonferroni")
print(meltdown$contrasts)


## Cohen's d ## Fig.3
library(effsize)

emm_grid <- emmeans(TBM, ~ Soil_conditioning | Species_status * Water_treatment)
psf_effect_sizes <- eff_size(emm_grid, 
                             sigma = sigma(TBM), 
                             edf = df.residual(TBM))

psf_results_df <- as.data.frame(psf_effect_sizes)
psf_results_clean <- psf_results_df %>%
  filter(
    grepl("Native heterospecific|Invasive heterospecific", contrast) &
      grepl("Conspecific", contrast)
  ) %>%
  mutate(
    contrast_new = case_when(
      grepl("Invasive heterospecific", contrast) ~ "Invasive heterospecific",
      grepl("Native heterospecific", contrast) ~ "Native heterospecific",
      TRUE ~ contrast
    )
  ) %>%
  mutate(
    contrast = factor(contrast_new, levels = c("Native heterospecific", "Invasive heterospecific"))
  ) %>%
  dplyr::select(-contrast_new)


ggplot(psf_results_clean, aes(x = contrast, y = effect.size, color = Species_status, group = Species_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.0, linewidth = 0.8, 
                position = position_dodge(width = 0.8)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~ Water_treatment, nrow = 2, ncol = 2) +
  scale_color_manual(values = c("Native" = "grey65", "Invasive" = "darkred")) +
  scale_x_discrete(labels = c(
    "Native heterospecific" = "Native heterospecific soil", 
    "Invasive heterospecific" = "Invasive heterospecific soil"
  )) +
  labs(x = "", 
       y = "Plant-soil feedback effect size (Cohen's d \u00B1 95% CI)",
       color = "Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 25, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "top", 
    strip.text = element_text(size = 10, face = "bold", margin = margin(t = 8, b = 8)),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )


## Root Shoot ratio ## Fig.S6
RS <-lmer(sqrt(RootShoot) ~ Species_status * Soil_conditioning * Water_treatment + 
            (1|Species) + (1|Soil_origin_site), data = data_mod)
anova(RS)

pairwise_RS <- emmeans(RS, pairwise ~ Water_treatment, 
                       adjust = "bonferroni")
print(pairwise_RS$contrasts)


## Link phylogenetic distance and traits with plant-soil feedback ##
eff <- read.csv("effect size.csv", header = TRUE)
str(eff)
eff <- eff %>%
  mutate(
    Species_status    = factor(Species_status, levels = c("Native", "Invasive")),
    Water_treatment   = factor(Water_treatment, levels = c("Control", "Drought", 
                                                           "Flood", "Drought + Flood"))
  )

# phylogenetic distance Fig.S10
eff_dis <- lmer(effect.size ~ Distance * Species_status * Water_treatment + (1|Species), 
                data = eff)
test(emtrends(eff_dis,~Species_status * Water_treatment,var = "Distance"))

# root diameter Fig.4A
eff_RD <- lmer(effect.size ~ RD * Species_status * Water_treatment + (1|Species), 
               data = eff)
test(emtrends(eff_RD,~Species_status * Water_treatment,var = "RD"))

# specific root length Fig.4B
eff_SRL <- lmer(effect.size ~ SRL * Species_status * Water_treatment + (1|Species), 
                data = eff)
test(emtrends(eff_SRL,~Species_status * Water_treatment, var = "SRL"))

## trait changes at the end of the conditioning phase
trait <- read.csv("Conditioning traits.csv", header = TRUE)
str(trait)
trait <- trait %>%
  mutate(
    Pot_nr            = as.numeric(Pot_nr),
    Species_status    = factor(Species_status, levels = c("Native", "Invasive")),
    Water_treatment   = factor(Water_treatment, levels = c("Control", "Drought", 
                                                           "Flood", "Drought + Flood")),
    Replicates        = as.numeric(Replicates),
    Soil_origin_site  = as.character(Soil_origin_site)
  )


# RD Fig.S9A 
RD <- lmer(log(RD) ~ Species_status * Water_treatment + (1|Species), 
           # random effect of Soil_origin_site was removed as the variance is 0
           data = trait)
anova(RD)

pairwise_RD <- emmeans(RD, pairwise ~ Water_treatment, 
                       adjust = "bonferroni")
print(pairwise_RD$contrasts)

# SRL Fig.S9B
SRL <- lmer(log(SRL) ~ Species_status * Water_treatment + (1|Species) + (1|Soil_origin_site), 
            data = trait)
anova(SRL)

pairwise_SRL <- emmeans(SRL, pairwise ~ Water_treatment, 
                        adjust = "bonferroni")
print(pairwise_SRL$contrasts)


## PSF on species level ##
species_psf <- eff %>%
  mutate(
    Species = factor(Species, levels = c(
      "Erigeron acris", "Erigeron annuus", 
      "Matricaria chamomilla", "Matricaria discoidea", 
      "Senecio jacobaea", "Senecio inaequidens", 
      "Solidago virgaurea", "Solidago canadensis"
    )),
    Origin_species = as.factor(Origin_species)
  )
str(species_psf)

mod_species <- lmer(effect.size ~ Species * Water_treatment + (1 | Origin_species), data = species_psf)
anova(mod_species)

pairwise <- emmeans(mod_species, pairwise ~ Water_treatment | Species, 
                        adjust = "bonferroni")
print(pairwise$contrasts)

emm_species <- emmeans(mod_species, ~ Species | Water_treatment)
emm_df <- as.data.frame(emm_species)

status_map <- species_psf %>% select(Species, Species_status) %>% distinct()
emm_df <- left_join(emm_df, status_map, by = "Species")

## Fig. S8
ggplot(emm_df, aes(x = Water_treatment, y = emmean, color = Species_status, group = Species_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.0, linewidth = 0.8, 
                position = position_dodge(width = 0.8)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~ Species, ncol = 2) +
  scale_color_manual(values = c("Native" = "grey65", "Invasive" = "darkred")) +
  labs(X = NULL,
       y = "Plant-soil feedback effect size (Cohen's d \u00B1 95% CI)",
       color = "Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 25, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "top", 
    strip.text = element_text(size = 10, face = "bold.italic", margin = margin(t = 8, b = 8)),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )
