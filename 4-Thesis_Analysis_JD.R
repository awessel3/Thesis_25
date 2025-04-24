library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)
library(loo)

library(bayesplot)
library(gridExtra)
library(stringr)
library(rlang)
library(ggridges)
library(forcats)
theme_set(theme_minimal())

# Analysis R Script 
library(brms)

# setwd("~/Desktop/Thesis_25")

trait_species <- read.csv("trait_species.csv")
# WVPT_species <- read.csv("WVPT_species3.3.csv")
traits_full <- read.csv('Traits.csv')

trait_species <- trait_species %>%
  mutate(Full.species.name = case_when(
    Full.species.name == "Clarkia purpurea ssp. quadrivulnera" ~ "Clarkia purpurea",
    TRUE ~ Full.species.name
  ))

trait_species <- rename(trait_species, Species.Name = Full.species.name)

# see_matches <- inner_join(WVPT_species, trait_species, by = 'Species.Name')
# see_matches 
## Prepping data for model 

#WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")
df_flr_final_summary <- read_rds("Data/df_flr_final_filtered_JD.rds")
 
unique(df_flr_final_summary$species)
print(df_flr_final_summary, n=2, width=Inf)

data <- df_flr_final_summary %>% dplyr::select(latitude, longitude, species, spring.temp, spring.precip, winter.temp, winter.precip, elevation, doy, life_history)
data <- na.omit(data)
unique(data$species)

#scaling 
dat.scaled <- data %>% mutate(
  doy_sc = as.vector(scale(doy)),
  latitude_sc = as.vector(scale(latitude)),
  ptemp_sc = as.vector(scale(spring.temp)),
  pprecip_sc = as.vector(scale(spring.precip)),
  elevation_sc = as.vector(scale(elevation))
)
print(dat.scaled,width=Inf,n=2)

# Notes ----
# save traceplots
#
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue", aes=.2) +
    geom_smooth(method = method, color = "red", ...)
  p
}


# Models ----
formula_full <- doy_sc ~ 1 + SPEI * latitude_sc * elevation_sc * life_history + 
  (1 + latitude_sc * elevation_sc * SPEI | species)

formula_full <- doy_sc ~ 1 + spring.temp * spring.precip * latitude_sc * elevation_sc * life_history + 
  (1 + latitude_sc * elevation_sc * spring.temp * spring.precip  | species)


formula_full <- doy_sc ~ 1 + ptemp_sc * latitude_sc * elevation_sc * pprecip_sc *life_history + 
  (1 + ptemp_sc * latitude_sc * elevation_sc * pprecip_sc | species)

formula0 <- doy_sc ~ 1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc + 
  (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

formula1 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + elevation_sc +
  pprecip_sc + life_history + (1 + latitude_sc + ptemp_sc + elevation_sc + pprecip_sc | species)

formula2 <- doy_sc ~ 1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc + life_history + (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

formula3 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + elevation_sc +
  pprecip_sc + life_history + (1 + ptemp_sc + latitude_sc +  elevation_sc + pprecip_sc | species)

formula4 <- doy_sc ~ 1 + (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

formula5 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + ptemp_sc * elevation_sc + pprecip_sc + (1 | species)

## fit0 ----

fit0 <- brm(
  formula = formula0,
  data = dat.scaled,
  family = gaussian(),  # Assuming DOY is approximately normally distributed
  # control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  init = "0"
)

summary(fit0)
get_variables(fit0)
#pp_check_fit#
pp_check(fit0, plotfun = "dens_overlay")

bayesplot::ppc_scatter_avg(y = data$doy_sc, yrep = posterior_predict(fit0))
plot(fit)
pairs(fit, np = nuts_params(fit))

#original model save
saveRDS(fit0, file = "../iNaturalist_SHARED/WVPT analysis/saved_models/fit0.RDS")
fit0 <- readRDS("../iNaturalist_SHARED/WVPT analysis/saved_models/fit0.RDS")  

pdf(file="../iNaturalist_SHARED/WVPT analysis/saved_models/fit0.traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit0, regex_pars = "^b_") + xlab("Post-warmup iteration"))
# print(mcmc_trace(fit0, regex_pars = ".*") + xlab("Post-warmup iteration"))
# print(mcmc_trace(fit0, regex_pars = "alpha") + xlab("Post-warmup iteration"))
dev.off()


# Get a table of parameter values
get_variables(fit0)
temp <- as.data.frame(signif(fixef(fit0), 3)) %>%
  rownames_to_column(var = "Parameter")
pdf(file = paste("../iNaturalist_SHARED/WVPT analysis/saved_models/fit0.TABLE.pdf", sep=''), height=15, width=15)
grid.table(temp)
dev.off()



#pass 1 save 
#saveRDS(fit, file = "Data/fit1.RDS")
fit1 <- readRDS("Data/fit1.RDS")  

# pass 2 save
#saveRDS(fit, file = "Data/fit2.RDS")
fit2 <- readRDS("Data/fit2.RDS")  

#pass 3 save 
#saveRDS(fit, file = "Data/fit3.RDS")
fit3 <- readRDS("Data/fit3.RDS") 

#pass4 save 
#saveRDS(fit, file = "Data/fit4.RDS")
fit4 <- readRDS("Data/fit4.RDS") 

#pass5 save 
#saveRDS(fit, file = "Data/fit5.RDS")
fit5 <- readRDS("Data/fit5.RDS") 

loo1 <- loo(fit0, fit_full)
loo1


loo1 <- loo(fit1, fit2)
loo1

 
as_draws_df(fit) %>% head(3)

lp_draws <- linpred_draws(fit, newdata = original_data)


# Plot model output ----

#fit <- fit3
spp <- unique(data$species)
spp

## Slopes ----

fixef_draws <- fit0 %>%
  spread_draws(`b_.*`, regex = TRUE)  # this uses a regex to match all fixed effects

# overall slopes
fixef_long <- fixef_draws %>%
  pivot_longer(
    cols = starts_with("b_"),
    names_to = "term",
    values_to = "estimate"
  ) %>%
  mutate(
    # Remove 'b_' prefix and replace colons with multiplication symbol for clarity
    term_clean = term %>%
      str_remove("^b_") %>%
      str_replace_all(":", " × ") %>%
      str_replace_all("_sc", "")  # remove "_sc" if you want to simplify variable names
  )
head(fixef_long)

ggplot(fixef_long, aes(x = estimate, fill = term_clean)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ term_clean, scales = "free") +
  geom_vline(xintercept = 0, linetype='dashed')+
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Distributions of Fixed Effects",
    x = "Estimate",
    y = "Density"
  ) +
  theme(legend.position = "none")


# Extract random slopes for each species

# Get the random slope term names
slope_terms <- get_variables(fit0) %>%
  str_subset("^r_species\\[.*?,(.*?)\\]$") %>%
  str_match("^r_species\\[.*?,(.*?)\\]$") %>%
  .[, 2] %>%
  unique()

# Build expressions like r_species[species, <term>]
random_exprs <- paste0("r_species[species,", slope_terms, "]") %>%
  parse_exprs()

# Build fixed effect names
fixed_exprs <- paste0("b_", slope_terms) %>% syms()


# Fixed effect draws
fixed_draws <- fit0 %>%
  spread_draws(!!!fixed_exprs)

# Random effect draws for species
random_draws <- fit0 %>%
  spread_draws(!!!random_exprs)

# combine and calculate species-specific slopes

# extract raw term names for all species
all_r_vars <- get_variables(fit0) %>%
  str_subset("^r_species\\[.*?,.*\\]$")
head(all_r_vars)

random_draws <- fit0 %>%
  spread_draws(!!!syms(all_r_vars))

random_long <- random_draws %>%
  pivot_longer(
    cols = any_of(all_r_vars),  # now ensures valid columns
    names_to = "term_raw",
    values_to = "rand"
  ) %>%
  mutate(
    term = str_match(term_raw, "\\[.*?,(.*?)\\]")[, 2],
    species = str_match(term_raw, "\\[(.*?),")[, 2],
    fixed_name = paste0("b_", term),
    rand = as.numeric(rand)  # <---- KEY FIX
  )

fixed_long <- fixed_draws %>%
  pivot_longer(
    cols = all_of(paste0("b_", slope_terms)),
    names_to = "fixed_name",
    values_to = "fix"
  ) %>%
  mutate(fix = as.numeric(fix))  # ensure numeric!

species_slopes <- random_long %>%
  left_join(fixed_long, by = c(".chain", ".iteration", ".draw", "fixed_name")) %>%
  mutate(
    species_slope = fix + rand,
    term_clean = term %>%
      str_replace_all(":", " × ") %>%
      str_replace_all("_sc", "")
  ) %>%
  filter(term != "Intercept")   # skip intercepts
  
head(species_slopes)

species_slopes %>%
  ggplot(aes(x = species_slope, fill = species)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ term_clean, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Posterior Distributions of Species-Specific Slopes",
    x = "Slope Estimate",
    y = "Density"
  )

# 1. Calculate per-species × term probabilities of direction
probs_df <- species_slopes %>%
  group_by(species, term_clean) %>%
  summarise(
    prob_gt0 = mean(species_slope > 0),
    prob_diff0 = if_else(prob_gt0 > 0.5, prob_gt0, 1 - prob_gt0),
    prob_label = scales::percent(prob_diff0, accuracy = 1),
    .groups = "drop"
  )

# 2. Join back to full draws to retain point_interval plotting
slopes_for_plot <- species_slopes %>%
  left_join(probs_df, by = c("species", "term_clean"))

# 3. Order species within each term by decreasing prob_diff0
slopes_for_plot <- slopes_for_plot %>%
  group_by(term_clean) %>%
  mutate(species = fct_reorder(species, prob_diff0)) %>%
  ungroup()

# 4. Plot using full posterior samples
ggplot(slopes_for_plot, aes(x = species_slope, y = species)) +
  stat_pointinterval(.width = 0.90, point_size = 2.5, fatten_point = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(
    data = probs_df,
    aes(x = 0.05, y = species, label = prob_label),
    inherit.aes = FALSE,
    size = 3, hjust = 0
  ) +
  facet_wrap(~ term_clean, scales = "free_x") +
  labs(
    title = "Species-specific Slopes with 95% Credible Intervals",
    subtitle = "Posterior P(slope ≠ 0) shown as % next to each species",
    x = "Slope Estimate",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

## Ridge Plot!! ----
# library(ggridges)

# reorder species for each term by probability
slopes_ridge <- slopes_for_plot %>%
  group_by(term_clean) %>%
  mutate(species = fct_reorder(species, prob_diff0)) %>%
  ungroup()

# ggplot(slopes_ridge, aes(x = species_slope, y = species, fill = prob_diff0)) +
#   geom_density_ridges(
#     scale = 1.1, alpha = 0.8,
#     quantile_lines = TRUE, quantiles = 0.5,
#     vline_color = "gray50"
#   ) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
#   facet_wrap(~ term_clean, scales = "free_x") +
#   scale_fill_viridis_c(option = "plasma", name = "P(≠0)") +
#   labs(
#     title = "Species-Specific Slope Distributions",
#     subtitle = "Colored by P(slope ≠ 0), with medians and density",
#     x = "Slope Estimate", y = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "bottom",
#     panel.spacing = unit(1, "lines")
#   )


# Binned version

slopes_binned <- species_slopes %>%
  group_by(term_clean, species) %>%
  summarise(
    prob_diff0 = mean(species_slope > 0) %>% {pmax(., 1 - .)},
    .groups = "drop"
  ) %>%
  mutate(prob_bin = case_when(
    prob_diff0 > 0.95 ~ "P > 0.95",
    prob_diff0 > 0.9 ~ "0.90 < P ≤ 0.95",
    prob_diff0 > 0.8 ~ "0.80 < P ≤ 0.90",
    TRUE ~ "P ≤ 0.80"
  )) %>%
  left_join(species_slopes, by = c("term_clean", "species")) %>%
  group_by(term_clean) %>%
  mutate(species_ordered = fct_reorder(species, prob_diff0)) %>%
  ungroup()


ggplot(slopes_binned, aes(x = species_slope, y = species_ordered, fill = prob_bin)) +
  geom_density_ridges(scale = 1.1, rel_min_height = 0.01, color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ term_clean, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "P > 0.95" = "#2166AC",
      "0.90 < P ≤ 0.95" = "#67A9CF",
      "0.80 < P ≤ 0.90" = "#D1E5F0",
      "P ≤ 0.80" = "#FDDBC7"
    ),
    name = "Probability ≠ 0"
  ) +
  labs(
    title = "Species-Specific Posterior Slopes",
    x = "Posterior Slope",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )
# 


# fit_full ----

# formula_full <- doy_sc ~ 1 + ptemp_sc * latitude_sc * elevation_sc * pprecip_sc *life_history + 
#  (1 + ptemp_sc * latitude_sc * elevation_sc * pprecip_sc | species)

fit_full <- brm(
  formula = formula_full,
  data = dat.scaled,
  family = gaussian(),  # Assuming DOY is approximately normally distributed
  # control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  init = "0"
)

summary(fit_full)
# saveRDS(fit_full, file = "Data/fit_full.RDS")
fit_full <- readRDS("Data/fit_full.RDS") 


fixef_draws <- fit_full %>%
  spread_draws(`b_.*`, regex = TRUE)  # this uses a regex to match all fixed effects

# overall slopes
fixef_long <- fixef_draws %>%
  pivot_longer(
    cols = starts_with("b_"),
    names_to = "term",
    values_to = "estimate"
  ) %>%
  mutate(
    # Remove 'b_' prefix and replace colons with multiplication symbol for clarity
    term_clean = term %>%
      str_remove("^b_") %>%
      str_replace_all(":", " × ") %>%
      str_replace_all("_sc", "")  # remove "_sc" if you want to simplify variable names
  )
head(fixef_long,2)

ggplot(fixef_long, aes(x = estimate, fill = term_clean)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ term_clean, scales = "free") +
  geom_vline(xintercept = 0, linetype='dashed')+
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Distributions of Fixed Effects",
    x = "Estimate",
    y = "Density"
  ) +
  theme(legend.position = "none")
ggsave("Figures/fit_full_main.pdf", width=12, height=12)

# do a ridge plot
fixef_probs <- fixef_long %>%
  group_by(term_clean) %>%
  summarise(
    prob_diff0 = mean(estimate > 0) %>% {pmax(., 1 - .)},
    .groups = "drop"
  ) %>%
  mutate(
    prob_label = paste0("P ≠ 0: ", round(prob_diff0, 2))
  )

# Step 2: Join with original data
fixef_ridges <- fixef_long %>%
  left_join(fixef_probs, by = "term_clean") %>%
  mutate(term_ordered = fct_reorder(term_clean, prob_diff0))

# Step 3: Plot
ggplot(fixef_ridges, aes(x = estimate, y = term_ordered, fill = term_ordered)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_text(
    data = distinct(fixef_ridges, term_ordered, prob_label),
    aes(x = -Inf, y = term_ordered, label = prob_label),
    hjust = -0.05,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  labs(
    title = "Posterior Distributions of Fixed Effects",
    x = "Posterior Estimate",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12)
  ) +
  coord_cartesian(xlim = c(min(fixef_ridges$estimate), max(fixef_ridges$estimate) + 0.5))
ggsave("Figures/fit_full_params.pdf", width=10, height=10)



# Step 1: Compute posterior probabilities
fixef_probs <- fixef_long %>%
  group_by(term_clean) %>%
  summarize(
    prob_diff0 = abs(mean(estimate > 0) - 0.5) * 2,
    .groups = "drop"
  ) %>%
  mutate(
    prob_category = case_when(
      prob_diff0 > 0.95 ~ "P > 0.95",
      prob_diff0 > 0.90 ~ "0.90 < P ≤ 0.95",
      prob_diff0 > 0.80 ~ "0.80 < P ≤ 0.90",
      TRUE ~ "≤ 0.80"
    ),
    prob_category = factor(prob_category, levels = c("P > 0.95", "0.90 < P ≤ 0.95", "0.80 < P ≤ 0.90", "≤ 0.80"))
  )

# Step 2: Join to full draws
fixef_annotated <- fixef_long %>%
  left_join(fixef_probs, by = "term_clean") %>%
  mutate(
    # Determine interaction order
    interaction_level = case_when(
      str_count(term_clean, " × ") == 0 ~ "1 Main Effect",
      str_count(term_clean, " × ") == 1 ~ "2 Two-way Interaction",
      str_count(term_clean, " × ") == 2 ~ "3 Three-way Interaction",
      str_count(term_clean, " × ") == 3 ~ "4 Four-way Interaction"
    ),
    interaction_level = factor(interaction_level, levels = c("1 Main Effect", "2 Two-way Interaction", "3 Three-way Interaction", "4 Four-way Interaction"))
  )

# Step 3: Order terms within their interaction class by prob_diff0
term_order <- fixef_annotated %>%
  distinct(term_clean, interaction_level, prob_diff0) %>%
  arrange(interaction_level, desc(prob_diff0)) %>%
  mutate(term_clean = factor(term_clean, levels = rev(term_clean)))

fixef_annotated <- fixef_annotated %>%
  mutate(term_clean = factor(term_clean, levels = levels(term_order$term_clean)))

# Step 4: Plot
ggplot(fixef_annotated, aes(x = estimate, y = term_clean, fill = prob_category)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ interaction_level, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_fill_manual(
    name = "P(β ≠ 0)",
    values = c(
      "P > 0.95" = "#1b7837",
      "0.90 < P ≤ 0.95" = "#a6dba0",
      "0.80 < P ≤ 0.90" = "#fdae61",
      "≤ 0.80" = "#d73027"
    )
  ) +
  labs(
    title = "Posterior Distributions of Fixed Effects",
    x = "Posterior Estimate",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0)
  )



# Extract random slopes for each species

# Get the random slope term names
slope_terms <- get_variables(fit_full) %>%
  str_subset("^r_species\\[.*?,(.*?)\\]$") %>%
  str_match("^r_species\\[.*?,(.*?)\\]$") %>%
  .[, 2] %>%
  unique()

# Build expressions like r_species[species, <term>]
random_exprs <- paste0("r_species[species,", slope_terms, "]") %>%
  parse_exprs()

# Build fixed effect names
fixed_exprs <- paste0("b_", slope_terms) %>% syms()


# Fixed effect draws
fixed_draws <- fit_full %>%
  spread_draws(!!!fixed_exprs)

# Random effect draws for species
random_draws <- fit_full %>%
  spread_draws(!!!random_exprs)

# combine and calculate species-specific slopes

# extract raw term names for all species
all_r_vars <- get_variables(fit_full) %>%
  str_subset("^r_species\\[.*?,.*\\]$")
head(all_r_vars)

random_draws <- fit_full %>%
  spread_draws(!!!syms(all_r_vars))

random_long <- random_draws %>%
  pivot_longer(
    cols = any_of(all_r_vars),  # now ensures valid columns
    names_to = "term_raw",
    values_to = "rand"
  ) %>%
  mutate(
    term = str_match(term_raw, "\\[.*?,(.*?)\\]")[, 2],
    species = str_match(term_raw, "\\[(.*?),")[, 2],
    fixed_name = paste0("b_", term),
    rand = as.numeric(rand)  # <---- KEY FIX
  )

fixed_long <- fixed_draws %>%
  pivot_longer(
    cols = all_of(paste0("b_", slope_terms)),
    names_to = "fixed_name",
    values_to = "fix"
  ) %>%
  mutate(fix = as.numeric(fix))  # ensure numeric!

species_slopes <- random_long %>%
  left_join(fixed_long, by = c(".chain", ".iteration", ".draw", "fixed_name")) %>%
  mutate(
    species_slope = fix + rand,
    term_clean = term %>%
      str_replace_all(":", " × ") %>%
      str_replace_all("_sc", "")
  ) %>%
  filter(term != "Intercept")   # skip intercepts

head(species_slopes)

species_slopes %>%
  ggplot(aes(x = species_slope, fill = species)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ term_clean, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Posterior Distributions of Species-Specific Slopes",
    x = "Slope Estimate",
    y = "Density"
  )

# 1. Calculate per-species × term probabilities of direction
probs_df <- species_slopes %>%
  group_by(species, term_clean) %>%
  summarise(
    prob_gt0 = mean(species_slope > 0),
    prob_diff0 = if_else(prob_gt0 > 0.5, prob_gt0, 1 - prob_gt0),
    prob_label = scales::percent(prob_diff0, accuracy = 1),
    .groups = "drop"
  )

# 2. Join back to full draws to retain point_interval plotting
slopes_for_plot <- species_slopes %>%
  left_join(probs_df, by = c("species", "term_clean"))

# 3. Order species within each term by decreasing prob_diff0
slopes_for_plot <- slopes_for_plot %>%
  group_by(term_clean) %>%
  mutate(species = fct_reorder(species, prob_diff0)) %>%
  ungroup()

# 4. Plot using full posterior samples
ggplot(slopes_for_plot, aes(x = species_slope, y = species)) +
  stat_pointinterval(.width = 0.90, point_size = 2.5, fatten_point = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(
    data = probs_df,
    aes(x = 0.05, y = species, label = prob_label),
    inherit.aes = FALSE,
    size = 3, hjust = 0
  ) +
  facet_wrap(~ term_clean, scales = "free_x") +
  labs(
    title = "Species-specific Slopes with 95% Credible Intervals",
    subtitle = "Posterior P(slope ≠ 0) shown as % next to each species",
    x = "Slope Estimate",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

## Ridge Plot!! ----
# library(ggridges)

# reorder species for each term by probability
slopes_ridge <- slopes_for_plot %>%
  group_by(term_clean) %>%
  mutate(species = fct_reorder(species, prob_diff0)) %>%
  ungroup()

# ggplot(slopes_ridge, aes(x = species_slope, y = species, fill = prob_diff0)) +
#   geom_density_ridges(
#     scale = 1.1, alpha = 0.8,
#     quantile_lines = TRUE, quantiles = 0.5,
#     vline_color = "gray50"
#   ) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
#   facet_wrap(~ term_clean, scales = "free_x") +
#   scale_fill_viridis_c(option = "plasma", name = "P(≠0)") +
#   labs(
#     title = "Species-Specific Slope Distributions",
#     subtitle = "Colored by P(slope ≠ 0), with medians and density",
#     x = "Slope Estimate", y = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "bottom",
#     panel.spacing = unit(1, "lines")
#   )


# Binned version

slopes_binned <- species_slopes %>%
  group_by(term_clean, species) %>%
  summarise(
    prob_diff0 = mean(species_slope > 0) %>% {pmax(., 1 - .)},
    .groups = "drop"
  ) %>%
  mutate(prob_bin = case_when(
    prob_diff0 > 0.95 ~ "P > 0.95",
    prob_diff0 > 0.9 ~ "0.90 < P ≤ 0.95",
    prob_diff0 > 0.8 ~ "0.80 < P ≤ 0.90",
    TRUE ~ "P ≤ 0.80"
  )) %>%
  left_join(species_slopes, by = c("term_clean", "species")) %>%
  group_by(term_clean) %>%
  mutate(species_ordered = fct_reorder(species, prob_diff0)) %>%
  ungroup()


ggplot(slopes_binned, aes(x = species_slope, y = species_ordered, fill = prob_bin)) +
  geom_density_ridges(scale = 1.1, rel_min_height = 0.01, color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ term_clean, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "P > 0.95" = "#2166AC",
      "0.90 < P ≤ 0.95" = "#67A9CF",
      "0.80 < P ≤ 0.90" = "#D1E5F0",
      "P ≤ 0.80" = "#FDDBC7"
    ),
    name = "Probability ≠ 0"
  ) +
  labs(
    title = "Species-Specific Posterior Slopes",
    x = "Posterior Slope",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )
# 


## Visualize predictions ----


# Step 1: Define focal variables and scaling parameters
focal_vars <- c("ptemp", "pprecip", "elevation", "latitude")

# Create a helper to unscale a variable
unscale <- function(scaled, original) {
  mean_val <- mean(original, na.rm = TRUE)
  sd_val   <- sd(original, na.rm = TRUE)
  return(scaled * sd_val + mean_val)
}

# Step 2: Create newdata for each focal variable on unscaled scale
newdata_list <- map(focal_vars, function(var) {
  var_sc <- paste0(var, "_sc")
  var_seq <- seq(
    from = quantile(data[[var]], 0.05, na.rm = TRUE),
    to   = quantile(data[[var]], 0.95, na.rm = TRUE),
    length.out = 100
  )
  
  # Scale it using original data mean/sd
  var_seq_sc <- (var_seq - mean(data[[var]], na.rm = TRUE)) / sd(data[[var]], na.rm = TRUE)
  
  # Create mean values for other predictors
  others <- setdiff(focal_vars, var)
  others_means <- map_dfc(others, ~ {
    val <- mean(data[[.x]], na.rm = TRUE)
    val_sc <- (val - mean(data[[.x]], na.rm = TRUE)) / sd(data[[.x]], na.rm = TRUE)
    tibble(!!paste0(.x, "_sc") := val_sc)
  })
  
  tibble(!!var_sc := var_seq_sc) %>%
    bind_cols(others_means) %>%
    mutate(!!var := var_seq)  # keep unscaled x variable for plotting
}) %>% set_names(focal_vars)

# Step 3: Get linpred_draws + convert to unscaled DOY
plots <- imap(newdata_list, function(new_df, var_name) {
  draws <- linpred_draws(fit0, newdata = new_df, ndraws = 1000, allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(x = new_df[[var_name]])
  
  # Unscale DOY (assuming DOY was scaled in model input)
  doy_mean <- mean(data$DOY, na.rm = TRUE)
  doy_sd   <- sd(data$DOY, na.rm = TRUE)
  
  draws <- draws %>%
    mutate(DOY_pred = .linpred * doy_sd + doy_mean)
  
  # Raw data layer
  raw_df <- data %>%
    select(DOY, all_of(var_name)) %>%
    drop_na()
  
  ggplot(draws, aes(x = x, y = DOY_pred)) +
    geom_point(data = raw_df, aes_string(x = var_name, y = "DOY"), color = "gray70", alpha = 0.5, inherit.aes = FALSE) +
    stat_lineribbon(.width = c(0.5, 0.8, 0.95), fill = "#3182bd", alpha = 0.6) +
    labs(
      x = var_name,
      y = "Day of Year",
      title = paste("Effect of", var_name, "on DOY")
    ) +
    theme_minimal(base_size = 14)
})

# Step 4: Combine 4 plots into 2x2 panel
(plots$ptemp | plots$pprecip) /
  (plots$latitude | plots$elevation) +
  plot_annotation(title = "Posterior Estimates with Raw Data")



## Predictions ----
## vs latitude
# Create a new dataset to predict over
data.predict <- crossing(
  elevation_sc = mean(data$elevation_sc), # predictions at mean elevation
  ptemp_sc = mean(data$ptemp_sc), # temperature range
  pprecip_sc = mean(data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = seq(min(data$latitude_sc), max(data$latitude_sc), length.out = 5),
  species = spp)

# make predictions using the fitted model
fitted.pred <-  linpred_draws(object = fit0, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

# Unscale the predicted DOY
fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * sd(data$doy)) + mean(data$doy)
         ,latitude = (latitude_sc * sd(data$latitude)) + mean(data$latitude)
         ,elevation = (elevation_sc * sd(data$elevation)) + mean(data$elevation)
         ,ptemp = (ptemp_sc * sd(data$spring.temp)) + mean(data$spring.temp)
         ,pprecip = (pprecip_sc * sd(data$spring.precip)) + mean(data$spring.precip)
         )  
head(fitted.pred)

# Plot
ggplot(fitted.pred, aes(x = latitude , y = DOY_pred, color=species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="Latitude") +
  scale_fill_brewer(palette = "Greys", guide = "none") 
ggsave("Figures/DOY_latitude.pdf", width=9, height=6)

ggplot(fitted.pred, aes(x = latitude , y = DOY_pred))+
  geom_jitter(alpha=.2)+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="Latitude") +
  facet_wrap(~ species, scales = "free")
  # scale_fill_brewer(palette = "Greys", guide = "none") 
ggsave("Figures/DOY_latitude_facet.pdf", width=12, height=8)


#Temp vs DOY
data.predict <- crossing(
  elevation_sc = mean(data$elevation_sc), # predictions at mean elevation
  ptemp_sc = seq(min(data$ptemp_sc), max(data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(data$latitude_sc, na.rm = TRUE),
  species = spp
  )


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring.temp = (ptemp_sc * ptemp_scale) + ptemp_center)



# Plot
# fit#_TempDOY_plot
ggplot(fitted.pred, aes(x = spring.temp , y = DOY_pred, color = species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="preceding temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()


#Temp|latitude vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = mean(data$elevation_sc, na.rm = TRUE), # predictions at mean elevation
  ptemp_sc = seq(min(data$ptemp_sc), max(data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = seq(min(data$latitude_sc), max(data$latitude_sc), length.out = 5),
  species = spp
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring.temp = (ptemp_sc * ptemp_scale) + ptemp_center,
         latitude = (latitude_sc * latitude_scale) + latitude_center)

# Plot
# fit#_TempLatDOY_plot
ggplot(fitted.pred, aes(x = spring.temp, y = DOY_pred, color = factor(round(latitude, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

# SS_fit#_TempLatDOY_plot
ggplot(fitted.pred, aes(x = spring.temp, y = DOY_pred, color = factor(round(latitude, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()

#Temp|Elevation vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = seq(min(data$elevation_sc), max(data$elevation_sc), length.out = 5), # predictions at mean elevation
  ptemp_sc = seq(min(data$ptemp_sc), max(data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(data$latitude_sc, na.rm = TRUE),
  species = spp
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring.temp = (ptemp_sc * ptemp_scale) + ptemp_center,
         elevation = (elevation_sc * elevation_scale) + elevation_center)

# Plot
#fit#_TempElevDOY_plot
ggplot(fitted.pred, aes(x = spring.temp, y = DOY_pred, color = factor(round(elevation, 0)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

# SS_fit#_TempElevDOY_plot
ggplot(fitted.pred, aes(x = spring.temp, y = DOY_pred, color = factor(round(elevation, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()



# Creating phenological sensitivity plot ----

## Random Slopes just for Temp 
# view variables that can be extracted 
get_variables(fit)

# extracr species specific slopes for temperature
sp_slopes_posterior_temp <- fit %>%
  spread_draws(r_species[species, ptemp_sc ])  
head(sp_slopes_posterior_temp)

unique(sp_slopes_posterior_temp$ptemp_sc)

ggplot(sp_slopes_posterior_temp, aes(x = species, y = r_species)) + 
  stat_summary(geom = "pointrange", fun.data = mean_hdi) +
  labs(title = "Species-Specific Slopes", x = "Species", y = "Slope") +
  theme_minimal()

species_temp_ps <- fit %>%
  spread_draws(r_species[species, ptemp_sc]) %>%
  group_by(species) %>%
  summarise(
    mean_sensitivity = mean(r_species),  
    lower_95 = quantile(r_species, 0.025),
    upper_95 = quantile(r_species, 0.975)
  ) 

# unscaling estimated parameters 
species_temp_ps <- species_temp_ps %>%
  mutate(sensitivity_unscaled = mean_sensitivity * (doy_scale / ptemp_scale),
         lower_95_unscaled = lower_95 * (doy_scale / ptemp_scale),
         upper_95_unscaled = upper_95 * (doy_scale / ptemp_scale))


#fixing species names containing . between 
species_temp_ps$species <- gsub("\\.", " ", species_temp_ps$species)

# Extract mean flowering
mean_doy <- data %>%
  group_by(species) %>%
  summarise(mean_doy = mean(((doy_sc * doy_scale) + doy_center), na.rm = TRUE))

# Merge data sets
temp_ps_plot_dat <- left_join(species_temp_ps, mean_doy, by = 'species')
temp_ps_plot_dat <- left_join(temp_ps_plot_dat, life_hist, by = "species")

# fit0_temp_PS_plot
ggplot(temp_ps_plot_dat, aes(x = mean_doy, y = sensitivity_unscaled, color = species)) +
  geom_point(size = 3, aes(color = species)) + 
  stat_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.75) + 
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower_95_unscaled, ymax = upper_95_unscaled), linewidth = 1, size = 0.7)


# removing Bidens frondosa  - why so high? 
temp_ps_woBF <- temp_ps_plot_dat %>% 
  filter(species != "Bidens frondosa")

ggplot(temp_ps_woBF, aes(x = mean_doy, y = sensitivity_unscaled, color = life_history)) +
  geom_point(size = 3, aes(color = life_history)) + 
  stat_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.75) + 
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower_95_unscaled, ymax = upper_95_unscaled), linewidth = 1, size = 0.7)

# life history comparison plot 

ps_life_hist <- temp_ps_woBF %>% 
  group_by(life_history) %>% 
  summarise(mean_sensitivity = mean(mean_sensitivity),
            sensitivity_unscaled = mean(sensitivity_unscaled))

ps_life_hist

ggplot(ps_life_hist, aes(x = life_history, y = sensitivity_unscaled, color = life_history)) +
  geom_point(size = 3, aes(color = life_history)) + 
  geom_hline(yintercept = 0) 
