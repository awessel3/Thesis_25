library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)
library(loo)


# Analysis R Script 
library(brms)
library(ggridges)

setwd("~/Desktop/Thesis_25")

df_flr_final_summary <- read_rds("Data/df_flr_final_summary.rds")
df_flr_final_filtered <- read_rds("Data/df_flr_final_filtered.rds")
trait_species <- read.csv("trait_species.csv")
traits_full <- read.csv('Traits.csv')
traits_full$SpName <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", traits_full$SpName)

traits_full <- rename(traits_full, species = SpName)


life_hist <- read.csv("Data/Flowering_WVPT_life_history.csv")
spp <- unique(df_flr_final_filtered$species)

df_traits_flr <- traits_full %>%
  filter(species %in% spp) 
df_traits_flr_final <- left_join(df_traits_flr, life_hist, by = "species")

length(unique(df_traits_flr_final$species))

trait_hist_count <- df_traits_flr_final %>% 
  group_by(life_history) %>% 
  summarise(count = n())
trait_hist_count <- count(df_traits_flr_final, c("annual"))
trait_hist_count

## Prepping data for model 

unique(df_flr_final_summary$species)

#Excluding species - either for low number of observations or skewing data 


data <- df_flr_final_filtered %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                               preceding_precip, elevation, doy, life_history)
data <- na.omit(data)
unique(data$species)

#scaling 
doy_num <- as.numeric(data$doy)
doy_center <- mean(doy_num)
doy_scale <- sd(doy_num)
data$doy_sc <- (doy_num - doy_center) / doy_scale

latitude_num <- as.numeric(data$latitude)
latitude_center <- mean(latitude_num)
latitude_scale <- sd(latitude_num)
data$latitude_sc <- (latitude_num - latitude_center) / latitude_scale

ptemp_num <- as.numeric(data$preceding_temp)
ptemp_center <- mean(ptemp_num)
ptemp_scale <- sd(ptemp_num)
data$ptemp_sc <- (ptemp_num - ptemp_center) / ptemp_scale

pprecip_num <- as.numeric(data$preceding_precip)
pprecip_center <- mean(pprecip_num)
pprecip_scale <- sd(pprecip_num)
data$pprecip_sc <- (pprecip_num - pprecip_center) / pprecip_scale     

elevation_num <- as.numeric(data$elevation)
elevation_center <- mean(elevation_num)
elevation_scale <- sd(elevation_num)
data$elevation_sc <- (elevation_num - elevation_center) / elevation_scale

# fit0
formula <- doy_sc ~ 1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc + 
  (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

#altering model to leave on out to recognize necessary complexity 

# pass 1: fit1
#formula1 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + elevation_sc +
 # pprecip_sc + life_history + (1 + latitude_sc + ptemp_sc + elevation_sc + pprecip_sc | species)

# pass 2: 
#fit2
formula2 <- doy_sc ~ 1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc + life_history + (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

# pass 3: 
#fit3 
formula3 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + elevation_sc +
  pprecip_sc + life_history + (1 + ptemp_sc + latitude_sc +  elevation_sc + pprecip_sc | species)

#pass 4:
formula4 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc + life_history + (1 | species)

#pass 5: 
formula_full <- doy_sc ~ 1 + ptemp_sc * latitude_sc * elevation_sc * pprecip_sc * life_history +
  (1 + latitude_sc * ptemp_sc * elevation_sc * pprecip_sc | species)

#pass 6: 



fit <- brm(
  formula = formula_full,
  data = data,
  family = gaussian(),  
  #control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  #init = "0"
)

summary(fit)
#pp_check_fit#
pp_check(fit, plotfun = "dens_overlay")

bayesplot::ppc_scatter_avg(y = data$doy_sc, yrep = posterior_predict(fit))
plot(fit)
pairs(fit, np = nuts_params(fit))

#original model save
#saveRDS(fit, file = "Data/fit0.RDS")
fit0 <- readRDS("Data/fit0.RDS")  

#pass 1 save 
#saveRDS(fit, file = "Data/fit1.RDS")
fit1 <- readRDS("Data/fit1.RDS")  

# pass 2 save
#saveRDS(fit, file = "Data/fit2.RDS")
fit2 <- readRDS("Data/fit2.RDS")  

#pass 3 save 
saveRDS(fit, file = "Data/fit3.RDS")
fit3 <- readRDS("Data/fit3.RDS") 


saveRDS(fit, file = "Data/fit_full.RDS")
fit_full <- readRDS("Data/fit_full.RDS") 

loo1 <- loo(fit_full)
loo1


#prep for plotting 
fit <- fit_full
original_data <- fit$data 
spp <- unique(original_data$species)
spp
life_history <- unique(original_data$life_history)
life_history


## Initial Plot Creation ----
# Create a new dataset to predict over
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc), # predictions at mean elevation
  ptemp_sc = mean(original_data$ptemp_sc), # temperature range
  pprecip_sc = mean(original_data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = seq(min(original_data$latitude_sc), max(original_data$latitude_sc), length.out = 5),
  species = spp, 
  life_history = life_history)


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object = fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

# Unscale the predicted DOY
fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         latitude = (latitude_sc * latitude_scale) + latitude_center)  


# Plot
#fit#_LatDOY_plot
ggplot(fitted.pred, aes(x = latitude , y = DOY_pred, color=species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="Latitude") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()
#ggsave("Figures/DOY_latitude.pdf", width=5, height=4)


#Temp vs DOY
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc), # predictions at mean elevation
  ptemp_sc = seq(min(original_data$ptemp_sc), max(original_data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(original_data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
  )


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         preceding_temp = (ptemp_sc * ptemp_scale) + ptemp_center)



# Plot
# fit#_TempDOY_plot
ggplot(fitted.pred, aes(x = preceding_temp , y = DOY_pred, color = species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="preceding temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()


#Temp|latitude vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc, na.rm = TRUE), # predictions at mean elevation
  ptemp_sc = seq(min(original_data$ptemp_sc), max(original_data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(original_data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = seq(min(original_data$latitude_sc), max(original_data$latitude_sc), length.out = 5),
  species = spp,
  life_history = life_history
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         preceding_temp = (ptemp_sc * ptemp_scale) + ptemp_center,
         latitude = (latitude_sc * latitude_scale) + latitude_center)

# Plot
# fit#_TempLatDOY_plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(round(latitude, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

# SS_fit#_TempLatDOY_plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(round(latitude, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()

#Temp|Elevation vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = seq(min(original_data$elevation_sc), max(original_data$elevation_sc), length.out = 5), # predictions at mean elevation
  ptemp_sc = seq(min(original_data$ptemp_sc), max(original_data$ptemp_sc), length.out = 100), # temperature range
  pprecip_sc = mean(original_data$pprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         preceding_temp = (ptemp_sc * ptemp_scale) + ptemp_center,
         elevation = (elevation_sc * elevation_scale) + elevation_center)

# Plot
#fit#_TempElevDOY_plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(round(elevation, 0)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

# SS_fit#_TempElevDOY_plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(round(elevation, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()



## Creating phenological sensitivity plot ----

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
mean_doy <- original_data %>%
  group_by(species) %>%
  summarise(mean_doy = mean(((doy_sc * doy_scale) + doy_center), na.rm = TRUE))

# Merge data sets
temp_ps_plot_dat <- left_join(species_temp_ps, mean_doy, by = 'species')
temp_ps_plot_dat <- left_join(temp_ps_plot_dat, life_hist, by = "species")
temp_ps_plot_dat


# fit#_temp_PS_plot
ggplot(temp_ps_plot_dat, aes(x = mean_doy, y = sensitivity_unscaled, color = life_history)) +
  geom_point(size = 3, aes(color = life_history)) + 
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

# Extract posterior draws for the life history (perennial effect)
life_history_post <- fit %>%
  spread_draws(b_life_historyperennial)

head(life_history_post)


ps_life_hist <- temp_ps_plot_dat %>%
  group_by(life_history) %>% 
  summarise(
    mean_sensitivity = mean(sensitivity_unscaled, na.rm = TRUE),  
    lower_95 = quantile(sensitivity_unscaled, 0.025, na.rm = TRUE),  
    upper_95 = quantile(sensitivity_unscaled, 0.975, na.rm = TRUE))

ps_life_hist

ggplot(ps_life_hist, aes(x = life_history, y = mean_sensitivity, color = life_history, group = life_history)) +
  geom_point(size = 3, aes(color = life_history)) +  
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size = 0.7) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  
  theme(legend.position = "none")



## from Jeffs scripts - ridge plots of posterior distributions ----
fixef_draws <- fit_full %>%
  spread_draws(`b_.*`, regex = TRUE)
fixef_draws 

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

fixef_probs <- fixef_long %>%
  group_by(term_clean) %>%
  summarise(
    prob_diff0 = mean(estimate > 0) %>% {pmax(., 1 - .)},
    .groups = "drop"
  ) %>%
  mutate(
    prob_label = paste0("P ≠ 0: ", round(prob_diff0, 2))
  )

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
