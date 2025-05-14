library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)
library(loo)
library(RColorBrewer)



# Analysis R Script 
library(brms)
library(ggridges)
library(ggdist)
library(rlang)
library(stringr)

setwd("~/Desktop/Thesis_25")

df_flr_final_summary <- read_rds("Data/df_flr_final_summary.rds")
df_flr_final_filtered <- read_rds("Data/df_flr_final_filtered.rds")
dim(df_flr_final_filtered)

## old trait data prep ---- 
#trait_species <- read.csv("trait_species.csv")
#traits_full <- read.csv('Traits.csv')
#traits_full$SpName <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", traits_full$SpName)

#traits_full <- rename(traits_full, species = SpName)


#life_hist <- read.csv("Data/Flowering_WVPT_life_history.csv")
#spp <- unique(df_flr_final_filtered$species)

#df_traits_flr <- traits_full %>%
# filter(species %in% spp) 
#df_traits_flr_final <- left_join(df_traits_flr, life_hist, by = "species")

#length(unique(df_traits_flr_final$species))

#trait_hist_count <- df_traits_flr_final %>% 
#  group_by(life_history) %>% 
#  summarise(count = n())
# trait_hist_count <- count(df_traits_flr_final, c("annual"))
# trait_hist_count

## Prepping data for model ---- 

unique(df_flr_final_filtered$species)

#Excluding species - either for low number of observations or skewing data 


data <- df_flr_final_filtered %>% dplyr::select(latitude, longitude, species, spring_temp,
                                               spring_precip, elevation, doy, life_history)
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

stemp_num <- as.numeric(data$spring_temp)
stemp_center <- mean(stemp_num)
stemp_scale <- sd(stemp_num)
data$stemp_sc <- (stemp_num - stemp_center) / stemp_scale


sprecip_num <- as.numeric(data$spring_precip)
sprecip_center <- mean(sprecip_num)
sprecip_scale <- sd(sprecip_num)
data$sprecip_sc <- (sprecip_num - sprecip_center) / sprecip_scale     

elevation_num <- as.numeric(data$elevation)
elevation_center <- mean(elevation_num)
elevation_scale <- sd(elevation_num)
data$elevation_sc <- (elevation_num - elevation_center) / elevation_scale

# fit0
formula <- doy_sc ~ 1 + stemp_sc * latitude_sc + stemp_sc * elevation_sc +
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

#FULL MODEL
formula_full <- doy_sc ~ 1 + latitude_sc * stemp_sc * elevation_sc * sprecip_sc * life_history +
  (1 + latitude_sc * stemp_sc * elevation_sc * sprecip_sc | species)

#
formula_select1 <- doy_sc ~ 1 + 
  latitude_sc * stemp_sc * sprecip_sc * life_history + elevation + 
  (1 + latitude_sc * stemp_sc * sprecip_sc + elevation_sc | species)

fit <- brm(
  formula = formula_select1,
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

#Pass 5 save - full fit 
saveRDS(fit, file = "Data/fit_full.RDS")
fit_full <- readRDS("Data/fit_full.RDS") 

saveRDS(fit, file = "Data/fit_fullselect1.RDS")
fit_fullselect1 <- readRDS("Data/fit_fullselect1.RDS") 
fit_fullselect1

loo1 <- loo(fit_full)
loo2 <- loo(fit_fullselect1)
loo <- loo(fit_full, fit_fullselect1)
loo


#prep for plotting 
fit <- fit_full
original_data <- fit$data 
colnames(original_data)
spp <- unique(original_data$species)
spp
life_history <- unique(original_data$life_history)
life_history

params <- get_variables(fit)
write_csv(tibble(parameter = params), file = "Analysis_Output/fit_full_params.csv")

## Initial Plot Creation ----
# Create a new data set to predict over
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc), # predictions at mean elevation
  stemp_sc = mean(original_data$stemp_sc), # temperature range
  sprecip_sc = mean(original_data$sprecip_sc, na.rm = TRUE), # mean precipitation
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
lat <- ggplot(fitted.pred, aes(x = latitude, y = DOY_pred)) +
  geom_point(data = data, aes(x = latitude, y = doy), color = "grey", alpha = 0.4) +
  stat_lineribbon(aes(y = DOY_pred), .width = c(0.5, 0.9), 
                  color = 'dodgerblue1', alpha = 0.5) +
  labs(y = NULL, x = "Latitude") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

lat

#Temp vs DOY
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc), # predictions at mean elevation
  stemp_sc = seq(min(original_data$stemp_sc), max(original_data$stemp_sc), length.out = 100), # temperature range
  sprecip_sc = mean(original_data$sprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
  )


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring_temp = (stemp_sc * stemp_scale) + stemp_center)



# Plot
# fit#_TempDOY_plot
temp <- ggplot(fitted.pred, aes(x = spring_temp , y = DOY_pred)) +
  geom_point(data = data, aes(x = spring_temp, y = doy), color = "grey", alpha = 0.4) +
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE, color = 'dodgerblue1', alpha = 0.5) +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  labs(x = "Spring Temperature", y= "Day of Flowering (DOY)") +
  theme_minimal()
temp

#Precip vs DOY
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc), 
  stemp_sc = mean(original_data$stemp_sc),
  sprecip_sc = seq(min(original_data$stemp_sc), max(original_data$stemp_sc), length.out = 100), 
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
)


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring_precip = (sprecip_sc * sprecip_scale) + sprecip_center)



# Plot
# fit#_precipDOY_plot
precip <- ggplot(fitted.pred, aes(x = spring_precip , y = DOY_pred))+
  geom_point(data = data, aes(x = spring_precip, y = doy), color = "grey", alpha = 0.4) +
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE, color = 'dodgerblue1', alpha = 0.5) +
  labs(x = "Spring Precipitation", y  = NULL) +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()
precip


#Temp|latitude vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc, na.rm = TRUE), # predictions at mean elevation
  stemp_sc = seq(min(original_data$stemp_sc), max(original_data$stemp_sc), length.out = 100), # temperature range
  sprecip_sc = mean(original_data$sprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = seq(min(original_data$latitude_sc), max(original_data$latitude_sc), length.out = 5),
  species = spp,
  life_history = life_history
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring_temp = (stemp_sc * stemp_scale) + stemp_center,
         latitude = (latitude_sc * latitude_scale) + latitude_center)

#color palette 
palette_blues <- colorRampPalette(colors = c("white", "#004b88"))(12)
scales::show_col(palette_blues[4:12])

# Plot
# fit#_TempLatDOY_plot
temp_lat <- ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(latitude, 1)))) + 
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE,  alpha = 0.5) + 
  scale_color_brewer(palette = "Dark2", name = "Latitude") +  
  scale_fill_brewer(palette = "Greys", guide = "none") +
  labs(x = "Spring Temperature", y = "Day of Flowering (DOY)") +
  theme_minimal()  
temp_lat


# SS_fit#_TempLatDOY_plot
SS_temp_lat <- ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(latitude, 1)))) + 
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE,  alpha = 0.5) + 
  scale_color_brewer(palette = "Dark2", name = "Latitude") +  
  scale_fill_brewer(palette = "Greys", guide = "none") +
  labs(x = "Spring Temperature", y = NULL) +
  facet_wrap(~species) + 
  theme_minimal()  
SS_temp_lat

#Temp|Elevation vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation_sc = seq(min(original_data$elevation_sc), max(original_data$elevation_sc), length.out = 5), # predictions at mean elevation
  stemp_sc = seq(min(original_data$stemp_sc), max(original_data$stemp_sc), length.out = 100), # temperature range
  sprecip_sc = mean(original_data$sprecip_sc, na.rm = TRUE), # mean precipitation
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring_temp = (stemp_sc * stemp_scale) + stemp_center,
         elevation = (elevation_sc * elevation_scale) + elevation_center)

# Plot
#fit#_TempElevDOY_plot
temp_elev <- ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(elevation, 1)))) + 
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE, alpha = 0.5) +  
  scale_color_brewer(palette = "Set1", name = "Elevation") +  
  scale_fill_brewer(palette = "Greys", guide = "none") +
  labs(x = "Spring Temperature", y = "Day of Flowering (DOY)") +
  theme_minimal()  
temp_elev

# SS_fit#_TempElevDOY_plot
ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(elevation, 1)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()

#Temp|Precip vs DOY
data.predict <- crossing(
  elevation_sc = mean(original_data$elevation_sc, na.rm = TRUE), 
  stemp_sc = seq(min(original_data$stemp_sc), max(original_data$stemp_sc), length.out = 100), 
  sprecip_sc =  seq(min(original_data$sprecip_sc), max(original_data$sprecip_sc), length.out = 5), 
  latitude_sc = mean(original_data$latitude_sc, na.rm = TRUE),
  species = spp,
  life_history = life_history
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 500, allow_new_levels = TRUE) %>%
  mutate(DOY_pred_sc = .linpred)

fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_sc * doy_scale) + doy_center,
         spring_temp = (stemp_sc * stemp_scale) + stemp_center,
         spring_precip = (sprecip_sc * sprecip_scale) + sprecip_center)

# Plot
#fit#_TempPrecipDOY_plot
temp_precip <- ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(spring_precip, 1)))) + 
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE, alpha = 0.5) +  
  scale_color_brewer(palette = "Set1", name = "Spring Precipitation") +  
  scale_fill_brewer(palette = "Greys", guide = "none") +
  labs(x = "Spring Temperature", y = NULL) +
  theme_minimal() 

# SS_fit#_TempPrecipDOY_plot
ggplot(fitted.pred, aes(x = spring_temp, y = DOY_pred, color = factor(round(spring_precip, 1)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  facet_wrap(~species) +
  theme_minimal()

#Saving large figure 1 
top_plots <- temp + lat + precip
top_plots
ggsave(plot=top_plots,"Analysis_Images/figure1_draft.pdf", width=6, height=4)


bottom_plots <- temp_elev + temp_lat + temp_precip
bottom_plots

all_fig1 <- top_plots / bottom_plots
all_fig1

#ggsave("Analysis_Images/figure_1_draft.pdf", all_fig1,width=7, height=5)


## Density plots - species specific ----

# comparison of species specific slopes to mean doy 
species_draws2 <- read_rds("Data/species_draws2.rds")
species_summary <- read.csv("Analysis_Output/species_summary.csv")
overall_summary <- read.csv("Analysis_Output/overall_summary.csv")
dim(species_summary)
species_summary$species <- gsub("\\.", " ", species_summary$species)
species_summary
colnames(species_summary)


avg_doy <- data %>% 
  group_by(species) %>% 
  summarise(doy = mean(doy))
avg_doy

temp_eff <- species_summary %>% 
  filter(term_lab == "Temperature")
dim(temp_eff)
temp_eff.plot <- left_join(temp_eff, avg_doy, by = "species") 
dim(temp_eff.plot)
temp_eff.plot

overall_model <- lm(mean_est ~ doy, data = temp_eff.plot)

# Get the R-squared value
overall_r_squared <- summary(overall_model)$r.squared

temp_eff.plot

ggplot(temp_eff.plot, aes( x = doy, y = mean_est, color = species)) + geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x, color = "black", size = 0.75) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size = 0.7) +
  annotate("text", x = 150, y = 0.05, 
           label = paste("R² =", round(overall_r_squared, 2)), size = 5, color = "black")
  theme_minimal()
  
## create general plots 

#1 term
colnames(species_summary)

stacked_species_draw.plot <-ggplot(species_draws2, aes(x = total, y = term_lab, fill = species)) +
  geom_density_ridges(
    scale = 1,
    alpha = 0.3,
    quantile_lines = FALSE,
    quantiles = c(0.025, 0.5, 0.975),
  ) +
  labs(y = NULL, 
       x = "Mean Estimate") +
  geom_vline(xintercept = 0, color = "red" ,linetype = "dashed", linewidth = 1) +
  theme_minimal() 
  
stacked_species_draw.plot <- stacked_species_draw.plot + theme(legend.position = "none")
stacked_species_draw.plot

ggsave(plot=stacked_species_draw.plot,"Analysis_Images/full_model/stacked_species_draw.pdf", width=7, height=3.5)
  




## Creating phenological sensitivity plot ----

## Random Slopes just for Temp 
# view variables that can be extracted 
get_variables(fit)

# extracr species specific slopes for temperature
sp_slopes_posterior_temp <- fit %>%
  spread_draws(r_species[species, stemp_sc ])  
head(sp_slopes_posterior_temp)

unique(sp_slopes_posterior_temp$stemp_sc)

ggplot(sp_slopes_posterior_temp, aes(x = r_species, y = species)) + 
  stat_summary(geom = "pointrange", fun.data = mean_hdi) +
  labs(title = "Species-Specific Slopes", y = "Species", x = "Slope") +
  theme_minimal()

species_temp_ps <- fit %>%
  spread_draws(r_species[species, stemp_sc]) %>%
  group_by(species) %>%
  summarise(
    mean_sensitivity = mean(r_species),  
    lower_95 = quantile(r_species, 0.025),
    upper_95 = quantile(r_species, 0.975)
  ) 

# unscaling estimated parameters 
species_temp_ps <- species_temp_ps %>%
  mutate(sensitivity_unscaled = mean_sensitivity * (doy_scale / stemp_scale),
         lower_95_unscaled = lower_95 * (doy_scale / stemp_scale),
         upper_95_unscaled = upper_95 * (doy_scale / stemp_scale))


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
ggplot(temp_ps_plot_dat, aes(x = mean_doy, y = sensitivity_unscaled, color = species)) +
  geom_point(size = 3, aes(color = species)) + 
  stat_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.75) + 
  geom_hline(yintercept = 0) 
  #geom_errorbar(aes(ymin = lower_95_unscaled, ymax = upper_95_unscaled), linewidth = 1, size = 0.7)

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

## making similar plots for other variables 

get_variables(fit)

# extracr species specific slopes for temperature
sp_slopes_posterior_lat <- fit %>%
  spread_draws(r_species[species, b_latitude_sc])  
head(sp_slopes_posterior_lat)

unique(sp_slopes_posterior_lat$b_latitude_sc)

ggplot(sp_slopes_posterior_lat, aes(x = r_species, y = species)) + 
  stat_summary(geom = "pointrange", fun.data = mean_hdi) +
  labs(title = "Species-Specific Slopes", y = "Species", x = "Slope") +
  theme_minimal()

species_lat_ps <- fit %>%
  spread_draws(r_species[species, latitude_sc]) %>%
  group_by(species) %>%
  summarise(
    mean_sensitivity = mean(r_species),  
    lower_95 = quantile(r_species, 0.025),
    upper_95 = quantile(r_species, 0.975)
  ) 

# unscaling estimated parameters 
species_lat_ps <- species_lat_ps %>%
  mutate(sensitivity_unscaled = mean_sensitivity * (doy_scale / stemp_scale),
         lower_95_unscaled = lower_95 * (doy_scale / stemp_scale),
         upper_95_unscaled = upper_95 * (doy_scale / stemp_scale))


#fixing species names containing . between 
species_lat_ps$species <- gsub("\\.", " ", species_lat_ps$species)

# Extract mean flowering
mean_doy <- original_data %>%
  group_by(species) %>%
  summarise(mean_doy = mean(((doy_sc * doy_scale) + doy_center), na.rm = TRUE))

# Merge data sets
lat_ps_plot_dat <- left_join(species_lat_ps, mean_doy, by = 'species')
lat_ps_plot_dat <- left_join(lat_ps_plot_dat, life_hist, by = "species")
lat_ps_plot_dat


# fit#_temp_PS_plot
ggplot(lat_ps_plot_dat, aes(x = mean_doy, y = sensitivity_unscaled, color = species)) +
  geom_point(size = 3, aes(color = species)) + 
  stat_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.75) + 
  geom_hline(yintercept = 0) 
#geom_errorbar(aes(ymin = lower_95_unscaled, ymax = upper_95_unscaled), linewidth = 1, size = 0.7)







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
#fit#_general_ridgeplot
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



# Extract random slopes for each species
fit

posterior::as_draws_df(fit) %>%
  colnames() %>%
  str_subset("^b_")


# Get the random slope term names
slope_terms <- get_variables(fit) %>%
  str_subset("^r_species\\[.*?,(.*?)\\]$") %>%
  str_match("^r_species\\[.*?,(.*?)\\]$") %>%
  .[, 2] %>%
  unique()

# Build expressions like r_species[species, <term>]
random_exprs <- paste0("r_species[species,", slope_terms, "]") %>%
  parse_exprs()

# Build fixed effect names

fixed_exprs <- paste0("b_", slope_terms) %>% syms()
# Extract posterior draws for those fixed effects ## really unsure why this is not working 
fixed_draws <- fit %>%
  spread_draws(!!!fixed_exprs)

# Random effect draws for species 
random_draws <- fit %>%
  spread_draws(!!!random_exprs)

# combine and calculate species-specific slopes

# extract raw term names for all species
all_r_vars <- get_variables(fit) %>%
  str_subset("^r_species\\[.*?,.*\\]$")
head(all_r_vars)

random_draws <- fit %>%
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

