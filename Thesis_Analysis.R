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

## Prepping data for model 

WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")

data <- WVPT_climate_summary %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                               preceding_precip, elevation, doy)
data <- na.omit(data)

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


#priors <- c(
  #set_prior("normal(0, 10)", class = "b"),          
 # set_prior("normal(0, 100)", class = "Intercept"))


# fit1
formula <- doy ~ 1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc + 
  (1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

#altering model to leave on out to reconginze necessary complexity 

# pass 1: remove cross-level interaction between temp & lat - made no difference in complexity
# fit2
formula1 <- doy ~ 1 + ptemp_sc + latitude_sc + ptemp_sc * elevation_sc +
  pprecip_sc +(1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)

# pass 2: remove cross-level interaction between temp & elev - pass 1 perfomed sig better. 
formula2 <- doy ~ 1 + ptemp_sc + latitude_sc + ptemp_sc + elevation_sc +
  pprecip_sc +(1 + ptemp_sc * latitude_sc + ptemp_sc * elevation_sc + pprecip_sc | species)


fit <- brm(
  formula = formula,
  formula = formula2,
  data = data,
  family = gaussian(),  # Assuming DOY is approximately normally distributed
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  #init = "0"
)

pairs(fit, np = nuts_params(fit))

#original model save
saveRDS(fit, file = "Data/fit1.RDS")
fit1 <- readRDS("Data/fit1.RDS")  

#pass 1 save 
saveRDS(fit, file = "Data/fit2.RDS")
fit2 <- readRDS("Data/fit2.RDS")  

#pass 2 save, saving issue, am not rerunning model, therefore not saved. 
#saveRDS(fit, file = "Data/fit3.RDS")
#fit3 <- readRDS("Data/fit3.RDS")  



loo1 <- loo(fit)
loo1

 
as_draws_df(fit) %>% head(3)

summary(fit)
plot(fit)
lp_draws <- linpred_draws(fit, newdata = original_data)


original_data <- fit$data 
spp <- unique(original_data$species)

## Initial Plot Creation ----
# Create a new dataset to predict over
data.predict <- crossing(
  elevation_sc = (mean(original_data$elevation, na.rm = TRUE) - elevation_center) / elevation_scale, 
  ptemp_sc = (mean(original_data$preceding_temp, na.rm = TRUE) - ptemp_center) / ptemp_scale,  
  pprecip_sc = (mean(original_data$preceding_precip, na.rm = TRUE) - pprecip_center) / pprecip_scale,
  latitude_sc = (seq_range(quantile(original_data$latitude, probs = c(.05, .5, .95)), n = 50) - latitude_center) / latitude_scale,
  species = spp  
)


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred_sc = .linpred )

# Unscale the predicted DOY
fitted.pred <- fitted.pred %>%
  mutate(DOY_pred = (DOY_pred_scaled * doy_scale) + doy_center,
         latitude = (latitude_sc * latitude_scale) + latitude_center)  


# Plot
ggplot(fitted.pred, aes(x = latitude , y = DOY_pred, color=species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="Latitude") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()
#ggsave("Figures/DOY_latitude.pdf", width=5, height=4)


#Temp vs DOY
data.predict <- crossing(
  elevation = mean(original_data$elevation) # predictions at mean elevation
  ,preceding_temp = seq(min(original_data$preceding_temp), max(original_data$preceding_temp), length.out = 100) # predictions at mean temp
  ,preceding_precip = mean(original_data$preceding_precip, na.rm=TRUE) # predictions at mean precip
  ,latitude = mean(original_data$latitude, na.rm = TRUE)
  ,species= spp 
)


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred = .linpred )



# Plot
ggplot(fitted.pred, aes(x = preceding_temp , y = DOY_pred, color=species))+
  stat_lineribbon(.width = c(.5,.9),show.legend=TRUE) +
  labs(y="Day of Year flowering", x="preceding temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()


#Temp|latitude vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation = mean(original_data$elevation), # predictions at mean elevation
  preceding_temp = seq(min(original_data$preceding_temp), max(original_data$preceding_temp), length.out = 100), # temperature range
  preceding_precip = mean(original_data$preceding_precip, na.rm = TRUE), # mean precipitation
  latitude = seq(min(original_data$latitude), max(original_data$latitude), length.out = 5),
  species = spp
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred = .linpred)

# Plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(round(latitude, 2)))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()

#Temp|Elevation vs DOY
# Define predictions dataset
data.predict <- crossing(
  elevation = seq(min(original_data$elevation), max(original_data$elevation), length.out = 5),
  preceding_temp = seq(min(original_data$preceding_temp), max(original_data$preceding_temp), length.out = 100), # temperature range
  preceding_precip = mean(original_data$preceding_precip, na.rm = TRUE), # mean precipitation
  latitude = mean(original_data$latitude, na.rm = TRUE),
  species = spp
)

# Make predictions using the fitted model
fitted.pred <- linpred_draws(object = fit, newdata = data.predict, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred = .linpred)

# Plot
ggplot(fitted.pred, aes(x = preceding_temp, y = DOY_pred, color = factor(elevation))) +
  stat_lineribbon(.width = c(0.5, 0.9), show.legend = TRUE) +
  labs(y = "Day of Year Flowering", x = "Preceding Temperature") +
  scale_fill_brewer(palette = "Greys", guide = "none") +
  theme_minimal()



# Creating phenological sensitivity plot ----

#create dataframe with all variable ranges to predict doy 
new_data <- crossing(
  elevation =  seq(min(original_data$elevation), max(original_data$elevation), length.out = 5),
  preceding_temp = seq(min(original_data$preceding_temp), max(original_data$preceding_temp), length.out = 100), # temperature range
  preceding_precip = seq(min(original_data$preceding_precip), max(original_data$preceding_precip), length.out = 5), # mean precipitation
  latitude = seq(min(original_data$latitude), max(original_data$latitude), length.out = 5),
  species = spp
)

#predict doy 
predictions <-  linpred_draws(object = fit, newdata = new_data, ndraws = 1000, allow_new_levels = TRUE) %>%
  mutate(DOY_pred = .linpred)

summary(fit)$fixed

fixef(fit)

ranef(fit)

coef(fit)

posterior <- as_draws_df(fit)
head(posterior)


#extract observed and predicted doy and calculate sensitivity 
sensitivity_draws <- predictions %>%
  group_by(species, .draw) %>%
  summarise(
    observed_doy = mean(data$doy),  
    pred_doy = mean(DOY_pred),    
    sensitivity = observed_doy - pred_doy
  )

#summarize draws and calculate ci
sensitivity_summary <- sensitivity_draws %>%
  group_by(species) %>%
  summarise(
    sensitivity = mean(sensitivity),
    pred_doy = mean(pred_doy),
    observed_doy = mean(observed_doy),
    lower_ci = quantile(sensitivity, 0.025),
    upper_ci = quantile(sensitivity, 0.975)
  )


# Issues with observed doy-- (Something changing) --- calculate doy from original data?
 ## Why does it change? 
ob_doy <- WVPT_climate_summary %>% 
  group_by(species) %>% 
  summarise(observed_doy = mean(doy))

# Join data with true observed doy 
#ps.zero <- full_join(observed_doy, sensitivity_summary, by = "species")
ps <- full_join(ob_doy, sensitivity_summary, by = "species")
ps <- ps %>% select(c(-observed_doy.y)) %>% rename(observed_doy = observed_doy.x)

str(ps)
ps$observed_doy <- as.numeric(ps$observed_doy)
ps$pred_doy <- as.numeric(ps$pred_doy)
ps$sensitivity <- as.numeric(ps$sensitivity)
ps$lower_ci <- as.numeric(ps$lower_ci)
ps$upper_ci <- as.numeric(ps$upper_ci)
str(ps)

ps$observed_doy <- round(ps$observed_doy, 0)


#Create plot
## How can I make accurate errorbars? 
##W What can I gather from this figure? 
ggplot(ps, aes(x = observed_doy, y = sensitivity, color = species)) +
  geom_point(size = 3) +  
  geom_hline(aes(yintercept = 0), color = "red") + 
  #geom_errorbar(aes(x = observed_doy, ymin = lower_ci, ymax = upper_ci), width = 0) + 
  labs(
    x = "Observed DOY (Mean)",
    y = "Phenological Sensitivity (DOY/C)",
    x
  ) +
  theme_minimal()


#old
observed_doy <- data %>% 
  group_by(species) %>% 
  summarise(observed_doy = mean(doy))


pred_doy <- predictions %>%  
  group_by(species) %>% 
  summarise(pred_doy = mean(DOY_pred))
  
ps <- full_join(observed_doy, pred_doy)

ps <- ps %>%
  mutate(sensitivity = observed_doy - pred_doy)

str(ps)
ps$observed_doy <- as.numeric(ps$observed_doy)
ps$pred_doy <- as.numeric(ps$pred_doy)
ps$sensitivity <- as.numeric(ps$sensitivity)
