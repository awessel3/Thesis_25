library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)


# Analysis R Script 
library(brms)

WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")


#WVPT_Annual_complete <- WVPT_Annual_complete %>%
 # rowwise() %>%
 # mutate(springTemp = mean(c_across(c("tmean_6", "tmean_7", "tmean_8", "tmean_9")), na.rm = TRUE)) %>%
 # ungroup()

priors <- c(
  set_prior("normal(0, 500)", class = "b"),           # Prior for fixed effects (slope)
  set_prior("normal(0, 5)", class = "Intercept")  # Prior for the intercept
)

formula <- doy ~ 1 + elevation + preceding_temp + preceding_precip + latitude + 
  preceding_temp * latitude + preceding_temp * elevation + (1 | species) 

formula2 <- doy ~ (1 + elevation + preceding_temp + preceding_precip + latitude | species) 

formula3 <-  doy ~ 1 + elevation + preceding_temp + preceding_precip + latitude + (1 | species) 


data <- WVPT_climate_summary %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                               preceding_precip, elevation, doy)

WVPT_climate_summary <- na.omit(WVPT_climate_summary)

fit <- brm(
  formula = formula2,
  data = WVPT_climate_summary,
  family = gaussian(),  # Assuming DOY is approximately normally distributed
  #prior = priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  init = "0"
)

 
as_draws_df(fit) %>% head(3)

fit
summary(fit)
plot(fit)

original_data <- fit$data 
lp_draws <- linpred_draws(fit, newdata = original_data)

# doy ~ 1 + elevation + preceding_temp + preceding_precip + latitude + preceding_temp 
#* latitude + preceding_temp * elevation + (1 | species) 


spp <- unique(original_data$species)

## Initial Plot Creation ----
# Create a new dataset to predict over
data.predict <- crossing(
  elevation = mean(original_data$elevation) # predictions at mean elevation
  ,preceding_temp = mean(original_data$preceding_temp, na.rm=TRUE) # predictions at mean temp
  ,preceding_precip = mean(original_data$preceding_precip, na.rm=TRUE) # predictions at mean precip
  ,latitude = round(as.numeric(seq_range(quantile(original_data$latitude, probs=c(.05,.5,.95)), n=50)), digits=3)
  ,species= spp 
)


# make predictions using the fitted model
fitted.pred <-  linpred_draws(object=fit, newdata=data.predict, ndraws=1000, allow_new_levels=TRUE) %>%
  mutate(DOY_pred = .linpred )



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
