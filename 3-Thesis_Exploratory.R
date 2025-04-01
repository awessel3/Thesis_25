library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)
library(ggplot2)
library(ggeffects)
library(GGally)


# Analysis R Script 
library(brms)
library(lme4)

setwd("~/Desktop/Thesis_25")

#WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")
df_flr_final_summary <- read_rds("Data/df_flr_final_summary.rds")

test.data <- df_flr_final_summary  %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                               preceding_precip, elevation, doy, precip, temp)
test.data <- na.omit(test.data)

ggplot(data = test.data, 
       aes(x = preceding_temp, y = doy)) + 
  geom_point() 

#Exploratory just plotting against each other
#print out for each species 

precip <- ggplot(data = test.data, 
                 aes(x = preceding_precip, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
precip

temp <- ggplot(data = test.data, 
                 aes(x = preceding_temp, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
temp

elev <- ggplot(data = test.data, 
               aes(x = elevation, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
elev

lat <- ggplot(data = test.data, 
               aes(x = latitude, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
lat

#All_exploreplots
lat + elev + temp + precip

#Against eachother
precip_temp <- ggplot(data = test.data, 
                 aes(x = preceding_precip, y = preceding_temp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


temp_lat <- ggplot(data = test.data, 
               aes(x = preceding_temp, y = latitude)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

elev_temp <- ggplot(data = test.data, 
               aes(x = elevation, y = preceding_temp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


lat_precip <- ggplot(data = test.data, 
              aes(x = preceding_precip, y = latitude)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

#All_cross_exploreplots
precip_temp + temp_lat + elev_temp + lat_precip

#Species Specific

species <- unique(test.data$species)

#temp
#SS_TempDOY_explore
ggplot(test.data, aes( x = preceding_temp, y = doy)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~species)
    
#precip
#SS_precipDOY_explore
ggplot(test.data, aes( x = preceding_precip, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#lat
#SS_latDOY_explore
ggplot(test.data, aes( x = latitude, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#elev
#SS_elevDOY_explore
ggplot(test.data, aes( x = elevation, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#Species Specific Cross-interaction
#SS_TempPrecip_explore
ggplot(test.data, aes( x = preceding_precip, y = preceding_temp)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)


#SS_templat_explore
ggplot(test.data, aes( x = preceding_temp, y = latitude)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#SS_elevTemp_explore
ggplot(test.data, aes( x = elevation, y = preceding_temp)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#SS_latPrecip_explore
ggplot(test.data, aes( x = latitude, y = preceding_precip)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)


# Exploratory with simpler regressions

model1 <- lm(doy ~ 1 + preceding_temp + preceding_precip + elevation + latitude + (1 | species),
             data = test.data)

model <- doy ~ 1 + preceding_temp + preceding_precip + elevation + latitude + (1|species)

fit1 <- lmer(model, test.data)


summary(fit1)$coefficients

predictors <- data.frame(
  pretemp = test.data$preceding_temp,
  preprecip = test.data$preceding_precip,
  precip = test.data$precip,
  temp = test.data$temp,
  elev = test.data$elevation,
  lat = test.data$latitude,
  long = test.data$longitude
)

cortest_results <- list()

for (i in colnames(predictors)) {
  for (j in colnames(predictors)) {
    if (i != j) {
    test_result <- cor.test(predictors[[i]], predictors[[j]], method = "pearson")
    cortest_results[[paste(i, j, sep = "_vs_")]] <- test_result
    }
  }
}

cortest_results


all.var <- ggpredict(fit1, terms = c("preceding_temp", "preceding_precip", "elevation", "latitude"))


ggplot(all.var, aes(x = x, y = predicted)) +
  geom_line() + geom_path(linewidth = 1, color = "seagreen") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "seagreen", alpha = 0.2) +
  labs( y = "doy", x = "temp + precip + elev + lat")

test.data$composite <- with(test.data, 
                            0.25 * preceding_temp + 0.25 * preceding_precip + 0.25 * elevation + 0.25 * latitude)

composite_pred <- ggpredict(fit1, terms = c("composite"))

ggplot(composite_pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1, color = "seagreen") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "seagreen", alpha = 0.2) +
  labs(y = "doy", x = "Composite Predictor") +
  theme_minimal()


model2 <- doy ~ 1 + preceding_temp + preceding_precip + elevation + latitude + 
  preceding_temp * latitude(1| species)






