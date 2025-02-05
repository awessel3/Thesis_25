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

WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")


priors <- c(
  set_prior("normal(0, 500)", class = "b"),           # Prior for fixed effects (slope)
  set_prior("normal(0, 5)", class = "Intercept")  # Prior for the intercept
)

test.data <- WVPT_climate_summary %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                               preceding_precip, elevation, doy)


#Exploratory just plotting against each other

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


precip_temp + temp_lat + elev_temp + lat_precip



# Exploratory with simpler regressions

model1 <- lm(doy ~ 1 + preceding_temp + preceding_precip + elevation + latitude,
             data = test.data)

model <- doy ~ 1 + preceding_temp + preceding_precip + elevation + latitude + (1|species)

fit1 <- lmer(model1, test.data)

test.data$predicted_doy <- predict(model1)

predict(model1)

summary(model1)

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






