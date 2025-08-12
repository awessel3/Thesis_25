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

setwd("~/Desktop/Thesis_25")


#df_flr_final_filtered <- read_rds("Data/df_flr_final_filtered.rds")

full_range_flr_filtered <- read_rds("Data/full_range_flr_filtered.rds")


# test.data <- full_range_flr_filtered  %>% dplyr::select(latitude, longitude, species, elevation, 
#                                                       doy, precip, temp, life_history,
#                                                       spring_temp, spring_precip)

test.data <- full_range_flr_filtered
test.data <- na.omit(test.data)
dim(test.data)
str(test.data)


#frequency of observations

#general precip and temp 
precip <- ggplot(data = test.data, 
                 aes(x = precip, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
precip

temp <- ggplot(data = test.data, 
               aes(x = temp, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
temp



#Total observations per species 

species_sum <- test.data %>% 
  group_by(species) %>% 
  summarise(total_observations = n()) 
species_sum

#total_observations_plot
ggplot(species_sum, aes(x=  total_observations, y = species)) + geom_col()

#Exploratory just plotting against each other
#print out for each species 

precip <- ggplot(data = test.data, 
                 aes(x = spring_precip, y = doy)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
precip

temp <- ggplot(data = test.data, 
                 aes(x = spring_temp, y = doy)) + 
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

#New - All_exploreplots_spring, replaced preceding_temp with spring_temp
lat + elev + temp + precip

#Against eachother
precip_temp <- ggplot(data = test.data, 
                 aes(x = spring_precip, y = spring_temp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


temp_lat <- ggplot(data = test.data, 
               aes(x = spring_temp, y = latitude)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

elev_temp <- ggplot(data = test.data, 
               aes(x = elevation, y = spring_temp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


lat_precip <- ggplot(data = test.data, 
              aes(x = spring_precip, y = latitude)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

#All_cross_exploreplots

#New: All_cross_exploreplots_spring
precip_temp + temp_lat + elev_temp + lat_precip

#Species Specific

species <- unique(test.data$species)

#temp
#SS_TempDOY_explore
#New: SS_TempDOY_explore_spring
ggplot(test.data, aes( x = spring_temp, y = doy)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~species)
    
#precip
#SS_precipDOY_explore
#New: SS_precipDOY_explore_spring
ggplot(test.data, aes( x = spring_precip, y = doy)) +
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
#New: SS_TempPrecip_explore_spring
ggplot(test.data, aes( x = spring_precip, y = spring_temp)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)


#SS_templat_explore
#New: SS_templat_explore_spring
ggplot(test.data, aes( x = spring_temp, y = latitude)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#SS_elevTemp_explore
#New: SS_elevTemp_explore_spring
ggplot(test.data, aes( x = elevation, y = spring_temp)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)

#SS_latPrecip_explore
#New: SS_latPrecip_explore_spring
ggplot(test.data, aes( x = latitude, y = spring_precip)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)






