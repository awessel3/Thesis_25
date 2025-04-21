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

theme_set(theme_minimal())

# setwd("~/Desktop/Thesis_25")

# (moved filtering to processing script)

# read data ----
#WVPT_climate_summary <- read_rds("Data/WVPT_climate_summary.rds")
dat <- read_rds("Data/df_flr_final_filtered_JD.rds")

print(dat,n=3,width=Inf)
     
# Total observations per species ----

species_sum <- dat %>% 
  group_by(species, life_history) %>% 
  summarise(total_observations = n()) 
species_sum

head(species_sum,2)

#total_observations_plot
ggplot(species_sum, aes(x=  total_observations, y = fct_reorder(species, total_observations), fill = life_history)) + 
  geom_col()+
  geom_text(aes(label = total_observations),
            hjust = -0.1,        # push them just to the right of the bar
            size  = 3) +         # adjust text size as needed
  # make room on the right so labels aren't clipped
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x    = "Total observations",
       y    = "Species",
       fill = "Life history") 
  ggsave(filename = "../Ava_iNat_SHARED/Figures/raw_obs.pdf", width = 6, height = 6)


  
  
# Pairs plots ----
print(dat,n=2,width=Inf)

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue", aes=.2) +
    geom_smooth(method = method, color = "red", ...)
  p
}

ggpairs(dat, columns=c("latitude","elevation","preceding_temp","preceding_precip","doy" ), lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "blue")),
        upper = list(continuous = wrap("cor", size = 10)))
ggsave(filename = "../Ava_iNat_SHARED/Figures/raw_correlations.pdf", width = 9, height = 9)


ggpairs(dat, columns=c("latitude","elevation","winter.temp","winter.precip","spring.temp","spring.precip","doy" ), lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "blue")),
        upper = list(continuous = wrap("cor", size = 10)))
ggsave(filename = "../Ava_iNat_SHARED/Figures/raw_correlations_seasonal.pdf", width = 9, height = 9)


#Species Specific

species <- unique(dat$species)

#temp
#SS_TempDOY_explore
ggplot(dat, aes( x = spring.temp, y = doy)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~species)
ggsave(filename = "../Ava_iNat_SHARED/Figures/spring_temp_sp.pdf", width = 9, height = 9)

#precip
#SS_precipDOY_explore
ggplot(dat, aes( x = spring.precip, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)
ggsave(filename = "../Ava_iNat_SHARED/Figures/spring_precip_sp.pdf", width = 9, height = 9)

#lat
#SS_latDOY_explore
ggplot(dat, aes( x = latitude, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)
ggsave(filename = "../Ava_iNat_SHARED/Figures/lat_sp.pdf", width = 9, height = 9)

#elev
#SS_elevDOY_explore
ggplot(dat, aes( x = elevation, y = doy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)
ggsave(filename = "../Ava_iNat_SHARED/Figures/elev_sp.pdf", width = 9, height = 9)

#Species Specific Cross-interaction
#SS_TempPrecip_explore
ggplot(dat, aes( x = preceding_precip, y = preceding_temp)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) + 
  facet_wrap(~species)
ggsave(filename = "../Ava_iNat_SHARED/Figures/temp_precip.pdf", width = 9, height = 9)




# Exploratory with simpler regressions

library(lme4) # allows mixed effects
library(lmerTest) # gets p-values from mixed effects models

model1 <- lmer(doy ~ 1 + preceding_temp * preceding_precip + elevation + latitude + (1 | species),
               data = dat)
summary(model1)

model2 <- lmer(doy ~ 1 + preceding_temp * preceding_precip + elevation + latitude + (preceding_temp * preceding_precip | species),
               data = dat)
summary(model2)



