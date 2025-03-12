library(lubridate)
library(tidyverse)
library(neonUtilities)
library(sp)
library(sf)
library(elevatr)
library(ggplot2)

setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

# Retrieving Elevation Data R Script 

#reading in All Species -- Test with all species in WVPT region 

#test set
#WVPT_Annual_wclimate <- read_rds("Data/WVPT_Annual_wclimate.rds")


#full set
df_flr_final_wClimate <- read_rds("Data/df_flr_final_wClimate.rds")


## selecting coordinates ----
elev <- df_flr_final_wClimate

## Temp convert to spatial to get coords ---- 


elev.test <- as.data.frame(dplyr::select(elev, longitude, latitude)) %>%
  rename(x = longitude, y = latitude)
elev.sf <- st_as_sf(elev.test, coords = c("x", "y"), crs = 4326)

elev2 <- get_elev_point(elev.sf, prj = "+proj=longlat +datum=WGS84", src = "aws")
elev_add <- cbind(elev2, elev.test)

elev_add <- elev_add %>% rename(latitude = y, longitude = x)

# Check for duplicates in WVPT_Annual_wclimate
duplicate_coords_x <- df_flr_final_wClimate %>%
  group_by(latitude, longitude) %>%
  filter(n() > 1)

# Check for duplicates in elev_add
duplicate_coords_y <- elev_add %>%
  group_by(latitude, longitude) %>%
  filter(n() > 1)

# Print results to see if there are duplicates
print(duplicate_coords_x)
print(duplicate_coords_y)

elev_add_unique <- elev_add %>%
  group_by(latitude, longitude) %>%
  summarize(elevation = mean(elevation, na.rm = TRUE)) %>%
  ungroup()

df_flr_final_complete <- right_join(df_flr_final_wClimate, elev_add_unique, by = c("latitude", "longitude"))

df_flr_final_complete <- df_flr_final_complete %>% dplyr::select(-geometry)
# -private_latitude, -private_longitude,-species_guess, -common_name, -iconic_taxon_name


saveRDS(df_flr_final_complete, file="Data/df_flr_final_complete.rds") 

write_csv(df_flr_final_complete, file="Data/df_flr_final_complete.csv")

