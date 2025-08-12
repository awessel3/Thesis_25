library(lubridate)
library(tidyverse)
library(neonUtilities)
library(sp)
library(sf)
library(elevatr)
library(ggplot2)

setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

############## PROCESS ELEVATION DATA ################### 


# Retrieving Elevation Data R Script 

#reading in All Species -- Test with all species in WVPT region 


#full set
full_range_flr <- read_rds("Data/dat.clim.rds")
head(full_range_flr)
unique(full_range_flr$species)
dim(full_range_flr)
length(unique(full_range_flr$species))

df_flr_final_wClimate <- read_rds("Data/df_flr_final_wClimate.rds")
head(df_flr_final_wClimate)


#df_flr_final_wClimate <- read_rds("Data/df_flr_final_wClimate_JD.rds")


## selecting coordinates ----
#elev <- df_flr_final_wClimate

elev <- full_range_flr

## Temp convert to spatial to get coords ---- 

elev.test <- as.data.frame(dplyr::select(elev, longitude, latitude)) %>%
  rename(x = longitude, y = latitude)
elev.sf <- st_as_sf(elev.test, coords = c("x", "y"), crs = 4326)

elev2 <- get_elev_point(elev.sf, prj = "+proj=longlat +datum=WGS84", src = "aws")
elev_add <- cbind(elev2, elev.test)

elev_add <- elev_add %>% rename(latitude = y, longitude = x)

# Check for duplicates in WVPT_Annual_wclimate

# duplicate_coords_x <- df_flr_final_wClimate %>%
#   group_by(latitude, longitude) %>%
#   filter(n() > 1)

duplicate_coords <- full_range_flr %>%
  group_by(latitude, longitude) %>%
  filter(n() > 1)

# Print results to see if there are duplicates
print(duplicate_coords)

elev_add_unique <- elev_add %>%
  group_by(latitude, longitude) %>%
  summarize(elevation = mean(elevation, na.rm = TRUE)) %>%
  ungroup()

#df_flr_final_complete <- right_join(df_flr_final_wClimate, elev_add_unique, by = c("latitude", "longitude"))

full_range_flr_complete <- right_join(full_range_flr, elev_add_unique, by = c("latitude", "longitude"))

# df_flr_final_complete <- df_flr_final_complete %>% dplyr::select(-geometry)

full_range_flr_complete <- full_range_flr_complete %>% dplyr::select(-geometry)

# dim(df_flr_final_complete)
# str(df_flr_final_complete)

dim(full_range_flr_complete)
str(full_range_flr_complete)


# saveRDS(df_flr_final_complete, file="Data/df_flr_final_complete.rds") 
# write_csv(df_flr_final_complete, file="Data/df_flr_final_complete.csv")

saveRDS(full_range_flr_complete, file="Data/full_range_flr_complete.rds") 
#write_csv(full_range_flr_completee, file="Data/full_range_flr_complete.csv")

#saveRDS(df_flr_final_complete, file="Data/df_flr_final_complete_JD.rds") 
#write_csv(df_flr_final_complete, file="Data/df_flr_final_complete_JD.csv")

