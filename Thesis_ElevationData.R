library(lubridate)
library(tidyverse)
library(neonUtilities)
library(sp)
library(sf)
library(elevatr)
library(ggplot2)

setwd("/Users/avawessel/Desktop/Thesis")
getwd()

# Retrieving Elevation Data R Script 

#reading in AllSpecies -- Test with all species in WVPT region 
WVPT_Annual_wclimate <- read_rds("Data/WVPT_Annual_wclimate.rds")

#Filtering Species that have at least __ observations ? future addition 


## selecting coordinates ----
elev <- WVPT_Annual_wclimate

## Temp convert to spatial to get coords ---- 


elev.test <- as.data.frame(dplyr::select(elev, longitude, latitude)) %>%
  rename(x = longitude, y = latitude)
elev.sf <- st_as_sf(elev.test, coords = c("x", "y"), crs = 4326)

elev2 <- get_elev_point(elev.sf, prj = "+proj=longlat +datum=WGS84", src = "aws")
elev_add <- cbind(elev2, elev.test)

elev_add <- elev_add %>% rename(latitude = y, longitude = x)

# Check for duplicates in WVPT_Annual_wclimate
duplicate_coords_x <- WVPT_Annual_wclimate %>%
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

WVPT_Annual_complete <- right_join(WVPT_Annual_wclimate, elev_add_unique, by = c("latitude", "longitude"))

WVPT_Annual_complete <- WVPT_Annual_complete %>% dplyr::select(-private_latitude, -private_longitude,
                                                        -species_guess, -common_name, -iconic_taxon_name,
                                                         )
saveRDS(WVPT_Annual_complete, file="Data/WVPT_Annual_complete.rds") 

write_csv(WVPT_Annual_complete, file="Data/WVPT_Annual_complete.csv")


## Exploratory SPecies Sum table for 11/20

species_sum <- WVPT_Annual_complete %>% 
  group_by(species) %>% 
  summarise(sum = n()) %>% 
  arrange(desc(sum))


ggplot(species_sum, aes(x = sum, y = reorder(species, sum), fill = species)) + geom_col() +
  labs(y = "species", x = "total observations")
