library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(SPEI)
library(CoordinateCleaner)


setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

############## DATA CLEANING ################### 

df_flr_final_complete <- read_rds("Data/df_flr_final_complete.rds")
full_range_flr_complete <- read_rds("Data/full_range_flr_complete.rds")

head(df_flr_final_complete)
colnames(df_flr_final_complete)
dim(df_flr_final_complete)

## Adding life history

life_hist <- read.csv("Data/Flowering_WVPT_life_history.csv")
life_hist
life_hist_species <- unique(life_hist$species)
length(life_hist_species)

full_range_species <- unique(full_range_flr_complete$species)
length(full_range_species)
full_range_species


## Part 1: Adding temp/precip variables 
# - not needed for full range dataset

# feb - may 
spring_temp <- rowMeans(df_flr_final_complete[,c("tmean_6","tmean_7", "tmean_8",
                                                 "tmean_9")])

#df_flr_final_complete$spring_temp <- spring_temp

spring_precip <- rowMeans(df_flr_final_complete[,c("ppt_6","ppt_7", "ppt_8",
                                                   "ppt_9")])

#df_flr_final_complete$spring_precip <- spring_precip

# plot check
ggplot(df_flr_final_complete, aes(x = spring_temp, y = doy)) + geom_point()
ggplot(df_flr_final_complete, aes(x = spring_precip, y = doy)) + geom_point()

## Part 2: Observation bias filtering 

# coords_filter <- df_flr_final_complete 
coords_filter <- as.data.frame(full_range_flr_complete)
str(coords_filter)
coords_filter$dataset <- "clean"

# geographic cleaning 
coords.test <- clean_coordinates(coords_filter, lat = "latitude", lon = "longitude",
                                 tests = c("capitals",
                                       "centroids","equal",
                                       "gbif",
                                       "outliers", "seas",
                                       "zeros"))
summary(coords.test)

# coordinate conversion cleaning 
clean_dataset(coords_filter, lat = "latitude", lon = "longitude",  ds = "dataset") #returning nothing problematic

tests <- c(".zer",".cap",".cen", ".otl", ".sea", ".gbf", ".summary", '.equ', '.val')

coords.test2 <- coords.test %>%
  filter(rowSums(dplyr::select(., all_of(tests)) == FALSE) == 0)
dim(coords.test)
dim(coords.test2)

# join data with removed coordinates

colnames(coords.test2)
# df_flr_final_filtered <- right_join(df_flr_final_complete, coords.test2, by = 'id')
# dim(df_flr_final_filtered)

full_range_flr_complete <- coords.test2
full_range_flr_complete <- full_range_flr_complete %>% dplyr::select(-.val, -.equ, -.zer, -.cap, -.cen, -.sea, -.otl,
                                                                     -.gbf, -.summary, -dataset)

full_range_flr_filtered <- left_join(full_range_flr_filtered, life_hist, by = "species")
head(full_range_flr_complete)
dim(full_range_flr_complete)
length(unique(full_range_flr_complete$species))


## Part 4: filtering for odd observations

excluded <- c("Sidalcea malviflora", 'Bidens frondosa')
specific_id <- c("126237225")
# sum(df_flr_final_complete$id == "126237225")

sum(full_range_flr_complete$id == "126237225")

# df_flr_final_filtered <- df_flr_final_complete %>%
#   filter(!(species %in% excluded)) %>% 
#   filter(id != specific_id)
# dim(df_flr_final_filtered)
# sum(df_flr_final_filtered$id == "126237225")
# sum(df_flr_final_filtered$species == 'Bidens frondosa')

#saving filtered data 
# dim(df_flr_final_complete)
# dim(df_flr_final_filtered)
# colnames(df_flr_final_complete)
# colnames(df_flr_final_filtered)
# 
# saveRDS(df_flr_final_filtered, file="Data/df_flr_final_filtered.rds") 
# write_csv(df_flr_final_filtered, file="Data/df_flr_final_filtered.csv") 

full_range_flr_filtered <- full_range_flr_complete %>%
  filter(!(species %in% excluded)) %>% 
  filter(id != specific_id)
dim(full_range_flr_filtered)
sum(full_range_flr_filtered$id == "126237225")
sum(full_range_flr_filtered$species == 'Bidens frondosa')

# renaming for consistency 

full_range_flr_filtered <- rename(full_range_flr_filtered, temp = mean_temp_fullyear)
full_range_flr_filtered <- rename(full_range_flr_filtered, precip = total_precip_fullyear)

full_range_flr_filtered <- rename(full_range_flr_filtered, spring_precip = total_precip_winterspring)
full_range_flr_filtered <- rename(full_range_flr_filtered, spring_temp = mean_temp_winterspring)

dim(full_range_flr_filtered)
dim(full_range_flr_complete)
dim(full_range_flr_filtered)
colnames(full_range_flr_complete)
colnames(df_flr_final_complete)
colnames(full_range_flr_filtered)


saveRDS(full_range_flr_filtered, file="Data/full_range_flr_filtered.rds")
write_csv(full_range_flr_filtered, file="Data/full_range_flr_filtered.csv")




## OLD -----

#Adding SPEI 

# creating function to reset month numbers 

remap_months <- setNames(1:12, c(5:12, 1:4))

df_flr_final_temp <- df_flr_final_summary %>%
  rename_with(.fn = function(cols) {
    is_tmean <- grepl("^tmean_\\d+$", cols)
    new_cols <- cols
    new_cols[is_tmean] <- paste0(
      "tmean_",
      remap_months[as.character(as.numeric(sub("tmean_", "", cols[is_tmean])))]
    )
    new_cols
  })
df_flr_final_temp

df_flr_final_precip <- df_flr_final_summary %>%
  rename_with(.fn = function(cols) {
    is_ppt <- grepl("^ppt_\\d+$", cols)
    new_cols <- cols
    new_cols[is_ppt] <- paste0(
      "ppt_",
      remap_months[as.character(as.numeric(sub("ppt_", "", cols[is_ppt])))]
    )
    new_cols
  })
df_flr_final_precip


# Adding temperature and precipitation together in long format

df_temp <- df_flr_final_temp %>%
  pivot_longer(
    cols = starts_with("tmean_"), 
    names_to = "temp_month", 
    names_prefix = "tmean_", 
    values_to = "temp_value"
  ) %>%
  mutate(month = as.numeric(gsub("tmean_", "", temp_month))) 



df_precip <- df_flr_final_precip %>%
  pivot_longer(
    cols = starts_with("ppt_"), 
    names_to = "precip_month", 
    names_prefix = "ppt_", 
    values_to = "precip_value"
  ) %>% 
  mutate(month = as.numeric(gsub("ppt_", "", precip_month)))  
dim(df_precip)


df_combined <- df_temp %>%
  left_join(df_precip %>% dplyr::select(id, month, precip_value), by = c("id", "month"))
dim(df_combined)

df_combined <- df_combined %>%  rename(observed_month = month) %>% 
  rename(month = temp_month)
colnames(df_combined)
str(df_combined)

df_combined$month <- as.numeric(df_combined$month)
str(df_combined)

# get monthly values for just first observation
monthly <- df_combined %>%
  filter(id == unique(id)[1]) %>%
  mutate(
    pet = thornthwaite(temp_value, latitude[1]),
    bal = precip_value - pet,  
    bal = as.vector(bal)       
  )
monthly <- monthly %>%
  arrange(observed_month)
monthly
dim(df_combined)
dim(monthly)
colnames(monthly)
str(monthly)


dat1 <- ts(monthly, start = c(2018, 12), frequency = 12)  
dat1 <- as.data.frame(dat1)
sum(is.na(monthly$bal))
dat1
str(dat1)

# Calculate the 1-month SPEI
spei1 <- spei(dat1[,'bal'], scale = 1)
dat1$spei1 <- as.vector(spei1$fitted)
spei1
# Calculate the 3-month SPEI
spei3 <- spei(dat1, 3)
dat1_spei3 <- as.data.frame(spei3$fitted)

# Calculate the 6-month SPEI
spei6 <- spei(dat1, 6)
dat1_spei6 <- as.data.frame(spei6$fitted)

# Add the SPEI values to the monthly data
monthly$spei1 <- dat1_spei1$fitted
monthly$spei3 <- dat1_spei3$fitted
monthly$spei6 <- dat1_spei6$fitted
head(monthly)

calculateClimate <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      observation_month = case_when(
        month(observed_on) == 9 ~ 1,  # September = 1
        month(observed_on) == 10 ~ 2, # October = 2
        month(observed_on) == 11 ~ 3, # November = 3
        month(observed_on) == 12 ~ 4, # December = 4
        month(observed_on) == 1  ~ 5, # January = 5
        month(observed_on) == 2  ~ 6, # February = 6
        month(observed_on) == 3  ~ 7, # March = 7
        month(observed_on) == 4  ~ 8, # April = 8
        month(observed_on) == 5  ~ 9, # May = 9
        month(observed_on) == 6  ~ 10, # June = 10
        month(observed_on) == 7  ~ 11, # July = 11
        month(observed_on) == 8  ~ 12  # August = 12
      ),
      
      preceding_months = list(
        c(
          observation_month, 
          ifelse(observation_month - 1 <= 0, 12 + (observation_month - 1), observation_month - 1),   
          ifelse(observation_month - 2 <= 0, 12 + (observation_month - 2), observation_month - 2),
          ifelse(observation_month -3 <= 0, 12 + (observation_month - 3), observation_month - 3)
        )
      ),
      
      preceding_temp = mean(
        unlist(lapply(preceding_months, function(m) {
          col_name <- paste0("tmean_", m)
          if (col_name %in% names(df)) df[[col_name]][cur_group_id()] else NA
        })),
        na.rm = TRUE
      ),
      
      preceding_precip = mean(
        unlist(lapply(preceding_months, function(m) {
          col_name <- paste0("ppt_", m)
          if (col_name %in% names(df)) df[[col_name]][cur_group_id()] else NA
        })),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>% 
    dplyr::select(-preceding_months, -observation_month)
}


unique(df_flr_final_complete$species)
length(unique(df_flr_final_complete$species))

unique(df_flr_final_complete$species)
length(unique(df_flr_final_complete$species))


#WVPT_climate_summary <- calculateClimate(WVPT_Annual_complete)
df_flr_final_summary <- calculateClimate(df_flr_final_complete)

ggplot(df_flr_final_complete, aes(x = preceding_precip, y = doy)) + geom_point()
ggplot(df_flr_final_complete, aes(x = preceding_temp, y = doy)) + geom_point()




