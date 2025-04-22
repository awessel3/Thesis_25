library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(SPEI)


setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

#WVPT_Annual_complete <- read_rds("Data/WVPT_Annual_complete.rds")

df_flr_final_complete <- read_rds("Data/df_flr_final_complete.rds")
head(df_flr_final_complete)
colnames(df_flr_final_complete)

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

# regular temp/precip no preceding, fixed window 
# feb - may 
spring_temp <- rowMeans(df_flr_final_complete[,c("tmean_6","tmean_7", "tmean_8",
                                                 "tmean_9")])

df_flr_final_summary$spring_temp <- spring_temp

spring_precip <- rowMeans(df_flr_final_complete[,c("ppt_6","ppt_7", "ppt_8",
                                                 "ppt_9")])

df_flr_final_summary$spring_precip <- spring_precip


ggplot(df_flr_final_summary, aes(x = spring_temp, y = doy)) + geom_point()
ggplot(df_flr_final_summary, aes(x = spring_precip, y = doy)) + geom_point()
ggplot(df_flr_final_summary, aes(x = preceding_precip, y = doy)) + geom_point()
ggplot(df_flr_final_summary, aes(x = preceding_temp, y = doy)) + geom_point()


# saving full files 
saveRDS(df_flr_final_summary, file="Data/df_flr_final_summary.rds") 
write_csv(df_flr_final_summary, file="Data/df_flr_final_summary.csv") 


## Data Cleaning --- 
df_flr_final_summary <- read_rds("Data/df_flr_final_summary.rds")
excluded <- c("Sidalcea malviflora", 'Bidens frondosa')
specific_id <- c("126237225")
sum(df_flr_final_summary$id == "126237225")

df_flr_final_filtered <- df_flr_final_summary %>%
  filter(!(species %in% excluded)) %>% 
  filter(id != specific_id)
dim(df_flr_final_filtered)
sum(df_flr_final_filtered$id == "126237225")
sum(df_flr_final_filtered$species == 'Bidens frondosa')

test.data <- df_flr_final_filtered  %>% dplyr::select(latitude, longitude, species, preceding_temp,
                                                      preceding_precip, elevation, doy, precip, temp, life_history,
                                                      spring_temp, spring_precip)
test.data <- na.omit(test.data)
dim(test.data)

#saving filtered data 
saveRDS(df_flr_final_filtered, file="Data/df_flr_final_filtered.rds") 
write_csv(df_flr_final_filtered, file="Data/df_flr_final_filtered.csv") 


#Adding SPEI 

# Adding temperature and precipitation together in long format

df_temp <- df_flr_final_summary %>%
  pivot_longer(
    cols = starts_with("tmean_"), 
    names_to = "temp_month", 
    names_prefix = "tmean_", 
    values_to = "temp_value"
  )


df_precip <- df_flr_final_summary %>%
  pivot_longer(
    cols = starts_with("ppt_"), 
    names_to = "precip_month", 
    names_prefix = "ppt_", 
    values_to = "precip_value"
  )

# Ensure that month columns have the same format for joining
df_temp <- df_temp %>%
  mutate(month = as.numeric(gsub("tmean_", "", temp_month))) 

df_precip <- df_precip %>%
  mutate(month = as.numeric(gsub("ppt_", "", precip_month)))  


df_combined <- df_temp %>%
  left_join(df_precip %>% select(id, month, precip_value), by = c("id", "month"))
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
    bal = precip_value - pet,  # Calculate the water balance
    bal = as.vector(bal)       # Convert water balance to a vector (if needed)
  )
dim(df_combined)
dim(monthly)
colnames(monthly)
str(monthly)


dat1 <- ts(monthly$bal, start = c(2018, 1), frequency = 12)  
head(dat1)

# Calculate the 1-month SPEI
spei1 <- spei(dat1, 1)
dat1_spei1 <- as.data.frame(spei1$fitted)
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



