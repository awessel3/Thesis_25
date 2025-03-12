library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)


setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

#WVPT_Annual_complete <- read_rds("Data/WVPT_Annual_complete.rds")

df_flr_final_complete <- read_rds("Data/df_flr_final_complete.rds")


calculateClimate <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      observation_month = month(observed_on),
      preceding_months = list((observation_month - 1):max(1, observation_month - 3)%%12),

      preceding_temp = mean(
        unlist(lapply(preceding_months, function(m) {
          col_name <- paste0("tmean_", m)
          if (col_name %in% names(df)) df[[col_name]][cur_group_id()] else 0
        })),
        na.rm = TRUE
      ),
    
      preceding_precip = mean(
        unlist(lapply(preceding_months, function(m) {
          col_name <- paste0("ppt_", m)
          if (col_name %in% names(df)) df[[col_name]][cur_group_id()] else 0
        })),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>% 
    dplyr::select(-preceding_months, -observation_month)
}

#WVPT_climate_summary <- calculateClimate(WVPT_Annual_complete)
df_flr_final_summary <- calculateClimate(df_flr_final_complete)


saveRDS(df_flr_final_summary, file="Data/df_flr_final_summary.rds") 

write_csv(df_flr_final_summary, file="Data/df_flr_final_summary.csv") 

