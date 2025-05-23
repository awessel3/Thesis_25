library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(prism)

# setwd("/Users/avawessel/Desktop/Thesis_25")
# getwd()



# Full Dataset-----

all_flr <- read.csv("Raw Full Dataset/full_apr725_1.csv")
length(unique(all_flr$scientific_name))
all_flr$scientific_name <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", all_flr$scientific_name)
length(unique(all_flr$scientific_name))
no_ann <- read.csv("Raw Full Dataset/full_apr725_2.csv")
no_ann$scientific_name <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", no_ann$scientific_name)
length(unique(no_ann$scientific_name))
dim(all_flr)
dim(no_ann)

match_check <- no_ann %>%
  filter(id %in% all_flr)
match_check

df_flr_nohist <- full_join(all_flr, no_ann)
unique(df_flr_nohist$scientific_name)
dim(df_flr_nohist)

life_hist <- read.csv("Data/Flowering_WVPT_life_history.csv")
life_hist
unique(life_hist$species)
dim(life_hist)

df_flr_nohist <- rename(df_flr_nohist, species = scientific_name)

df_flr_final <- left_join(df_flr_nohist, life_hist, by = "species")
head(df_flr_final)
dim(df_flr_final)
length(unique(df_flr_final$species))

# (Test Dataset) ------

# flowering_WVPT <- read.csv("Data/flowering_WVPT.csv")
# flowering_WVPT
# 
# selected_species = c("Plectritis congesta", "Collinsia grandiflora", 
#                      "Plagiobothrys figuratus", "Clarkia purpurea", 
#                      "Epilobium densiflorum", "Collomia grandiflora",
#                      "Navarretia squarrosa")
# 
# flowering_WVPT <- flowering_WVPT %>% filter(scientific_name %in% selected_species)

# -----


month <- seq(as.Date("2020-01-01"), 
             as.Date("2020-12-01"), 
             by = "1 month")
month_label <- lubridate::month(month, label = TRUE)

winter.months <- seq(as.Date("2020-09-01"), 
                     as.Date("2021-08-01"), 
                     by = "1 month")
winter.month_label <- lubridate::month(winter.months, label = TRUE)




df_flr_final <- df_flr_final %>%
  #replace(is.na(.),"0") %>%  # not working - why?  because it's a character field and was trying to put number
  mutate(observed_on=ymd(observed_on)
         ,year=as.numeric(year(observed_on))
         ,month=month(observed_on)
         ,doy = yday(observed_on)
  ) %>%
  arrange(year,doy) 
dim(df_flr_final)
head(df_flr_final)

ggplot(df_flr_final, aes(x = latitude, y = doy)) + geom_point()



df_flr_final <- mutate(df_flr_final,
                 winter.year = as.numeric(ifelse(month>8, year+1, year))
                 ,winter.month = ifelse(month>8, month-8, month+4)
                 ,winter.day = ifelse(month>8, doy-243, doy+120)
)
head(df_flr_final ,2)
str(df_flr_final)
table(df_flr_final$winter.month)

# df_flr_final <- df_flr_final  %>%  mutate(year.factor=as.factor(year), 
#                              winter.year=as.factor(winter.year),
#                              field.flowering.phenology = 1)

df_flr_final <- droplevels(filter(df_flr_final , year >= 2018, year<=2024))

#explore 
ggplot(df_flr_final, aes(x = latitude, y = doy)) + geom_point()

#printing out weird, processing error occuring during these steps? 
####


# PRISM (DOES NOT CONTAIN CANADA) -------

prism_set_dl_dir("Data/prism_data")

get_prism_monthlys(
  type = "tmean", 
  year=2017:2024,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "tmin", 
  year=2017:2024,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "tmax", 
  year=2017:2024,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "ppt", 
  year=2017:2024,
  mon=1:12,
  keepZip = FALSE
)

# get PRISM data
sub.tmean <- prism_archive_subset("tmean", "monthly", mon = 1:12)

RS <- pd_stack(sub.tmean) ##raster file   # (prism_stack and ls_prism_data deprecated)

flr.test <- df_flr_final
flr.spdf <-   SpatialPointsDataFrame(coords=flr.test[,c('longitude','latitude')], 
                                    data=flr.test, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))

flr.clim <- raster::extract(RS, flr.spdf,  fun=mean, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim,2); dim(flr.clim)
table(flr.clim$winter.month)
table(flr.clim$winter.year)

flr_long <- flr.clim %>% pivot_longer(cols=starts_with("PRISM")) %>%
  separate(col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year.prism=as.numeric(str_sub(dd, 1, 4)),
         month.prism=as.numeric(str_sub(dd, 5, 6)),
         winter.year = as.numeric(winter.year)) %>%
  dplyr::select(-c(r1,r2,latitude.1, longitude.1, dd)) # field.flowering.phenology

print(flr_long, width=Inf, n=2); dim(flr_long)
table(flr_long$month.prism)
table(flr_long$var)

flr_long <- flr_long %>%  
  # filter(class != "provisional")%>%
  mutate(winter.year.prism = ifelse(month.prism>8, year.prism+1, year.prism)
         ,winter.month.prism = ifelse(month.prism>8, month.prism-8, month.prism+4)
         ,tmean.month = paste(var, winter.month.prism, sep="_")
  )%>%
  filter(winter.year.prism == winter.year, winter.year < 2025) # filter to only the prism data from the same winter-year as the observation

print(flr_long, width=Inf, n=4)
table(flr_long$winter.month.prism)
table(flr_long$winter.year.prism)

# take out the columns we don't need anymore
flr_long <- flr_long %>% 
  dplyr::select(-c(res,var,class, month.prism,winter.month.prism,year.prism)) # need to remove for the pivot_wider?

print(flr_long, width=Inf, n=2); dim(flr_long)
table(flr_long$winter.month)


print(filter(flr_long, id == 9375178), width=Inf, n=Inf)

# pivot back to wide
# ! need to add ID Cols?  Was getting lists for values (multiple months for each obs)  ----
#OLD:  flr.dat.tmean <- pivot_wider(flr_long, names_from=tmean.month, values_from=value) 
flr.dat.tmean <- flr_long %>% pivot_wider(names_from=tmean.month, values_from=value) #%>%

# flr.dat.tmean <- flr_long %>% pivot_wider(id_cols=c(id,observed_on,latitude, longitude,species,life_history, year, month, doy, winter.year, winter.month, winter.day, year.factor, winter.month.prism), 
#                                           names_from=tmean.month, values_from=value) #%>%
# dplyr::select(-c(var))
print(flr.dat.tmean, width=Inf, n=3); dim(flr.dat.tmean)

# diagnose:
# test <- flr_long %>%
#   count(across(c(
#     id, observed_on, latitude, longitude, species,
#     life_history, year, month, doy, winter.year,
#     winter.month, winter.day, year.factor, tmean.month
#   ))) %>%
#   filter(n > 1)

print(test, width=Inf, n=5); dim(test)
print(filter(flr_long, id == 239248340), width=Inf, n=Inf)



  

## Get precip ----
sub.ppt <- prism_archive_subset("ppt", "monthly", mon = 1:12)
RS <- pd_stack(sub.ppt) 

flr.clim <- raster::extract(RS, flr.spdf,  fun=sum, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim,2)

flr_long <- pivot_longer(flr.clim, cols=starts_with("PRISM")) %>%
  separate(col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year.prism=as.numeric(str_sub(dd, 1, 4)),
         month.prism=as.numeric(str_sub(dd, 5, 6)))%>%
  dplyr::select(-c(r1,r2,latitude.1, longitude.1, dd)) # field.flowering.phenology
print(flr_long, width=Inf, n=2)

flr_long <- flr_long %>%  
  # filter(class != "provisional")%>%
  mutate(winter.year.prism = ifelse(month.prism>8, year.prism+1, year.prism)
         ,winter.month.prism = ifelse(month.prism>8, month.prism-8, month.prism+4)
         ,ppt.month = paste(var, winter.month.prism, sep="_")
  )%>%
  filter(winter.year.prism == winter.year, winter.year < 2025) # filter to only the prism data from # filter(winter.year.prism == winter.year) # filter to only the prism data from the same winter-year as the observation
print(flr_long, width=Inf, n=2)

flr_long <- flr_long %>% 
  dplyr::select(-c(month.prism,winter.month.prism, year.prism,class)) # need to remove for the pivot_wider?
print(flr_long, width=Inf, n=2)
table(flr_long$winter.month)


flr.dat.ppt <- pivot_wider(flr_long, names_from=ppt.month, values_from=value) %>%
  dplyr::select(-c(var))

print(flr.dat.ppt, width=Inf, n=2)
print(flr.dat.tmean, width=Inf, n=2)

# Join temp and precip
flr.dat <- left_join(flr.dat.tmean, flr.dat.ppt)
print(flr.dat, width=Inf, n=2)

flr.dat <-  flr.dat %>%
  rowwise()%>%
  mutate(precip=sum(c_across(ppt_1:ppt_9))
         ,temp=mean(c_across(tmean_1:tmean_9)))
print(flr.dat, width=Inf, n=2)

dim(df_flr_final)
dim(flr.dat)

flr.dat <- na.omit(flr.dat)
ggplot(flr.dat, aes(x = temp, y = doy)) + geom_point()
ggplot(flr.dat, aes(x = temp, y = precip)) + geom_point()

## take out weird outlier temp ---- 
flr.dat <- flr.dat %>% filter(temp > 2)
ggplot(flr.dat, aes(x = temp, y = precip)) + geom_point()

# 
#saveRDS(dat, file="data/dat_noClimate.rds") 
saveRDS(flr.dat, file="Data/df_flr_final_wClimate_JD.rds") 



