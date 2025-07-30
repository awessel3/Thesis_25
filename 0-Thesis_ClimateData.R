library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(prism)

setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

############## PROCESS CLIMATE DATA ################### 

### Load raw data 

# file naming: 
# full_monthdayyr_1 = all flowering
# full_monthdayyr_2 = no annotation

## Part 1: All flowering Dataset
all_flr <- read.csv("Raw Full Dataset/full_apr725_1.csv")
length(unique(all_flr$scientific_name))

# removing subspecies
all_flr$scientific_name <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", all_flr$scientific_name)
length(unique(all_flr$scientific_name))

#Part 2: No Annotation Dataset
no_ann <- read.csv("Raw Full Dataset/full_apr725_2.csv")

# removing subspecies
no_ann$scientific_name <- gsub("^([A-Za-z]+(?:\\s+[A-Za-z]+){1}).*", "\\1", no_ann$scientific_name)
length(unique(no_ann$scientific_name))
dim(all_flr)
dim(no_ann)

# check for matches - should be 0
match_check <- no_ann %>%
  filter(id %in% all_flr)
match_check

# Join full dataset 
df_flr_nohist <- full_join(all_flr, no_ann)
unique(df_flr_nohist$scientific_name)
dim(df_flr_nohist)

## Part 3: Load in Life history strategies data
life_hist <- read.csv("Data/Flowering_WVPT_life_history.csv")
life_hist
unique(life_hist$species)
dim(life_hist)

# rename so columns match 
df_flr_nohist <- rename(df_flr_nohist, species = scientific_name)

# join life history and full dataset 
df_flr_final <- left_join(df_flr_nohist, life_hist, by = "species")
head(df_flr_final)
dim(df_flr_final)
length(unique(df_flr_final$species))


## Part 4: Set up Dates for PRISM extraction 

month <- seq(as.Date("2020-01-01"), 
             as.Date("2020-12-01"), 
             by = "1 month")
month_label <- lubridate::month(month, label = TRUE)

winter.months <- seq(as.Date("2020-09-01"), 
                     as.Date("2021-08-01"), 
                     by = "1 month")
winter.month_label <- lubridate::month(winter.months, label = TRUE)

# breaking up dates into each column, year, month, doy
df_flr_final <- df_flr_final %>%
  #replace(is.na(.),"0") %>%  # not working - why?  because it's a character field and was trying to put number
  mutate(observed_on=ymd(observed_on)
         ,year=as.numeric(year(observed_on))
         ,month=month(observed_on)
         ,doy = yday(observed_on)
  ) %>%
  arrange(year,doy) 
dim(df_flr_final)

# calculating winter year, month, doy
df_flr_final <- mutate(df_flr_final,
                 winter.year = ifelse(month>8, year+1, year)
                 ,winter.month = ifelse(month>8, month-8, month+4)
                 ,winter.day = ifelse(month>8, doy-243, doy+120)
)
head(df_flr_final ,2)
str(df_flr_final)
table(df_flr_final$winter.month)

df_flr_final <- df_flr_final  %>%  mutate(year.factor=as.factor(year), 
                             winter.year=as.factor(winter.year),
                             field.flowering.phenology = 1)

# filter to selected data timeframe 
df_flr_final <- droplevels(filter(df_flr_final , year >= 2018, year<=2024))

### Applying PRISM Climate Data

# Month count begin at September. (September = ppt_1, tmean_1, October = ppt_2, tmean_2...)

## Part 5: Download PRISM files for given year, metric

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

## Part 6: Extract temperature metrics for given coordinates in dataset

sub.tmean <- prism_archive_subset("tmean", "monthly", mon = 1:12)

RS <- pd_stack(sub.tmean) ##raster file   # (prism_stack and ls_prism_data deprecated)

flr.test <- df_flr_final
flr.spdf <- SpatialPointsDataFrame(coords=flr.test[,c('longitude','latitude')], 
                                    data=flr.test, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))

flr.clim <- raster::extract(RS, flr.spdf,  fun=mean, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim,2); dim(flr.clim)
str(flr.clim)
table(flr.clim$winter.month)
table(flr.clim$winter.year)

# Unlist if any of the columns are lists
flr.clim <- flr.clim %>% 
  mutate(across(where(is.list), ~unlist(.)))

flr_long <- pivot_longer(flr.clim, cols=starts_with("PRISM")) %>%
  separate(col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year.prism=as.numeric(str_sub(dd, 1, 4)),
         month.prism=as.numeric(str_sub(dd, 5, 6)))%>%
  dplyr::select(-c(r1,r2,latitude.1, longitude.1, dd,field.flowering.phenology))
print(flr_long, width=Inf, n=2)
table(flr_long$month.prism)
str(flr_long)

# Align prism months with months from dataset accounting for same winter-year. 
flr_long <- flr_long %>%  
  # filter(class != "provisional")%>%
  mutate(winter.year.prism = ifelse(month.prism>8, year.prism+1, year.prism)
         ,winter.month.prism = ifelse(month.prism>8, month.prism-8, month.prism+4)
         ,tmean.month = paste(var, winter.month.prism, sep="_")
  )%>%
  filter(winter.year.prism == winter.year) # filter to only the prism data from the same winter-year as the observation

print(flr_long, width=Inf, n=2)
table(flr_long$winter.month.prism)
table(flr_long$winter.year.prism)
str(flr_long)

flr_long <- flr_long %>% 
  dplyr::select(-c(month.prism,winter.month.prism, year.prism,class))%>% # need to remove for the pivot_wider?
  print(flr_long, width=Inf, n=2)
table(flr_long$winter.month)
str(flr_long)

#addition - currently temp 
flr_long <- flr_long %>%
  group_by(id, observed_on, latitude, longitude, species, life_history, year, 
           month, doy, winter.year, winter.month, winter.day, year.factor, var, 
           res, winter.year.prism, tmean.month) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

flr.dat.tmean <- pivot_wider(flr_long, names_from = tmean.month, values_from = value) %>%
  dplyr::select(-c(var))

# pivot back to wide
flr.dat.tmean <- pivot_wider(flr_long, names_from=tmean.month, values_from=value) %>%
  dplyr::select(-c(var))


print(flr.dat.tmean, width=Inf, n=2)
str(flr.dat.tmean)


## Part 7: Extract precipitation metrics for given coordinates in dataset

sub.ppt <- prism_archive_subset("ppt", "monthly", mon = 1:12)
RS <- pd_stack(sub.ppt) 

flr.clim <- raster::extract(RS, flr.spdf,  fun=sum, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim,2)

flr_long <- pivot_longer(flr.clim, cols=starts_with("PRISM")) %>%
  separate(col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year.prism=as.numeric(str_sub(dd, 1, 4)),
         month.prism=as.numeric(str_sub(dd, 5, 6)))%>%
  dplyr::select(-c(r1,r2,latitude.1, longitude.1, dd,field.flowering.phenology))
print(flr_long, width=Inf, n=2)

flr_long <- flr_long %>%  
  # filter(class != "provisional")%>%
  mutate(winter.year.prism = ifelse(month.prism>8, year.prism+1, year.prism)
         ,winter.month.prism = ifelse(month.prism>8, month.prism-8, month.prism+4)
         ,ppt.month = paste(var, winter.month.prism, sep="_")
  )%>%
  filter(winter.year.prism == winter.year) # filter to only the prism data from the same winter-year as the observation
print(flr_long, width=Inf, n=2)

flr_long <- flr_long %>% 
  dplyr::select(-c(month.prism,winter.month.prism, year.prism,class)) # need to remove for the pivot_wider?
print(flr_long, width=Inf, n=2)
table(flr_long$winter.month)

# new addition - temporary 
flr_long <- flr_long %>%
  group_by(id, observed_on, latitude, longitude, species, life_history, year, 
           month, doy, winter.year, winter.month, winter.day, year.factor, var, 
           res, winter.year.prism, ppt.month) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

flr.dat.ppt <- pivot_wider(flr_long, names_from = ppt.month, values_from = value) %>%
  dplyr::select(-c(var))


flr.dat.ppt <- pivot_wider(flr_long, names_from=ppt.month, values_from=value) %>%
  dplyr::select(-c(var))


print(flr.dat.ppt, width=Inf, n=2)
str(flr.dat.ppt)

## Part 8: Join temperature and precipitation data

flr.dat <- left_join(flr.dat.tmean, flr.dat.ppt)
str(flr.dat)
dim(flr.dat)
print(flr.dat, width=Inf, n=2)

flr.dat <-  flr.dat %>%
  rowwise()%>%
  mutate(precip=sum(c_across(ppt_1:ppt_9))
         ,temp=mean(c_across(tmean_1:tmean_9))) %>% 
  mutate(across(where(is.list), ~unlist(.)))
print(flr.dat, width=Inf, n=2)
str(flr.dat)

dim(df_flr_final)
dim(flr.dat)
str(flr.dat)

flr.dat <- na.omit(flr.dat)
dim(flr.dat)
str(flr.dat)

# test plot 

ggplot(flr.dat, aes(x = temp, y = doy)) + geom_point()

## Part 9: Saving data

#saveRDS(dat, file="data/dat_noClimate.rds") 
saveRDS(flr.dat, file="Data/df_flr_final_wClimate.rds") 

print(flr.dat, n=2, width=Inf)




