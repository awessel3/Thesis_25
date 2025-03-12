library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(prism)

setwd("/Users/avawessel/Desktop/Thesis_25")
getwd()

# Full Dataset-----

all_flr <- read.csv("Raw Full Dataset/full_mar1125_1.csv")
no_ann <- read.csv("Raw Full Dataset/full_mar1225_2.csv")
dim(all_flr)
dim(no_ann)

match_check <- no_ann %>%
  filter(id %in% all_flr)

df_flr_final <- full_join(all_flr, no_ann)

dim(df_flr_final)

# Test Dataset ------

flowering_WVPT <- read.csv("Data/flowering_WVPT.csv")

selected_species = c("Plectritis congesta", "Collinsia grandiflora", 
                     "Plagiobothrys figuratus", "Clarkia purpurea", 
                     "Epilobium densiflorum", "Collomia grandiflora",
                     "Navarretia squarrosa")

flowering_WVPT <- flowering_WVPT %>% filter(scientific_name %in% selected_species)

### -----


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
  arrange(year,doy) %>%
  rename(species=scientific_name)

ggplot(df_flr_final, aes(x = latitude, y = doy)) + geom_point()



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

flr.clim <- extract(RS, flr.spdf,  fun=mean, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim,2); dim(flr.clim)
table(flr.clim$winter.month)
table(flr.clim$winter.year)

flr_long <- pivot_longer(flr.clim, cols=starts_with("PRISM")) %>%
  separate(col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year.prism=as.numeric(str_sub(dd, 1, 4)),
         month.prism=as.numeric(str_sub(dd, 5, 6)))%>%
  dplyr::select(-c(r1,r2,latitude.1, longitude.1, dd,field.flowering.phenology))
print(flr_long, width=Inf, n=2)
table(flr_long$month.prism)

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

flr_long <- flr_long %>% 
  dplyr::select(-c(month.prism,winter.month.prism, year.prism,class))%>% # need to remove for the pivot_wider?
  print(flr_long, width=Inf, n=2)
table(flr_long$winter.month)


# pivot back to wide
flr.dat.tmean <- pivot_wider(flr_long, names_from=tmean.month, values_from=value) %>%
  dplyr::select(-c(var))
print(flr.dat.tmean, width=Inf, n=2)


## Get precip ----
sub.ppt <- prism_archive_subset("ppt", "monthly", mon = 1:12)
RS <- pd_stack(sub.ppt) 

flr.clim <- extract(RS, flr.spdf,  fun=sum, na.rm=TRUE, sp=TRUE) 
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


flr.dat.ppt <- pivot_wider(flr_long, names_from=ppt.month, values_from=value) %>%
  dplyr::select(-c(var))
print(flr.dat.ppt, width=Inf, n=2)

flr.dat <- left_join(flr.dat.tmean, flr.dat.ppt)
print(flr.dat, width=Inf, n=2)

flr.dat <-  flr.dat %>%
  rowwise()%>%
  mutate(precip=sum(c_across(ppt_1:ppt_9))
         ,temp=mean(c_across(tmean_1:tmean_9)))
print(flr.dat, width=Inf, n=2)

dim(df_flr_final)
dim(flr.dat)

# 
#saveRDS(dat, file="data/dat_noClimate.rds") 
saveRDS(flr.dat, file="Data/df_flr_final_wClimate.rds") 

print(flr.dat, n=2, width=Inf)






