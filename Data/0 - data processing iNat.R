library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(raster)
library(sp)
library(prism)

prism_set_dl_dir("data/prism_data")

# Useful for plotting dates 
month <- seq(as.Date("2020-01-01"), 
             as.Date("2020-12-01"), 
             by = "1 month")
month_label <- lubridate::month(month, label = TRUE)

winter.months <- seq(as.Date("2020-09-01"), 
                     as.Date("2021-08-01"), 
                     by = "1 month")
winter.month_label <- lubridate::month(winter.months, label = TRUE)


# WVPT species
sp.all <- read.csv("Data/plant_species_on_iNat.csv")
sp <- read.csv("Data/sp.list.WVPT.csv")
# sp = species we want to use because they've been scored for phenology on iNat
sp <- left_join(sp, sp.all, by=c("species"="scientificName"))
head(sp)

# get IDs for querying iNat
# write.csv(t(sp$id), "sp_IDs_WVPT.csv",  row.names = FALSE,  quote=FALSE)

## read data ----

#var1 <- read.csv("Data/veg_wvpt.csv")
var4 <- read.csv("data/flowering_WVPT.csv") 

# phases <- c("Green","Buds","Flowers","Fruits")
# 
# head(var1)
# head(var2)
# head(var3)
head(var4,3); dim(var4)
str(var4)
head(sp,2); dim(sp)

var4$taxon_id
sp$id

# dat <- rbind(var1,var2,var3,var4) %>%
# Process raw data ----
dat <- var4 %>%
  filter(taxon_id %in% sp$id)%>%
  # pivot_wider(names_from=phenol, values_from = phenol) %>%
  replace(is.na(.),"0") %>%  # not working - why?  because it's a character field and was trying to put number
  mutate(flowering = 1
         ,observed_on=mdy(observed_on)
         # ,observed_on=ymd(observed_on)
         ,year=as.numeric(year(observed_on))
         ,month=month(observed_on)
         ,doy = yday(observed_on)
) %>%
arrange(year,doy) %>%
rename(species=scientific_name)

head(dat,2); dim(dat)
head(var4,2); dim(var4)

N.yearly <- dat %>% group_by(year) %>%
    summarize(N = length(species)) %>%
    arrange(desc(year))
print(N.yearly, n=Inf)
  
  pdf("Figures/N_over_time_allSpecies.pdf",width=6,height=6)
  ggplot(N.yearly, aes(y=N, x=year))+
    geom_bar(stat = "identity")+
    labs(y="Number of observations")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  dev.off()  
  
  
## take out early years with few observations----  
  dat <- droplevels(filter(dat, year > 2018, year<2024))
  head(dat); dim(dat)
  
  table(dat$year)

  # so here in dat we have observations from 2019 - 2023, for our subset of focal species in WVPT
  
  # subset to one species
  collin <- dat %>% filter(species=="Collinsia grandiflora")

  sp1 <- dat %>% filter(species=="Plectritis congesta")
  sp1 <- filter(dat, species=="Plectritis congesta")
  
# Get winter months ----
dat <- mutate(dat,
                  winter.year = ifelse(month>8, year+1, year)
                  ,winter.month = ifelse(month>8, month-8, month+4)
                  ,winter.day = ifelse(month>8, doy-243, doy+120)
)
head(dat,2)
str(dat)
table(dat$winter.month)


dat <- dat %>%  mutate(year.factor=as.factor(year), 
                       winter.year=as.factor(winter.year))

# link to PRISM data ----

# get PRISM data
sub.tmean <- prism_archive_subset("tmean", "monthly", mon = 1:12)
RS <- pd_stack(sub.tmean) ##raster file   # (prism_stack and ls_prism_data deprecated)

flr.test<- dat
# flr.test<- dat[1:10,]
table(flr.test$winter.month)
flr.spdf<-   SpatialPointsDataFrame(coords=flr.test[,c('longitude','latitude')], 
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

# 
saveRDS(dat, file="data/dat_noClimate.rds") 
saveRDS(flr.dat, file="data/dat_withClimate.rds") 

print(flr.dat, n=2, width=Inf)

