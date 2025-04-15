library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(raster)
library(sp)

install.packages("prism")
library(prism)
getwd(
)
prism_set_dl_dir("Data/prism_data")

get_prism_monthlys(
  type = "tmean", 
  year=2018:2023,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "tmin", 
  year=2018:2023,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "tmax", 
  year=2018:2023,
  mon=1:12,
  keepZip = FALSE
)
get_prism_monthlys(
  type = "ppt", 
  year=2018:2023,
  mon=1:12,
  keepZip = FALSE
)

sub.tmean <- prism_archive_subset("tmean", "monthly", mon = 1:12)
RS <- pd_stack(sub.tmean) ##raster file   # (prism_stack and ls_prism_data deprecated)

flr.test<- var4[1:4,]
flr.spdf<-   SpatialPointsDataFrame(coords=flr.test[,c('longitude','latitude')], 
                                    data=flr.test, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))

flr.clim <- extract(RS, flr.spdf,  fun=mean, na.rm=TRUE, sp=TRUE) 
flr.clim <- as.data.frame(flr.clim)
head(flr.clim)

flr_long <- pivot_longer(flr.clim, cols=starts_with("PRISM"))
print(flr_long, width=Inf, n=3)


test <- separate(flr_long, col=name, sep="_", into=c("r1","var","class","res","dd","r2")) %>%
  mutate(year=as.numeric(str_sub(dd, 1, 4)),
         month=as.numeric(str_sub(dd, 5, 6))
  )
print(test, width=Inf, n=3)





# old - extra - figuring stuff out ----
to_slice <- prism_archive_subset("tmean", "monthly", mon = 1)
boulder <- c(-105.2797, 40.0176)
p <- pd_plot_slice(to_slice, boulder)

# add a linear average and title
p + 
  stat_smooth(method="lm", se = FALSE) + 
  theme_bw() + 
  ggtitle("Average January temperature in Boulder, CO 1982-2014")


# Try to get datafrome of climate data ----
new_file<-c(1) ##change to corresponding file numbers
RS <- prism_archive_subset("tmean", "monthly", mon = 1)
df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
RS <- pd_stack(ls_prism_data()[new_file,1]) ##raster file
to_slice <- grep("_201607",RS[,1],value=T)##search through names

df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
month.df <- melt(df, c("x", "y"))
names(month.df)[1:2] <- c("lon", "lat") #rename columns


# example from https://rpubs.com/collnell/get_prism

sub.tmean <- prism_archive_subset("tmean", "monthly", mon = 1:12)
RS <- pd_stack(sub.tmean) ##raster file   # (prism_stack and ls_prism_data deprecated)
# to_slice <- grep("_201607",RS[,1],value=T)##search through names

cities<-data.frame(cities =c('Irvine','LA', 'SD', 'Palm Springs'),
                   Lat = c(33.6694649, 34.0522342, 32.7153292, 33.8302961),
                   Long = c(-117.8231107, -118.2436849, -117.1572551, -116.5452921))

# cities.spdf<-SpatialPointsDataFrame(coords=cities[,c('Long','Lat')], 
#                                     data=cities, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))
# 
# city.clim <- extract(RS, cities.spdf,  fun=mean, na.rm=TRUE, sp=TRUE)
# as.data.frame(city.clim)

