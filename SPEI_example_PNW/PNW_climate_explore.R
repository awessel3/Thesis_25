# Exploring precip and temp data from climate stations

library(brms)
library(lme4)
library(quantreg)

library(ggplot2)
library(RColorBrewer) 
library(gridExtra) # for grid.arrange 

library(tidyverse)
library(prettydoc)
library(tidybayes)
library(dplyr)
library(reshape2)

library(GGally)
library(SPEI)
library(viridis)
library(ggrepel)
library(patchwork)


## ggplot controls
theme_set(theme_minimal())

all.raw.data <- read.table("all.stations.csv", sep=',', header=T )
table(all.raw.data$NAME)


# Climate data ----

head(all.raw.data)
cdata0 <- all.raw.data %>% 
  mutate(year = strptime(DATE, format="%m/%d/%y")$year + 1900
         , month = strptime(DATE, format="%m/%d/%y")$mon + 1
         , day.year = strptime(DATE, format="%m/%d/%y")$yday + 1
  ) %>%
  drop_na(PRCP) %>%
  drop_na(TMAX) %>%
  drop_na(TMIN) 
head(cdata0)
table(cdata0$NAME)

## Get winter & spring months ----
cdata0 <- mutate(cdata0,
                 tmean =  (TMAX+TMIN) / 2
                 # tmean = as.numeric(ifelse(!is.na(TMAX), (TMAX+TMIN) / 2,  "NA"))
                 ,winter.year = ifelse(month>8, year+1, year)
                 ,winter.month = ifelse(month>8, month-8, month+4)
                 ,winter.day = ifelse(month>8, day.year-243, day.year+120)
                 #,PRCP = PRCP / 10 # convert mm to cm LATER because Thornthwaite equation needs mm
) %>%
  filter(winter.year > 1979) %>%
  # filter(winter.year > 1999) %>%
  drop_na(tmean)

head(cdata0); dim(cdata0)


cdata <- cdata0
# get cumulative Precip values
cdata <- arrange(cdata, STATION, winter.year, winter.day)

cdata <- cdata %>%  group_by(STATION, winter.year) %>%
  mutate(cum_PRCP = cumsum(PRCP)) %>%
  # ungroup()%>%
  as.data.frame()

head(cdata)

dim(cdata); table(cdata$NAME)
cdata$NAME <- recode_factor(cdata$NAME, `EUGENE MAHLON SWEET FIELD, OR US` = "Eugene",
                            `MEDFORD INTERNATIONAL AIRPORT, OR US` = "Medford",
                            `OLYMPIA AIRPORT, WA US` = "Olympia",
                            `PENDLETON E OR REGIONAL AIRPORT, OR US` = "Pendleton", 
                            `PORTLAND INTERNATIONAL AIRPORT, OR US` = "Portland",
                            `HJA` = "HJ Andrews",
                            `SEATTLE TACOMA AIRPORT, WA US` = "Seattle", 
                            `SEXTON SUMMIT, OR US` = "Sexton", 
                            `WHIDBEY ISLAND NAS, WA US` = "Whidbey", 
                            `SALEM AIRPORT MCNARY FIELD, OR US` = "Salem")

dim(cdata); table(cdata$NAME)
cdata <- droplevels(cdata)
dim(cdata); table(cdata$NAME)
cdata <- rename(cdata, station=NAME) %>%
  select(-STATION)
dim(cdata); table(cdata$station)

cdata$station <- fct_relevel(cdata$station,  "Medford" , "Sexton","Eugene"  ,   "Salem", "Portland"  , "Olympia", "Seattle" ,"Whidbey" )

# Get monthly summaries
monthly0 <- cdata %>% group_by(station, winter.year, winter.month) %>%
  summarise( tmean = mean(tmean, na.rm=TRUE)
             , precip = sum(PRCP)
             , lat=mean(LATITUDE)
  ) %>%
  #filter(winter.year>1938) %>%
  as.data.frame()
head(monthly0)

stations.lat <- group_by(cdata, station) %>%
  summarise(station=unique(station),
            lat = median(LATITUDE))%>%
  arrange(lat)

## Calculate SPEI ----

# calculate a few things needed for SPEI
# PET - potential evapotranspiration
# bal - water balance

# Start with first station
monthly <-filter(monthly0, station==stations.lat$station[1]) %>%
  mutate(station=unique(station),
         pet = thornthwaite(tmean, stations.lat$lat[1]),
         bal = precip - pet,
         bal = as.vector(bal)
  )

## Calculate SPEI
dat1 <- ts(monthly[,2:8], end=c(2020, 12), frequency=12)  # make monthly data a time series object (ts)
dat1 <- as.data.frame(dat1)
# there are different "'"scales" over which you calculate SPEI; i have mostly used 1, but calculate here 3-month and 6-month too.  i think of it like a moving window, but i haven't looked into these enough
spei1 <- spei(dat1[,'bal'], 1) # calculate SPEI from water balance
dat1$spei1 <- as.vector(spei1$fitted)
spei3 <- spei(dat1[,'bal'], 3) 
dat1$spei3 <- as.vector(spei3$fitted)
spei6 <- spei(dat1[,'bal'], 6)
dat1$spei6 <- as.vector(spei6$fitted)
monthly <- bind_cols(select(monthly,station), dat1)
head(monthly)

# Do same for the other stations
for(i in 2:dim(stations.lat)[1]) {
  s1 <- filter(monthly0, station==stations.lat$station[i])
  s2 <- mutate(s1, pet = thornthwaite(tmean, stations.lat$lat[i]))
  s2$bal <- s2$precip - s2$pet
  s2$bal <- as.vector(s2$bal)
  s2 <- ts(s2[,2:8], end=c(2020, 12), frequency=12)
  s2 <- as.data.frame(s2)
  spei1 <- spei(s2[,'bal'], 1)
  s2$spei1 <- as.vector(spei1$fitted)
  spei3 <- spei(s2[,'bal'], 3)
  s2$spei3 <- as.vector(spei3$fitted)
  spei6 <- spei(s2[,'bal'], 6)
  s2$spei6 <- as.vector(spei6$fitted)
  s2 <- bind_cols(select(s1,station), s2)
  
  monthly <- rbind(monthly, s2)
}
head(monthly); dim(monthly)
table(monthly$station); names(monthly)

# ! convert PRCP to cm ----
cdata <- mutate(cdata, 
                PRCP=PRCP/10 # convert to cm
                # ,cum_PRCP = cum_PRCP/10)
)

monthly <- mutate(monthly, 
                  precip=precip/10 # convert to cm
)


# get annual climate ----

annual.climate <- group_by(cdata, station, winter.year) %>%
  summarise(precip = sum(PRCP, na.rm=TRUE) 
            ,mean.temp = signif(mean(tmean, na.rm=TRUE), 4)) %>%
  filter(winter.year < 2025) %>%# not full 2025 data so is biased downward
  as.data.frame()


annual.spei <- group_by(monthly, station, winter.year) %>%
  summarize(spei1 = mean(spei1, na.rm=TRUE)
            ,spei3 = mean(spei3, na.rm=TRUE)
            ,spei6= mean(spei6, na.rm=TRUE))%>%
  filter(winter.year < 2025) # not full 2025 data so is biased downward
  head(annual.spei)
head(annual.climate)

annual.climate <- left_join(annual.climate, annual.spei)


# ! take out Eugene < 2000 ----
cdata <- cdata %>%
  filter(!(station == "Eugene" & winter.year < 2000))

annual.climate <- annual.climate %>%
  filter(!(station == "Eugene" & winter.year < 2000))


# all processing up until here ----

# plot monthly----

month <- seq(as.Date("2020-01-01"), 
             as.Date("2020-12-01"), 
             by = "1 month")
month_label <- lubridate::month(month, label = TRUE)

winter.months <- seq(as.Date("2020-09-01"), 
                     as.Date("2021-08-01"), 
                     by = "1 month")
winter.month_label <- lubridate::month(winter.months, label = TRUE)

monthly.melt <- pivot_longer(monthly, cols=c(tmean:spei6), names_to = "variable",values_to = "value")
head(monthly.melt)

monthly.melt$winter.month <- as.factor(monthly.melt$winter.month)

temp.max.grand<-max(monthly$tmean)
ppt.max.grand<-max(monthly$precip)

t.m <- ggplot(monthly, aes(x=winter.month, y=tmean, color=station)) +
  geom_smooth(se=FALSE)+  # 
  scale_x_continuous(breaks = 1:12, labels = winter.month_label)+
  ylab("Mean Temperature (C)") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_minimal() + labs(title = "", subtitle = "", caption = "", tag = "A") 

p.m <- ggplot(monthly, aes(x=winter.month, y=precip, color=station)) +
  geom_smooth(se=FALSE)+  # 
  scale_x_continuous(breaks = 1:12, labels = winter.month_label) +
  ylab("Mean Precipitation (cm)") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_minimal() + labs(title = "", subtitle = "", caption = "Figure 1.", tag = "B") +
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))

t.m/p.m
ggsave("Figures/means.pdf", width=10, height=6)

# ggplot(filter(monthly, winter.year > 2000), aes(x=tmean, y=precip, color=as.factor(winter.month))) +
ggplot(monthly, aes(x=tmean, y=precip, color=as.factor(winter.month))) +
  # geom_smooth()+  # 
  # scale_x_continuous(breaks = 1:12, labels = winter.month_label) +
  geom_point(alpha=.4, shape=16)+
  facet_wrap(~station)+
  ylab("Mean Precipitation (cm)") + xlab("Mean Temp (C)") 
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   theme_minimal() + labs(title = "", subtitle = "", caption = "Figure 1.", tag = "B") +
# theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
ggsave("Figures/climate_space_monthly.pdf", width=10, height=10)

## Climate space ----

pw <- ggplot(annual.climate, aes(x=mean.temp, y=precip, color=station)) +
  geom_point(alpha=.4)+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2))+
  labs(x='Mean temperature (C)', y='Annual precipitation (cm)',title = "", subtitle = "", caption = "Figure 2. Climate space; each point is a year", tag = "") +
  theme_minimal() + theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))

pw
ggsave("Figures/climate_space_annual.pdf", width=10, height=10)


## Trends over time

# ! add lables by stations -----
p1<-ggplot(data=annual.climate, aes(x=winter.year, y=precip, color=station)) +
  geom_point(alpha=.5)+
  geom_hline(yintercept=mean(annual.climate$precip), linetype='dashed')+
  geom_smooth(method=lm,se=FALSE)+
  labs(y="Annual precipitation (cm)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


t1<-ggplot(data=annual.climate, aes(x=winter.year, y=mean.temp, color=station)) +
  geom_point(alpha=.5)+ 
  geom_smooth(method=lm,se=FALSE)+
  geom_hline(yintercept=mean(annual.climate$mean.temp), linetype='dashed')+
  labs(y="Mean temperature (C)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal() + 
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

s1<-ggplot(data=annual.climate, aes(x=winter.year, y=spei1, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm,se=FALSE)+
  geom_hline(yintercept=mean(annual.climate$spei1), linetype='dashed') +
  labs(y="SPEI", x = "Year", title = "", subtitle = "",  tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 

s12<-ggplot(data=annual.climate, aes(x=winter.year, y=spei12, color=station)) +
  geom_point(alpha=.5)+
  ylab("Drought index (SPEI12)") + xlab("Year")+
  geom_smooth()

p1/t1/s1+ plot_layout(guides = "collect")
# p1+t1+s1 + plot_layout(guides = "collect")
ggsave("Figures/climate trends.pdf", width=8, height=12)


# add month name to cdata
month.names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
month.names.winter <- c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug")
#month.names2 <- c(1="Jan",2="Feb",3="Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
lookup <- tibble(month.name = month.names, month = 1:12)
cdata <- inner_join(cdata,lookup) 
cdata$month.name <- as.factor(cdata$month.name)
cdata$month.name <- fct_relevel(cdata$month.name,  "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

## Have the dry seasons been getting longer?

### Trends in monthly precip

# Calculate how long the growing season is in a year
# thresh 

winter.month.names <- c(
  "1" = "September", "2" = "October", "3" = "November",
  "4" = "December", "5" = "January", "6" = "February",
  "7" = "March", "8" = "April", "9" = "May",
  "10" = "June", "11" = "July", "12" = "August"
)

for(i in 1:dim(stations.lat)[1]) {
  #i<-1
  site.data <- filter(monthly.melt, station==stations.lat$station[i])
  
  # precip trends 
  ggplot(filter(site.data, variable=="precip"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm, se = FALSE)+
    labs(y="Precipitation (cm)", x = "", title = "Precip", subtitle = unique(site.data$station), caption = "", tag = "") +
    geom_point(alpha=.4)+
    scale_color_discrete(labels = winter.month.names)
  ggsave(paste("Figures/",unique(site.data$station),".monthly precip trends.pdf", sep=''), width=8, height=6)
  
  ggplot(filter(site.data, variable=="precip"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm)+
    labs(y="Precipitation (cm)", x = "", title = "Precip", subtitle = unique(site.data$station), caption = "", tag = "") +
    geom_point()+
    scale_color_discrete(labels = winter.month.names)+
    facet_wrap(~winter.month, ncol = 3, labeller = labeller(winter.month = winter.month.names)) 
  # geom_text_repel(data=filter(monthly.cal, year=="2019"), aes(x=2023, y=precip, label=month.names),size=4,segment.alpha=0, box.padding=0)+
  # expand_limits(x = 2030)
  ggsave(paste("Figures/",unique(site.data$station),".monthly precip separate.pdf", sep=''), width=8, height=6)
  
  # Temp trends 
  ggplot(filter(site.data, variable=="tmean"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm, se = FALSE)+
    labs(y="Temperature (Mean C)", x = "", title = "Temperature", subtitle = unique(site.data$station), caption = "", tag = "") +
    geom_point(alpha=.4)+
    scale_color_discrete(labels = winter.month.names)
  ggsave(paste("Figures/",unique(site.data$station),".monthly temp trends.pdf", sep=''), width=8, height=6)
  
  ggplot(filter(site.data, variable=="tmean"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm)+
    labs(y="Temperature (Mean C)", x = "", title = "Temperature", subtitle = unique(site.data$station), caption = "", tag = "") +
    geom_point()+
    scale_color_discrete(labels = winter.month.names)+
    facet_wrap(~winter.month, ncol = 3, labeller = labeller(winter.month = winter.month.names)) 
  # geom_text_repel(data=filter(monthly.cal, year=="2019"), aes(x=2023, y=precip, label=month.names),size=4,segment.alpha=0, box.padding=0)+
  # expand_limits(x = 2030)
  ggsave(paste("Figures/",unique(site.data$station),".monthly temp separate.pdf", sep=''), width=8, height=6)
  
  # SPEI trends
  ggplot(filter(site.data, variable=="spei1"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm, se = FALSE)+
    labs(y="Standardised Precipitation-Evapotranspiration Index \n(SPEI)", x = "", title = "SPEI", subtitle = unique(site.data$station), caption = "", tag = "") +  geom_point(alpha=.4)+
    geom_hline(yintercept=0, linetype='dashed')+
    scale_color_discrete(labels = winter.month.names)
  ggsave(paste("Figures/",unique(site.data$station),".monthly SPEI trends.pdf", sep=''), width=8, height=6)
  
  ggplot(filter(site.data, variable=="spei1"), aes(x=winter.year, y=value, color=winter.month)) +
    geom_smooth(method=lm)+
    labs(y="Standardised Precipitation-Evapotranspiration Index \n(SPEI)", x = "", title = "SPEI", subtitle = unique(site.data$station), caption = "", tag = "") +
    geom_point()+
    geom_hline(yintercept=0, linetype='dashed')+
    scale_color_discrete(labels = winter.month.names)+
    facet_wrap(~winter.month, ncol = 3, labeller = labeller(winter.month = winter.month.names)) 
  # geom_text_repel(data=filter(monthly.cal, year=="2019"), aes(x=2023, y=SPEI, label=month.names),size=4,segment.alpha=0, box.padding=0)+
  # expand_limits(x = 2030)
  ggsave(paste("Figures/",unique(site.data$station),".monthly SPEI separate.pdf", sep=''), width=8, height=6)
  
}

ggplot(filter(monthly.melt, variable=="spei1"), aes(x=winter.year, y=value, color=station)) +
  geom_point(alpha=.2)+
  geom_smooth(method=lm, se=FALSE)+
  labs(y="Standardised Precipitation-Evapotranspiration Index \n(SPEI)", x = "", title = "SPEI", subtitle = "all stations", caption = "", tag = "") +
  geom_hline(yintercept=0, linetype='dashed')+
  scale_color_discrete(labels = winter.month.names)+
  facet_wrap(~winter.month, ncol = 3, labeller = labeller(winter.month = winter.month.names)) 
ggsave("Figures/monthly SPEI separate.pdf", width=8, height=8)


### Trends in precip from spring to fall

seasonal.climate <- group_by(filter(cdata, station != "HJ Andrews"), station, winter.year) %>%
  summarise(precip38 = sum(PRCP[month > 2 & month<7], na.rm=TRUE),
            temp38 = mean(tmean[month > 2 & month<7], na.rm=TRUE),
            precip_winter = sum(PRCP[month > 9 | month<3], na.rm=TRUE),
            temp_winter = mean(tmean[month > 9 | month<3], na.rm=TRUE),
            precip49 = sum(PRCP[month > 3 & month<10], na.rm=TRUE),
            temp49 = mean(tmean[month > 3 & month<10], na.rm=TRUE)) %>%
  as.data.frame()


#seasonal.climate$station <- fct_relevel(seasonal.climate$station,  "Medford" ,"Eugene"  ,  "Salem", "Portland"  ,"Pendleton", "Olympia"  )
p1 <- ggplot(seasonal.climate, aes(x=winter.year, y=precip38, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+
  labs(y="Seasonal Precipitation (cm)", x = "Year", title = "March-June Total Precip", subtitle = "", caption = "", tag = "") +
  theme_minimal() + 
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 

p2 <- ggplot(seasonal.climate, aes(x=winter.year, y=precip49, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+
  labs(y="Seasonal Precipitation (cm)", x = "Year", title = "April-September Total Precip", subtitle = "", caption = "", tag = "") +
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 

p3 <- ggplot(seasonal.climate, aes(x=winter.year, y=precip_winter, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+
  labs(y="Seasonal Precipitation (cm)", x = "Year", title = "Oct-Feb Total Precip", subtitle = "", caption = "", tag = "") +
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) 

p1+p2 +p3+ plot_layout(guides = "collect")
ggsave("Figures/Seasonal Precip trends.pdf", width=14, height=6)


p1 <- ggplot(seasonal.climate, aes(x=winter.year, y=temp38, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+ 
  ylim(c(10,20))+
  ylab("Seasonal Mean Temp (C)") + xlab("") + ggtitle("March-June Temp")
p2 <- ggplot(seasonal.climate, aes(x=winter.year, y=temp49, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+
  ylab("Seasonal Mean Temp (C)") + xlab("") + ggtitle("April-September Temp")
p3 <- ggplot(seasonal.climate, aes(x=winter.year, y=temp_winter, color=station)) +
  geom_smooth(method=lm)+  # 
  geom_point(alpha=.2)+
  ylab("Seasonal Mean Temp (C)") + xlab("") + ggtitle("Oct-Feb Temp")

p1+p2+p3+ plot_layout(guides = "collect")
ggsave("Figures/Seasonal Temp trends.pdf", width=14, height=6)

# Has cloud cover been changing? -----

# ACSC = Average cloudiness sunrise to sunset from 30-second ceilometer data (percent) 
# ACSH = Average cloudiness sunrise to sunset from manual observations (percent)
# ACMH = Average cloudiness midnight to midnight from manual observations (percent)

monthly.clouds <- cdata %>% group_by(station, year, month) %>%
  summarise( acsh = mean(ACSH, na.rm=TRUE),
             precip = sum(PRCP),
             month.name = unique(month.name)
  ) %>%
  drop_na(acsh) %>%
  #filter(year>1938) %>%
  as.data.frame()
head(monthly.clouds)

# annual.clouds <- cdata %>% group_by(station, year) %>%
#   summarise( acsh = mean(ACSH, na.rm=TRUE),
#              #             month.name = unique(month.name)
#   ) %>%
#   drop_na(acsh) %>%
#   #filter(year>1938) %>%
#   as.data.frame()
# 
# ggplot(monthly.clouds, aes(x=year, y=acsh, color=month.name)) +
#   geom_smooth(formula = y ~ poly(x, 2), se=FALSE)+
#   facet_wrap(~station, ncol = 3) 
# ggsave("Figures/montly clouds trends.pdf", width=8, height=6)
# 
# 
# ggplot(annual.clouds, aes(x=year, y=acsh, color=station)) +
#   geom_smooth(formula = y ~ poly(x, 2), se=FALSE) #+
# # facet_wrap(~station, ncol = 3) 
# ggsave("Figures/annual clouds trends.pdf", width=8, height=6)
# 
# 
# ggplot(monthly.clouds, aes(x=precip, y=acsh, color=station)) +
#   geom_point(alpha=.2)+
#   geom_smooth(formula = y ~ I(x^2)) #+
# ggsave("Figures/precip_clouds.pdf", width=8, height=6)
# 


# Cumulative plots ----

head(cdata)
head(filter(cdata, winter.month==1))
filter(cdata, winter.month==7, winter.year==2025)$DATE

# need to do breaks on real months?
breaks.use = c(1,32,60,91,121,152,182,213,244,274,305,335)#seq(1,365,length.out=13), 
labels.use = winter.month_label

site.data <- site.data %>%
  mutate(highlight_group = case_when(
    winter.year == 2025 ~ "2025",
    winter.year %in% 2021:2024 ~ "2021–2024",
    TRUE ~ NA_character_
  ))


for(i in 1:dim(stations.lat)[1]) {
  site.data <- filter(cdata, station==stations.lat$station[i])
  
  max_data <- group_by(site.data, winter.year) %>%
    summarise(year=median(winter.year),
              max_cum_PRCP = max(cum_PRCP))
  
  mean_precip <- mean(max_data$max_cum_PRCP)
  
  ggplot(site.data, aes(x=winter.day, y=cum_PRCP, group=winter.year)) +
    geom_line(aes(color=winter.year)) +
    scale_color_viridis(option = "D") +  # Other options: "C", "B", "A", "E"
    geom_hline(yintercept=mean_precip, linetype="dashed") +
    geom_line(data=filter(site.data, winter.year==2021), color='forestgreen', size=2)+
    scale_x_continuous(breaks = breaks.use, labels = labels.use) +
    geom_line(data=filter(site.data, winter.year==2022), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2023), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2024), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2025), color="darkorange", size=2)+
    
    geom_text_repel(data=filter(max_data, max_cum_PRCP < 75 | max_cum_PRCP > 125), aes(x=380, y=max_cum_PRCP, label=winter.year),size=2,segment.alpha=0, box.padding=0,max.overlaps=15)+
    expand_limits(x = 400)+
    ggtitle(unique(site.data$station)) + ylab("Cumulative precip (cm)") + xlab("")+
    labs(color='Year') 
  
  ggsave(paste("Figures/",unique(site.data$station),".cumulative precip.pdf", sep=""), width=8, height=6)
  
  
  sub.data <- filter(site.data, winter.year>2000)
  max_data2 <- group_by(sub.data, winter.year) %>%
    summarise(year=median(winter.year),
              max_cum_PRCP = max(cum_PRCP))
  
  mean_precip2 <- mean(max_data2$max_cum_PRCP)
  
  ggplot(sub.data, aes(x=winter.day, y=cum_PRCP, group=winter.year)) +
    geom_line(aes(color=winter.year)) +
    scale_color_viridis(option = "D") +  # Other options: "C", "B", "A", "E"
    scale_x_continuous(breaks = breaks.use, labels = labels.use) +
    geom_text_repel(data=max_data2, aes(x=380, y=max_cum_PRCP, label=winter.year),
                    size=2, segment.alpha=0, box.padding=0, max.overlaps=25) +
    expand_limits(x = 400) +
    ggtitle(unique(site.data$station)) +
    geom_line(data=filter(site.data, winter.year==2021), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2022), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2023), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2024), color='forestgreen', size=2)+
    geom_line(data=filter(site.data, winter.year==2025), color="darkorange", size=2)+
    ylab("Cumulative precip (cm)") + xlab("") +
    labs(color = 'Year')
  
  ggsave(paste("Figures/",unique(site.data$station),".recent cumulative precip.pdf", sep=""), width=8, height=6)
  
} # end station loop


# Eugene only ----
# add labels by the lines
annual.EUG <- filter(annual.climate, station=="Eugene")
## Trends ----
p1<-ggplot(data=annual.EUG, aes(x=winter.year, y=precip, color=station)) +
  geom_point(alpha=.5)+
  geom_hline(yintercept=mean(annual.climate$precip), linetype='dashed')+
  geom_smooth(method=lm)+ # geom_smooth()+
  labs(y="Annual precipitation (cm)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


t1<-ggplot(data=annual.EUG, aes(x=winter.year, y=mean.temp, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$mean.temp), linetype='dashed')+
  labs(y="Mean temperature (C)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal() + 
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
t1
s1<-ggplot(data=annual.EUG, aes(x=winter.year, y=spei1, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$spei1), linetype='dashed') +
  # labs(y="Standardised Precipitation-Evapotranspiration Index", x = "Year", title = "", subtitle = "", caption = "Figure 5. Trends", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
s1

s12<-ggplot(data=annual.EUG, aes(x=winter.year, y=spei12, color=station)) +
  geom_point(alpha=.5)+
  ylab("Drought index (SPEI12)") + xlab("Year")+
  geom_smooth(method=lm) # geom_smooth()+
s12
p1
t1 
s1

p1/t1/s1 + plot_layout(guides = "collect")
ggsave("Figures/Eugene_climate trends.pdf", width=5, height=10)


## cumulative----
eug.data <- filter(cdata, station=="Eugene")
head(eug.data,2)

max_data <- group_by(eug.data, winter.year) %>%
  summarise(year=median(winter.year),
            max_cum_PRCP = max(cum_PRCP))
mean_precip <- mean(max_data$max_cum_PRCP)

# Create a dummy legend key
col2 <- "darkseagreen1"
legend.data <- data.frame(
  x = 0, y = 0,
  group = c("2025", "2021–2024"),
  color = c("darkorange", col2)
)

ggplot(eug.data, aes(x=winter.day, y=cum_PRCP, group=winter.year)) +
  geom_line(aes(color=winter.year))+
  geom_hline(yintercept=mean_precip, linetype="dashed") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_viridis(option = "D") +  # Other options: "C", "B", "A", "E"
  geom_text_repel(data=filter(max_data,winter.year>2020,winter.year<2025), aes(x=380, y=max_cum_PRCP, label=winter.year),
                  size=2, segment.alpha=0, box.padding=0, max.overlaps=25) +
  geom_text_repel(data=filter(max_data,winter.year==2025), aes(x=380, y=max_cum_PRCP, label=winter.year), color="darkorange",
                  size=3, segment.alpha=0, box.padding=0, max.overlaps=25) +
  ggtitle(unique(site.data$station)) +
  # geom_line(data=filter(eug.data, winter.year==2021), color=col2, size=2)+
  # geom_line(data=filter(eug.data, winter.year==2022), color=col2, size=2)+
  # geom_line(data=filter(eug.data, winter.year==2023), color=col2, size=2)+
  # geom_line(data=filter(eug.data, winter.year==2024), color=col2, size=2)+
  geom_line(data=filter(eug.data, winter.year==2025), color="darkorange", size=2)+
  ggtitle(unique(eug.data$station)) + ylab("Cumulative precip (cm)") + xlab("")+
  labs(color='Year') +
  theme_minimal()
ggsave("Figures/Eugene_precip_cumulative.pdf", width=8, height=5)



## monthly----
sub1<-filter(monthly.melt, station=="Eugene")
ppp <- ggplot(filter(sub1,variable=="tmean"), aes(x=winter.month, y=value)) +
  # geom_violin()+
  geom_jitter( size=0.4, alpha=0.9,color='darkred') +
  geom_boxplot(color='darkred', notch=TRUE)+  # 
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("degrees (C)")+xlab("")+
  ylim(0,temp.max.grand)+
  ggtitle(unique(sub1$station))+
  ggtitle("Eugene, Oregon")
ppp

ttt <- ggplot(filter(sub1,variable=="precip"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='blue', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='blue', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("Precip (cm)") + xlab("Month")
# ylim(0,ppt.max.grand)
ttt
spei <- ggplot(filter(sub1,variable=="spei1"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='brown', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='brown', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("SPEI") + xlab("Month")
# ylim(0,ppt.max.grand)
spei


ppp/ttt/spei
ggsave("Figures/eugene mean monthly climate.pdf", width=5, height=10)


# ! IDEA: simulate what would look like with different shelter scenarios.  e.g. week on, week off.... ----
head(eug.data,2)

# %% is called the modulus operator, or sometimes just "modulo."
#  It returns the remainder after division

eug.data_drought <- eug.data %>%
  mutate(PRCP50 = if_else((floor((winter.day - 1) / 7) %% 2) == 0, PRCP, 0),
         PRCP_2wk = if_else((floor((winter.day - 1) / 7) %% 4) < 2, PRCP, 0),
         PRCP_2_1wk = if_else((floor((winter.day - 1) / 7) %% 3) ==0, PRCP, 0))%>%
  mutate(PRCP_Jan_2 = if_else(
    winter.day < 121, 
    PRCP,  # before day 60: normal precipitation
    if_else(
      (floor((winter.day - 121) / 7) %% 4) < 2, 
      0,       # drought weeks
      PRCP     # non-drought weeks
    )
  ),
  PRCP_Jan_1 = if_else(
    winter.day < 121, 
    PRCP,  # before day 60: normal precipitation
    if_else(
      (floor((winter.day - 121) / 7) %% 2) ==0, 
      0,       # drought weeks
      PRCP     # non-drought weeks
    )
  ))

# First, pivot to long form
eug_long <- eug.data_drought %>%
  select(winter.day, winter.year, PRCP, PRCP50) %>%
  pivot_longer(cols = c(PRCP, PRCP50),
               names_to = "type",
               values_to = "precip") %>%
  group_by(winter.year, type) %>%
  arrange(winter.day, .by_group = TRUE) %>%
  mutate(cum_precip = cumsum(precip)) %>%
  ungroup()

ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type)) +
  geom_line() +
  facet_wrap(~ winter.year, scales = "free_x") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP50" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP50" = "50% Precip Reduction")
  ) +
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
ggsave("Figures/eugene half_drought1_separateYears.pdf", width=10, height=10)

main_plot<-ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type, group=paste(winter.year,type))) +
  geom_line() +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP50" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP50" = "50% Precip Reduction")
  ) +
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )+
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       subtitle = '1 week on; 1 week off',
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
main_plot
ggsave("Figures/eugene half_drought1.pdf", width=10, height=8)


annual_totals <- eug_long %>%
  group_by(winter.year, type) %>%
  summarise(total_precip = max(cum_precip), .groups = "drop")

max_precip <- max(annual_totals$total_precip, na.rm = TRUE)

hist_plot <- ggplot(annual_totals, aes(x = total_precip, fill = type)) +
  # geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, boundary = 0, closed = "left") +
  scale_fill_manual(
    values = c("PRCP" = "blue", "PRCP50" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP50" = "50% Precip Reduction")
  ) +
  coord_flip() +
  xlim(0, max_precip) +  # Force same scale as y-axis in main plot
  labs(x = "Annual Total Precipitation (cm)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
hist_plot

# Combine with histogram on the right
combined_plot <- main_plot + hist_plot + plot_layout(widths = c(3, 1))
combined_plot
ggsave("Figures/eugene half_drought.pdf", width=12, height=7)

## Percent----

eug_long_pct <- eug_long %>%
  group_by(winter.year, type) %>%
  mutate(total_annual_precip = max(cum_precip, na.rm = TRUE),
         cum_precip_pct = cum_precip / total_annual_precip * 100) %>%
  ungroup()

ggplot(eug_long_pct, aes(x = winter.day, y = cum_precip_pct, color = type, group = paste(winter.year, type))) +
  geom_line(alpha = 0.6) +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP50" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP50" = "50% Precip Reduction")
  ) +
  labs(
    title = "Cumulative Percent of Annual Precipitation",
    x = "Winter Day",
    y = "Percent of Annual Precipitation (%)"
  ) +
  geom_vline(xintercept=213, linetype='dashed')+
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )
ggsave("Figures/eugene half_drought_percent.pdf", width=7, height=6)


## 2-week droughts ----
# First, pivot to long form
eug_long <- eug.data_drought %>%
  select(winter.day, winter.year, PRCP, PRCP_2wk) %>%
  pivot_longer(cols = c(PRCP, PRCP_2wk),
               names_to = "type",
               values_to = "precip") %>%
  group_by(winter.year, type) %>%
  arrange(winter.day, .by_group = TRUE) %>%
  mutate(cum_precip = cumsum(precip)) %>%
  ungroup()

main_plot<-ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type, group=paste(winter.year,type))) +
  geom_line() +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_2wk" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_2wk" = "50% Precip Reduction")
  ) +
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )+
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       subtitle="2 weeks on, 2 weeks off",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
main_plot

annual_totals <- eug_long %>%
  group_by(winter.year, type) %>%
  summarise(total_precip = max(cum_precip), .groups = "drop")

max_precip <- max(annual_totals$total_precip, na.rm = TRUE)

hist_plot <- ggplot(annual_totals, aes(x = total_precip, fill = type)) +
  # geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, boundary = 0, closed = "left") +
  scale_fill_manual(
    values = c("PRCP" = "blue", "PRCP_2wk" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_2wk" = "50% Precip Reduction")
  ) +
  coord_flip() +
  xlim(0, max_precip) +  # Force same scale as y-axis in main plot
  labs(x = "Annual Total Precipitation (cm)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
hist_plot

# Combine with histogram on the right
combined_plot <- main_plot + hist_plot + plot_layout(widths = c(3, 1))
combined_plot

ggsave("Figures/eugene half_drought2wk.pdf", width=12, height=7)

## 2-week drought, 1 precip ----
# First, pivot to long form
eug_long <- eug.data_drought %>%
  select(winter.day, winter.year, PRCP, PRCP_2_1wk) %>%
  pivot_longer(cols = c(PRCP, PRCP_2_1wk),
               names_to = "type",
               values_to = "precip") %>%
  group_by(winter.year, type) %>%
  arrange(winter.day, .by_group = TRUE) %>%
  mutate(cum_precip = cumsum(precip)) %>%
  ungroup()

main_plot<-ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type, group=paste(winter.year,type))) +
  geom_line() +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_2_1wk" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_2_1wk" = "67% Precip Reduction")
  ) +
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )+
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       subtitle="2 weeks on, 1 weeks off",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
main_plot

annual_totals <- eug_long %>%
  group_by(winter.year, type) %>%
  summarise(total_precip = max(cum_precip), .groups = "drop")

max_precip <- max(annual_totals$total_precip, na.rm = TRUE)

hist_plot <- ggplot(annual_totals, aes(x = total_precip, fill = type)) +
  # geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, boundary = 0, closed = "left") +
  scale_fill_manual(
    values = c("PRCP" = "blue", "PRCP_2_1wk" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_2_1wk" = "67% Precip Reduction")
  ) +
  coord_flip() +
  xlim(0, max_precip) +  # Force same scale as y-axis in main plot
  labs(x = "Annual Total Precipitation (cm)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
hist_plot

# Combine with histogram on the right
combined_plot <- main_plot + hist_plot + plot_layout(widths = c(3, 1))
combined_plot

ggsave("Figures/eugene half_drought2_1wk.pdf", width=12, height=7)


## start Jan 1; then 2 weeks / 2 weeks ----

# First, pivot to long form
eug_long <- eug.data_drought %>%
  select(winter.day, winter.year, PRCP, PRCP_Jan_2) %>%
  pivot_longer(cols = c(PRCP, PRCP_Jan_2),
               names_to = "type",
               values_to = "precip") %>%
  group_by(winter.year, type) %>%
  arrange(winter.day, .by_group = TRUE) %>%
  mutate(cum_precip = cumsum(precip)) %>%
  ungroup()

ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type)) +
  geom_line() +
  facet_wrap(~ winter.year, scales = "free_x") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_Jan_2" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_2" = "50% Precip Reduction")
  ) +
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 

main_plot<-ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type, group=paste(winter.year,type))) +
  geom_line() +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_Jan_2" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_2" = "Starting Jan 1 \n50% Precip Reduction")
  ) +
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )+
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       subtitle="Starting Jan 1: 2 weeks on 2 weeks off",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
main_plot
ggsave("Figures/eugene half_drought2_Jan1.pdf", width=10, height=10)


annual_totals <- eug_long %>%
  group_by(winter.year, type) %>%
  summarise(total_precip = max(cum_precip), .groups = "drop")

max_precip <- max(annual_totals$total_precip, na.rm = TRUE)

hist_plot <- ggplot(annual_totals, aes(x = total_precip, fill = type)) +
  # geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, boundary = 0, closed = "left") +
  scale_fill_manual(
    values = c("PRCP" = "blue", "PRCP_Jan_2" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_2" = "Start Jan 1 \n50% Precip Reduction")
  ) +
  coord_flip() +
  xlim(0, max_precip) +  # Force same scale as y-axis in main plot
  labs(x = "Annual Total Precipitation (cm)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
hist_plot

# Combine with histogram on the right
combined_plot <- main_plot + hist_plot + plot_layout(widths = c(3, 1))
combined_plot
ggsave("Figures/eugene half_drought_Jan1_2wk.pdf", width=12, height=7)


## start Jan 1; then 1 weeks / 1 weeks ----

# First, pivot to long form
eug_long <- eug.data_drought %>%
  select(winter.day, winter.year, PRCP, PRCP_Jan_1) %>%
  pivot_longer(cols = c(PRCP, PRCP_Jan_1),
               names_to = "type",
               values_to = "precip") %>%
  group_by(winter.year, type) %>%
  arrange(winter.day, .by_group = TRUE) %>%
  mutate(cum_precip = cumsum(precip)) %>%
  ungroup()

ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type)) +
  geom_line() +
  facet_wrap(~ winter.year, scales = "free_x") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_Jan_1" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_1" = "50% Precip Reduction")
  ) +
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 

main_plot<-ggplot(eug_long, aes(x = winter.day, y = cum_precip, color = type, group=paste(winter.year,type))) +
  geom_line() +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_manual(
    name = "Scenario",
    values = c("PRCP" = "blue", "PRCP_Jan_1" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_1" = "Starting Jan 1 \n50% Precip Reduction")
  ) +
  theme(
    legend.position = c(0.05, 0.95),         # Position in panel coordinates
    legend.position.inside = "panel",        # << Key line to place it inside the plot
    legend.justification = c(0, 1),          # Align top-left corner of the legend box
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.just = "left"                 # Align content inside the legend box
  )+
  labs(title = "Cumulative Precipitation: Observed vs 50% Reduction",
       subtitle="Starting Jan 1: 1 week on 1 week off",
       x = "Winter Day",
       y = "Cumulative Precipitation (cm)") 
main_plot
ggsave("Figures/eugene half_drought1_Jan1.pdf", width=10, height=10)


annual_totals <- eug_long %>%
  group_by(winter.year, type) %>%
  summarise(total_precip = max(cum_precip), .groups = "drop")

max_precip <- max(annual_totals$total_precip, na.rm = TRUE)

hist_plot <- ggplot(annual_totals, aes(x = total_precip, fill = type)) +
  # geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, boundary = 0, closed = "left") +
  scale_fill_manual(
    values = c("PRCP" = "blue", "PRCP_Jan_1" = "red"),
    labels = c("PRCP" = "Observed Precipitation", 
               "PRCP_Jan_1" = "Start Jan 1 \n50% Precip Reduction")
  ) +
  coord_flip() +
  xlim(0, max_precip) +  # Force same scale as y-axis in main plot
  labs(x = "Annual Total Precipitation (cm)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
hist_plot

# Combine with histogram on the right
combined_plot <- main_plot + hist_plot + plot_layout(widths = c(3, 1))
combined_plot
ggsave("Figures/eugene half_drought_Jan1_1wk.pdf", width=12, height=7)




# Portland only ----

annual.PDX <- filter(annual.climate, station=="Portland")

p1<-ggplot(data=annual.PDX, aes(x=winter.year, y=precip, color=station)) +
  geom_point(alpha=.5)+
  geom_hline(yintercept=mean(annual.climate$precip), linetype='dashed')+
  geom_smooth(method=lm)+ # geom_smooth()+
  labs(y="Annual precipitation (cm)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


t1<-ggplot(data=annual.PDX, aes(x=winter.year, y=mean.temp, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$mean.temp), linetype='dashed')+
  labs(y="Mean temperature (C)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal() + 
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
t1
s1<-ggplot(data=annual.PDX, aes(x=winter.year, y=spei1, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$spei1), linetype='dashed') +
  # labs(y="Standardised Precipitation-Evapotranspiration Index", x = "Year", title = "", subtitle = "", caption = "Figure 5. Trends", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
s1

p1/t1/s1 + plot_layout(guides = "collect")
ggsave("Figures/PDX_climate trends.pdf", width=5, height=10)


## cumulative----
PDX.data <- filter(cdata, station=="Portland")
head(PDX.data,2)

max_data <- group_by(PDX.data, winter.year) %>%
  summarise(year=median(winter.year),
            max_cum_PRCP = max(cum_PRCP))
mean_precip <- mean(max_data$max_cum_PRCP)

ggplot(PDX.data, aes(x=winter.day, y=cum_PRCP, group=winter.year)) +
  geom_line(aes(color=winter.year))+
  geom_hline(yintercept=mean_precip, linetype="dashed") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_viridis(option = "D") +  # Other options: "C", "B", "A", "E"
  geom_text_repel(data=filter(max_data,winter.year>2020,winter.year<2025), aes(x=380, y=max_cum_PRCP, label=winter.year),
                  size=2, segment.alpha=0, box.padding=0, max.overlaps=25) +
  geom_text_repel(data=filter(max_data,winter.year==2025), aes(x=380, y=max_cum_PRCP, label=winter.year), color="darkorange",
                  size=3, segment.alpha=0, box.padding=0, max.overlaps=25) +
  geom_line(data=filter(PDX.data, winter.year==2025), color="darkorange", size=2)+
  ggtitle(unique(PDX.data$station)) + ylab("Cumulative precip (cm)") + xlab("")+
  labs(color='Year') +
  theme_minimal()
ggsave("Figures/PDX_precip_cumulative.pdf", width=8, height=5)




## monthly----
sub1<-filter(monthly.melt, station=="Portland")
ppp <- ggplot(filter(sub1,variable=="tmean"), aes(x=winter.month, y=value)) +
  # geom_violin()+
  geom_jitter( size=0.4, alpha=0.9,color='darkred') +
  geom_boxplot(color='darkred', notch=TRUE)+  # 
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("degrees (C)")+xlab("")+
  ylim(0,temp.max.grand)+
  ggtitle(unique(sub1$station))+
  ggtitle("PDXene, Oregon")
ppp

ttt <- ggplot(filter(sub1,variable=="precip"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='blue', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='blue', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("Precip (cm)") + xlab("Month")
# ylim(0,ppt.max.grand)
ttt
spei <- ggplot(filter(sub1,variable=="spei1"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='brown', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='brown', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("SPEI") + xlab("Month")
# ylim(0,ppt.max.grand)
spei


ppp/ttt/spei
ggsave("Figures/PDX mean monthly climate.pdf", width=5, height=10)



# Salem only ----

annual.Salem <- filter(annual.climate, station=="Salem")

p1<-ggplot(data=annual.Salem, aes(x=winter.year, y=precip, color=station)) +
  geom_point(alpha=.5)+
  geom_hline(yintercept=mean(annual.climate$precip), linetype='dashed')+
  geom_smooth(method=lm)+ # geom_smooth()+
  labs(y="Annual precipitation (cm)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


t1<-ggplot(data=annual.Salem, aes(x=winter.year, y=mean.temp, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$mean.temp), linetype='dashed')+
  labs(y="Mean temperature (C)", x = "Year", title = "", subtitle = "", caption = "", tag = "") +
  theme_minimal() + 
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
t1
s1<-ggplot(data=annual.Salem, aes(x=winter.year, y=spei1, color=station)) +
  geom_point(alpha=.5)+
  geom_smooth(method=lm)+ # geom_smooth()+
  geom_hline(yintercept=mean(annual.climate$spei1), linetype='dashed') +
  # labs(y="Standardised Precipitation-Evapotranspiration Index", x = "Year", title = "", subtitle = "", caption = "Figure 5. Trends", tag = "") +
  theme_minimal()+
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
s1

p1/t1/s1 + plot_layout(guides = "collect")
ggsave("Figures/Salem_climate trends.pdf", width=5, height=10)


## cumulative----
Salem.data <- filter(cdata, station=="Salem")
head(Salem.data,2)

max_data <- group_by(Salem.data, winter.year) %>%
  summarise(year=median(winter.year),
            max_cum_PRCP = max(cum_PRCP))
mean_precip <- mean(max_data$max_cum_PRCP)

ggplot(Salem.data, aes(x=winter.day, y=cum_PRCP, group=winter.year)) +
  geom_line(aes(color=winter.year))+
  geom_hline(yintercept=mean_precip, linetype="dashed") +
  scale_x_continuous(breaks = breaks.use, labels = labels.use) +
  scale_color_viridis(option = "D") +  # Other options: "C", "B", "A", "E"
  geom_text_repel(data=filter(max_data,winter.year>2020,winter.year<2025), aes(x=380, y=max_cum_PRCP, label=winter.year),
                  size=2, segment.alpha=0, box.padding=0, max.overlaps=25) +
  geom_text_repel(data=filter(max_data,winter.year==2025), aes(x=380, y=max_cum_PRCP, label=winter.year), color="darkorange",
                  size=3, segment.alpha=0, box.padding=0, max.overlaps=25) +
  geom_line(data=filter(Salem.data, winter.year==2025), color="darkorange", size=2)+
  ggtitle(unique(Salem.data$station)) + ylab("Cumulative precip (cm)") + xlab("")+
  labs(color='Year') +
  theme_minimal()
ggsave("Figures/Salem_precip_cumulative.pdf", width=8, height=5)




## monthly----
sub1<-filter(monthly.melt, station=="Salem")
ppp <- ggplot(filter(sub1,variable=="tmean"), aes(x=winter.month, y=value)) +
  # geom_violin()+
  geom_jitter( size=0.4, alpha=0.9,color='darkred') +
  geom_boxplot(color='darkred', notch=TRUE)+  # 
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("degrees (C)")+xlab("")+
  ylim(0,temp.max.grand)+
  ggtitle(unique(sub1$station))+
  ggtitle("Salem, Oregon")
ppp

ttt <- ggplot(filter(sub1,variable=="precip"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='blue', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='blue', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("Precip (cm)") + xlab("Month")
# ylim(0,ppt.max.grand)
ttt
spei <- ggplot(filter(sub1,variable=="spei1"), aes(x=winter.month, y=value)) +
  geom_boxplot(color='brown', notch=TRUE)+  # 
  geom_jitter( size=0.4, alpha=0.9,color='brown', width=.1) +
  scale_x_discrete(labels = winter.month_label)+ 
  ylab("SPEI") + xlab("Month")
# ylim(0,ppt.max.grand)
spei


ppp/ttt/spei
ggsave("Figures/Salem mean monthly climate.pdf", width=5, height=10)





