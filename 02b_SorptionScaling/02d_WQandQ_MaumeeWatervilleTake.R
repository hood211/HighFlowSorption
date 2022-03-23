# JMH Jan 2022
# GEts maumee river discharge and water quality data
# fills in gaps for SRP, SS, and TP load and concentration
# combines all data then calculates daily loads
# Exports daily load data

# libraries ----
library(tidyverse)
library(mgcv)
library(zoo)
library(MuMIn)
library(mgcViz)

# Data ----
## Discharge Data ----

# daily data 
# downloaded on 29 Dec 21 from the USGS website
# daily Q for Maumee at Waterville
# can't access these data with the function used for other streams
# no missing days in this
qMaum <- read.csv(file.path(here::here("01_RawData"), "MaumeeWatervilleDAilyQ1975to2019.csv")) %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%y", tz = "EST"),
         Qm3s_Daily = Q_cfs/35.31467) %>% 
  select(Date, Qm3s_Daily) 

### Percentile flow ----
# these are similar to what is on the gaging station site which has 87 water years = 212.4 m3/s
# Percentile flow - daily data
MaumTargetFlow_Daily <- quantile(qMaum$Qm3s_Daily, probs = 0.75, na.rm = TRUE)[[1]]  # in m3/s


## Water quality data ----
# downloaded from Heidleberg then converted to csv
  wqMaum <- read.csv("01_RawData/00wq_maumeedata.csv")[,1:15] %>% 
    select(-Days741001, -Future, -Month) %>% 
    mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M", tz = "EST"),
           Date = as.Date(dateTime, tz = "EST"),
           M = as.numeric(strftime(Date, format = "%m")),
           Y = as.numeric(strftime(Date, format = "%Y")),
           Qm3s = as.numeric(QCFS)/35.31467,
           SS_mgL = ifelse(SS_mgL == -9, as.numeric("NA"), SS_mgL),
           SRP_mgL = ifelse(SRP_mgL == -9, as.numeric("NA"), SRP_mgL),
           TP_mgL = ifelse(TP_mgL == -9, as.numeric("NA"), TP_mgL),
           # set <= 0 value to 1/2 minimum non-zero value
           SS_mgL = ifelse(SS_mgL <= 0, 0.200/2, SS_mgL), # >0 min = 0.2
           SRP_mgL = ifelse(SRP_mgL <= 0, 0.001/2, SRP_mgL), # >0 min = 0.001
           TP_mgL = ifelse(TP_mgL <= 0, 0.025/2, TP_mgL)) %>%   # >0 min = 0.025
    # remove days between big gap
    filter(Date <= as.POSIXct("1978-09-30",format = "%Y-%m-%d", tz = "EST") |
             Date >= as.POSIXct("1981-10-13",format = "%Y-%m-%d", tz = "EST")) %>% 
    mutate(# set one negative Q to positive value as 1/2 min
           Qm3s = ifelse(Qm3s <= 0, 0.368/2,Qm3s)) %>% 
    select(dateTime, Date, Y, M, SS_mgL:SRP_mgL, Qm3s) 
    
  # fix wrong time
    wqMaum[wqMaum$dateTime == as.POSIXct("0016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S", tz = "EST"),]$dateTime <-  as.POSIXct("2016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S", tz = "EST")
    
    wqMaum[wqMaum$Date == as.Date("0016-06-27", format= "%Y-%m-%d", tz = "EST"),]$Date <-  as.Date("2016-06-27", format= "%Y-%m-%d", tz = "EST")
    wqMaum$Y <- ifelse(wqMaum$Y == "16", 2016, wqMaum$Y)
    
# file goes from 10 Jan 1976 to 4 Oct 2019 with no missing days
# note that I can not calculate annual loads for 2019
  wqMaum2 <- wqMaum %>%
        # join this to fill in missing Q values, but left joining so I don't add missing dates
        left_join(qMaum, by = "Date") %>% 
        # fills in 12 rows  missing Q data
        mutate(Qm3s = ifelse(is.na(Qm3s), Qm3s_Daily, Qm3s)) %>% 
        select(-Qm3s_Daily)


# calculate daily loads ----
## deal with missing data ----
### split into separate df and calc hourly loads ----

# Fill in missing data - GAMS
  # used model 4 to evaluate several family options for each response including 
  # Gamma-log, Gamma-identity, Gaussian-log, Gaussian-identity, scat-log, and scat-identity
  # also tw with sqrt, log, identity, power
  
## SS load ----
### prepare data ----
wqMaum_SS <- wqMaum2 %>% 
  select(dateTime:SS_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60,
         # units g/h - mg/L * 1 = g/m3
    HrLoad_gSSh = SS_mgL * Qm3h) %>% 
  filter(!is.na(SS_mgL)) %>% 
    mutate(SmpTimeWindowDay = "NA") %>% 
    arrange(dateTime) 
  
### gam models ----
# check out the data
ggplot(wqMaum_SS, aes(y = HrLoad_gSSh, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_SS$HrLoad_gSSh); qqline(wqMaum_SS$HrLoad_gSSh)
qqnorm(log(wqMaum_SS$HrLoad_gSSh)); qqline(log(wqMaum_SS$HrLoad_gSSh))
hist(wqMaum_SS$HrLoad_gSSh)


# Build models
gmSS1 <- gam(HrLoad_gSSh ~ te(Qm3h,Y) + s(M,bs = "cc"),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS2 <- gam(HrLoad_gSSh ~ te(Qm3h,Y) + M,
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS3 <- gam(HrLoad_gSSh ~ te(Qm3h,Y),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS4 <- gam(HrLoad_gSSh ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS5 <- gam(HrLoad_gSSh ~ s(Qm3h) + s(Y),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS6 <- gam(HrLoad_gSSh ~ s(Qm3h) + s(M,bs = "cc"),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS7 <- gam(HrLoad_gSSh ~ s(Qm3h),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS8 <- gam(HrLoad_gSSh ~ s(Qm3h) + s(Y) + M,
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

gmSS9 <- gam(HrLoad_gSSh ~ s(Qm3h) + M,
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "ML")

### model selection/eval
model.sel(gmSS1,gmSS2,gmSS3,gmSS4,gmSS5,gmSS6,gmSS7,gmSS8,gmSS9)
check(getViz(gmSS1))
plot(gmSS1)
summary(gmSS1)

gmSS_best <- gam(HrLoad_gSSh ~ te(Qm3h,Y) + s(M,bs = "cc"),
             data = wqMaum_SS,
             family = Gamma(link = "log"),
             method = "REML")


## SS concentration ----
### prepare data ----
wqMaum_SSc <- wqMaum2 %>% 
  select(dateTime:SS_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60) %>% 
  filter(!is.na(SS_mgL)) %>% 
  arrange(dateTime) 

### gam models ----
# check out the data
ggplot(wqMaum_SSc, aes(y = SS_mgL, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_SSc$SS_mgL); qqline(wqMaum_SSc$SS_mgL)
qqnorm(log(wqMaum_SSc$SS_mgL)); qqline(log(wqMaum_SSc$SS_mgL))
hist(wqMaum_SSc$SS_mgL)


# model selection
gmSSc1 <- gam(SS_mgL ~ te(Qm3h,Y) + s(M,bs = "cc"),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc2 <- gam(SS_mgL ~ te(Qm3h,Y) + M,
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc3 <- gam(SS_mgL ~ te(Qm3h,Y),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc4 <- gam(SS_mgL ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc5 <- gam(SS_mgL ~ s(Qm3h) + s(Y),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc6 <- gam(SS_mgL ~ s(Qm3h) + s(M,bs = "cc"),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc7 <- gam(SS_mgL ~ s(Qm3h),
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc8 <- gam(SS_mgL ~ s(Qm3h) + s(Y) + M,
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

gmSSc9 <- gam(SS_mgL ~ s(Qm3h) + M,
              data = wqMaum_SSc,
              family = tw(link = "log"),
              method = "ML")

model.sel(gmSSc1,gmSSc2,gmSSc3,gmSSc4,gmSSc5,gmSSc6,gmSSc7,gmSSc8,gmSSc9)
check(getViz(gmSSc4))
plot(gmSSc4)
summary(gmSSc4)

gmSSc_best <- gam(SS_mgL ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
                  data = wqMaum_SSc,
                  family = tw(link = "log"),
                  method = "REML")


## SRP load ----
### prepare data ----
wqMaum_SRP <- wqMaum2 %>% 
  select(dateTime:M, SRP_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60,
         # units g/h - mg/L * 1 = g/m3
         HrLoad_gSRPh = SRP_mgL * Qm3h) %>% 
  filter(!is.na(SRP_mgL)) %>% 
  arrange(dateTime) 


### gam models ----
# check out the data
ggplot(wqMaum_SRP, aes(y = HrLoad_gSRPh, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_SRP$HrLoad_gSRPh); qqline(wqMaum_SRP$HrLoad_gSRPh)
qqnorm(log(wqMaum_SRP$HrLoad_gSRPh)); qqline(log(wqMaum_SRP$HrLoad_gSRPh))
hist(wqMaum_SRP$HrLoad_gSRPh)


# model selection
gmSRP1 <- gam(HrLoad_gSRPh ~ te(Qm3h,Y) + s(M,bs = "cc"),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP2 <- gam(HrLoad_gSRPh ~ te(Qm3h,Y) + M,
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP3 <- gam(HrLoad_gSRPh ~ te(Qm3h,Y),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP4 <- gam(HrLoad_gSRPh ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP5 <- gam(HrLoad_gSRPh ~ s(Qm3h) + s(Y),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP6 <- gam(HrLoad_gSRPh ~ s(Qm3h) + s(M,bs = "cc"),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP7 <- gam(HrLoad_gSRPh ~ s(Qm3h),
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP8 <- gam(HrLoad_gSRPh ~ s(Qm3h) + s(Y) + M,
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

gmSRP9 <- gam(HrLoad_gSRPh ~ s(Qm3h) + M,
             data = wqMaum_SRP,
             family = tw(link = "sqrt"),
             method = "ML")

# model selection/eval
model.sel(gmSRP1,gmSRP2,gmSRP3,gmSRP4,gmSRP5,gmSRP6,gmSRP7,gmSRP8,gmSRP9)
check(getViz(gmSRP1))
plot(gmSRP1)
summary(gmSRP1)

gmSRP_best <- gam(HrLoad_gSRPh ~ te(Qm3h,Y) + s(M,bs = "cc"),
              data = wqMaum_SRP,
              family = tw(link = "sqrt"),
              method = "REML")

## SRP concentration ----
### prepare data ----
wqMaum_SRPc <- wqMaum2 %>% 
  select(dateTime:M, SRP_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60) %>% 
  filter(!is.na(SRP_mgL)) %>% 
  arrange(dateTime) 


### gam models ----
# check out the data
ggplot(wqMaum_SRPc, aes(y = SRP_mgL, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_SRPc$SRP_mgL); qqline(wqMaum_SRPc$SRP_mgL)
qqnorm(log(wqMaum_SRPc$SRP_mgL)); qqline(log(wqMaum_SRPc$SRP_mgL))
hist(wqMaum_SRPc$SRP_mgL)


# model selection
gmSRPc1 <- gam(SRP_mgL ~ te(Qm3h,Y) + s(M,bs = "cc"),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc2 <- gam(SRP_mgL ~ te(Qm3h,Y) + M,
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc3 <- gam(SRP_mgL ~ te(Qm3h,Y),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc4 <- gam(SRP_mgL ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc5 <- gam(SRP_mgL ~ s(Qm3h) + s(Y),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc6 <- gam(SRP_mgL ~ s(Qm3h) + s(M,bs = "cc"),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc7 <- gam(SRP_mgL ~ s(Qm3h),
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc8 <- gam(SRP_mgL ~ s(Qm3h) + s(Y) + M,
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

gmSRPc9 <- gam(SRP_mgL ~ s(Qm3h) + M,
               data = wqMaum_SRPc,
               family = tw(link = "sqrt"),
               method = "ML")

model.sel(gmSRPc1,gmSRPc2,gmSRPc3,gmSRPc4,gmSRPc5,gmSRPc6,gmSRPc7,gmSRPc8,gmSRPc9)
check(getViz(gmSRPc1))
plot(gmSRPc1)
summary(gmSRPc1)

# model selection eval
gmSRPc_best <- gam(SRP_mgL ~ te(Qm3h,Y) + s(M,bs = "cc"),
                   data = wqMaum_SRPc,
                   family = tw(link = "sqrt"),
                   method = "REML")


## TP load ----
### prepare data ----
wqMaum_TP <- wqMaum2 %>% 
  select(dateTime:M, TP_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60,
         # units g/h - mg/L * 1 = g/m3
         HrLoad_gTPh = TP_mgL * Qm3h) %>% 
  filter(!is.na(TP_mgL)) %>% 
  arrange(dateTime) 


### gam models ----
# check out the data
ggplot(wqMaum_TP, aes(y = HrLoad_gTPh, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_TP$HrLoad_gTPh); qqline(wqMaum_TP$HrLoad_gTPh)
qqnorm(log(wqMaum_TP$HrLoad_gTPh)); qqline(log(wqMaum_TP$HrLoad_gTPh))
hist(wqMaum_TP$HrLoad_gTPh)


# model selection
gmTP1 <- gam(HrLoad_gTPh ~ te(Qm3h,Y) + s(M,bs = "cc"),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP2 <- gam(HrLoad_gTPh ~ te(Qm3h,Y) + M,
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP3 <- gam(HrLoad_gTPh ~ te(Qm3h,Y),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP4 <- gam(HrLoad_gTPh ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP5 <- gam(HrLoad_gTPh ~ s(Qm3h) + s(Y),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP6 <- gam(HrLoad_gTPh ~ s(Qm3h) + s(M,bs = "cc"),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP7 <- gam(HrLoad_gTPh ~ s(Qm3h),
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP8 <- gam(HrLoad_gTPh ~ s(Qm3h) + s(Y) + M,
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")

gmTP9 <- gam(HrLoad_gTPh ~ s(Qm3h) + M,
              data = wqMaum_TP,
             family = tw(link = "sqrt"),
              method = "ML")


model.sel(gmTP1,gmTP2,gmTP3,gmTP4,gmTP5,gmTP6,gmTP7,gmTP8,gmTP9)
check(getViz(gmTP1))
plot(gmTP1)
summary(gmTP1)

gmTP_best <- gam(HrLoad_gTPh ~ te(Qm3h,Y) + s(M,bs = "cc"),
                  data = wqMaum_TP,
                 family = tw(link = "sqrt"),
                  method = "REML")


## TP concentrations ----
### prepare data ----
wqMaum_TPc <- wqMaum2 %>% 
  select(dateTime:M, TP_mgL, Qm3s) %>% 
  mutate(Qm3h = Qm3s*60*60) %>% 
  filter(!is.na(TP_mgL)) %>% 
  arrange(dateTime) 


### gam models ----
# check out the data
ggplot(wqMaum_TPc, aes(y = TP_mgL, x = Qm3h)) +
  geom_point()

qqnorm(wqMaum_TPc$TP_mgL); qqline(wqMaum_TPc$TP_mgL)
qqnorm(log(wqMaum_TPc$TP_mgL)); qqline(log(wqMaum_TPc$TP_mgL))
hist(wqMaum_TPc$TP_mgL)


# model selection
gmTPc1 <- gam(TP_mgL ~ te(Qm3h,Y) + s(M,bs = "cc"),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc2 <- gam(TP_mgL ~ te(Qm3h,Y) + M,
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc3 <- gam(TP_mgL ~ te(Qm3h,Y),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

# used this to evaluate several family options including 
# Gamma-log, Gamma-identity, Gaussian-log, Gaussian-identity, scat-log, and scat-identity
gmTPc4 <- gam(TP_mgL ~ s(Qm3h) + s(Y) + s(M,bs = "cc"),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc5 <- gam(TP_mgL ~ s(Qm3h) + s(Y),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc6 <- gam(TP_mgL ~ s(Qm3h) + s(M,bs = "cc"),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc7 <- gam(TP_mgL ~ s(Qm3h),
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc8 <- gam(TP_mgL ~ s(Qm3h) + s(Y) + M,
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")

gmTPc9 <- gam(TP_mgL ~ s(Qm3h) + M,
              data = wqMaum_TPc,
              family = tw(link = "log"),
              method = "ML")


model.sel(gmTPc1,gmTPc2,gmTPc3,gmTPc4,gmTPc5,gmTPc6,gmTPc7,gmTPc8,gmTPc9)
check(getViz(gmTPc1))
plot(gmTPc1)
summary(gmTPc1)

gmTPc_best <- gam(TP_mgL ~ te(Qm3h,Y) + s(M,bs = "cc"),
                  data = wqMaum_TPc,
                  family = tw(link = "log"),
                  method = "REML")



# Recombine files ----
wqMaum3 <- wqMaum_SRP %>% 
          select(-SRP_mgL) %>% 
          full_join(wqMaum_SS %>% 
                      select(dateTime, HrLoad_gSSh), by = "dateTime") %>% 
          full_join(wqMaum_TP %>% 
                      select(dateTime, HrLoad_gTPh), by = "dateTime") %>% 
          full_join(wqMaum_SSc %>% 
                      select(dateTime, SS_mgL), by = "dateTime") %>% 
          full_join(wqMaum_SRPc %>% 
                      select(dateTime, SRP_mgL), by = "dateTime") %>% 
          full_join(wqMaum_TPc %>% 
                      select(dateTime, TP_mgL), by = "dateTime") %>% 
          mutate(Date = as.Date(dateTime, format = "%Y-%m-%d")) %>% 
          # fill in missing dates and Q's
          full_join(qMaum, by = "Date", tz = "EST") %>% 
          # fill in missing M and Y
          mutate(M = as.numeric(strftime(Date, format = "%m")),
                 Y = as.numeric(strftime(Date, format = "%Y"))) %>% 
          select(-Qm3s) %>% 
          # just save this
          mutate(Qm3h_Backup = Qm3h,
                 # fill in NA's with missing values, note that this rewrites Qm3h, for predict
                 Qm3h = ifelse(is.na(Qm3h_Backup), Qm3s_Daily*60*60, Qm3h_Backup),
                 # fill in missing dateTime's with Date and noon
                 dateTime = ifelse(is.na(dateTime), 
                                   as.POSIXct(paste0(Date, "12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST")),
                                   dateTime),
                 dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00", tz = "EST")) %>% 
          mutate(WqSmpGapDay_SS = "NA",
                 WqSmpGapDay_SRP = "NA",
                 WqSmpGapDay_TP = "NA",
                 WqSmpGapDay_SSc = "NA",
                 WqSmpGapDay_SRPc = "NA",
                 WqSmpGapDay_TPc = "NA",
                 # this is daily at the coarsest
                 WqSmpGapDay_Q = "NA") %>% 
          arrange(dateTime)

# Estimate NA's
  # If gap betwen samples is <= 5 than na.approx ELSE gam
  ## Model predictions ----
    wqMaum3$HrLoad_gSSh_gamPred = exp(predict(gmSS_best, wqMaum3)) # log link
    wqMaum3$HrLoad_gSRPh_gamPred = (predict(gmSRP_best, wqMaum3))^2 #sqrt link
    wqMaum3$HrLoad_gTPh_gamPred = (predict(gmTP_best, wqMaum3))^2 #sqrt link
    wqMaum3$SS_mgL_gamPred = exp(predict(gmSSc_best, wqMaum3)) # log link
    wqMaum3$SRP_mgL_gamPred = (predict(gmSRPc_best, wqMaum3))^2 #sqrt link
    wqMaum3$TP_mgL_gamPred = exp(predict(gmTPc_best, wqMaum3)) # log link
    
  
  
  ## NA approx ----
    wqMaum3$HrLoad_gSSh_aproxPred = na.approx(wqMaum3$HrLoad_gSSh)
    wqMaum3$HrLoad_gSRPh_aproxPred = na.approx(wqMaum3$HrLoad_gSRPh)
    wqMaum3$HrLoad_gTPh_aproxPred = na.approx(wqMaum3$HrLoad_gTPh)
    wqMaum3$SS_mgL_aproxPred = na.approx(wqMaum3$SS_mgL)
    wqMaum3$SRP_mgL_aproxPred = na.approx(wqMaum3$SRP_mgL)
    wqMaum3$TP_mgL_aproxPred = na.approx(wqMaum3$TP_mgL)

  ## Calc window b/w water quality samples ----
  ## function that calcs gaps
    SmpGapCalcFun <- function(WQdata, SoluteLoadRaw, SmpGapDay_name){
      # WQdata <- wqMaum3
      # SoluteLoadRaw <- "HrLoad_gSSh"
      # SmpGapDay_name <- "SmpGapDay_SS"
      for(i in 2:dim(WQdata)[1]){
        # i = 22
        dateTime_i = WQdata[i,]$dateTime

        dateTime_im1 = as.POSIXct(ifelse(!is.na(WQdata[i-1,SoluteLoadRaw]),WQdata[i-1,]$dateTime,
                                         WQdata[max(which(which(!is.na(WQdata[,SoluteLoadRaw])) < i)),]$dateTime ),
                                  origin = "1970-01-01 00:00.00", tz = "EST")
        
        SmpGapDay_i = (as.numeric(dateTime_i) - as.numeric(dateTime_im1))/60/60/24
        WQdata[i,SmpGapDay_name] = SmpGapDay_i
      }
      
      WQdata[,SmpGapDay_name] <- as.numeric(WQdata[,SmpGapDay_name])
      WQdata[1,SmpGapDay_name] <- 1
      WQdata
    }
  
  ## calc gaps ----
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "HrLoad_gSSh", SmpGapDay_name = "WqSmpGapDay_SS")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "HrLoad_gSRPh", SmpGapDay_name = "WqSmpGapDay_SRP")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "HrLoad_gTPh", SmpGapDay_name = "WqSmpGapDay_TP")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "SS_mgL", SmpGapDay_name = "WqSmpGapDay_SSc")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "SRP_mgL", SmpGapDay_name = "WqSmpGapDay_SRPc")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "TP_mgL", SmpGapDay_name = "WqSmpGapDay_TPc")
  wqMaum3 <-  SmpGapCalcFun(WQdata = wqMaum3, SoluteLoadRaw = "Qm3h", SmpGapDay_name = "WqSmpGapDay_Q")
  
  # make first Q after big gap 1 day
  wqMaum3[wqMaum3$dateTime == as.POSIXct("1981-10-13 19:00:00",format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$WqSmpGapDay_Q <- 1
  
  ## Clean up data ----
  wqMaum4 <-  wqMaum3 %>% 
              # remove dates between gap
              filter(Date <= as.POSIXct("1978-09-30",format = "%Y-%m-%d", tz = "EST") |
                       Date >= as.POSIXct("1981-10-13 19:00:00",format = "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
              # select best est for missing data
              # gap <= 5 days na.appox, gap >5 days
              mutate(HrLoad_gSSh_bestPred = ifelse(WqSmpGapDay_SS <= 5, HrLoad_gSSh_aproxPred, HrLoad_gSSh_gamPred),
                     HrLoad_gSSh2 = ifelse(is.na(HrLoad_gSSh), HrLoad_gSSh_bestPred, HrLoad_gSSh),
                     HrLoad_gSRPh_bestPred = ifelse(WqSmpGapDay_SRP <= 5, HrLoad_gSRPh_aproxPred, HrLoad_gSRPh_gamPred),
                     HrLoad_gSRPh2 = ifelse(is.na(HrLoad_gSRPh), HrLoad_gSRPh_bestPred, HrLoad_gSRPh),
                     HrLoad_gTPh_bestPred = ifelse(WqSmpGapDay_TP <= 5, HrLoad_gTPh_aproxPred, HrLoad_gTPh_gamPred),
                     HrLoad_gTPh2 = ifelse(is.na(HrLoad_gTPh), HrLoad_gTPh_bestPred, HrLoad_gTPh),
                     SS_mgL_bestPred = ifelse(WqSmpGapDay_SSc <= 5, SS_mgL_aproxPred, SS_mgL_gamPred),
                     SS_mgL2 = ifelse(is.na(SS_mgL), SS_mgL_bestPred, SS_mgL),
                     SRP_mgL_bestPred = ifelse(WqSmpGapDay_SRPc <= 5, SRP_mgL_aproxPred, SRP_mgL_gamPred),
                     SRP_mgL2 = ifelse(is.na(SRP_mgL), SRP_mgL_bestPred, SRP_mgL),
                     TP_mgL_bestPred = ifelse(WqSmpGapDay_TPc <= 5, TP_mgL_aproxPred, TP_mgL_gamPred),
                     TP_mgL2 = ifelse(is.na(TP_mgL), TP_mgL_bestPred, TP_mgL)) 
  

  ## look at predictions----
    # missing values I'm predicting
    # SS load
    ggplot(wqMaum4, aes(y = HrLoad_gSSh_bestPred, x = HrLoad_gSSh)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()
  
    # SS conc
    ggplot(wqMaum4, aes(y = SS_mgL_bestPred, x = SS_mgL)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()
    
    # SRP load
    ggplot(wqMaum4, aes(y = HrLoad_gSRPh_bestPred, x = HrLoad_gSRPh)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()
    
    # SRP conc
    # these remain hard to predict
    # good thing is that predictions for high flow aren't bad, which is what we care about
    ggplot(wqMaum4 %>% 
             mutate(LfHf = ifelse(Qm3h/60/60 >= MaumTargetFlow_Daily, "HF", "LF")), aes(y = SRP_mgL_bestPred, x = SRP_mgL, color = LfHf)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()

    # TP load
    ggplot(wqMaum4, aes(y = HrLoad_gTPh_bestPred, x = HrLoad_gTPh)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()
    
    # TP conc
    ggplot(wqMaum4, aes(y = TP_mgL_bestPred, x = TP_mgL)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      scale_x_log10() +
      scale_y_log10()

# Calculate daily loads ----
    ## Blank hourly data frame ----
    HourlyTS <- as.data.frame(seq(as.POSIXct("1975-01-10 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),
                                  as.POSIXct("2019-10-03 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"), by = "hours"))
    names(HourlyTS) <- "dateTime"
    
    ## combine with Q and loads ----
    wqMaum_hours <- HourlyTS %>% 
                    left_join(wqMaum4 %>% 
                                # only with gaps filled
                                select(dateTime, Qm3h, HrLoad_gSSh2, HrLoad_gSRPh2, HrLoad_gTPh2, SS_mgL2, SRP_mgL2, TP_mgL2),
                              by = "dateTime")
    
    ## sum to daily loads ----
    # fill in last time point
    wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$Qm3h <- 
      wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$Qm3h
    wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gSSh2 <- 
      wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gSSh2
    wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gSRPh2 <- 
      wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gSRPh2
    wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gTPh2 <- 
      wqMaum_hours[wqMaum_hours$dateTime == as.POSIXct("2019-10-03 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "EST"),]$HrLoad_gTPh2

    
    
    ## na.approx & sum ----
    ### Loads ----
    wqMaum_DailyLoads <- wqMaum_hours %>% 
                    # if NA then na.approx
                    mutate(Qm3h = ifelse(is.na(Qm3h), na.approx(Qm3h), Qm3h),
                           HrLoad_gSSh2 = ifelse(is.na(HrLoad_gSSh2), na.approx(HrLoad_gSSh2), HrLoad_gSSh2),
                           HrLoad_gSRPh2 = ifelse(is.na(HrLoad_gSRPh2), na.approx(HrLoad_gSRPh2), HrLoad_gSRPh2),
                           HrLoad_gTPh2 = ifelse(is.na(HrLoad_gTPh2), na.approx(HrLoad_gTPh2), HrLoad_gTPh2)) %>% 
                    mutate(Date = as.Date(dateTime)) %>% 
                    # sum to daily loads
                    group_by(Date) %>% 
                    summarise(across(Qm3h:HrLoad_gTPh2, sum, na.rm = T)) %>% 
                    rename(Qm3day = Qm3h, HrLoad_gSSday = HrLoad_gSSh2, HrLoad_gSRPday = HrLoad_gSRPh2, HrLoad_gTPday = HrLoad_gTPh2)
                    
    ### Conc ----
    wqMaum_DailyConc <- wqMaum_hours %>% 
                    # if NA then na.approx
                    mutate(Qm3h = ifelse(is.na(Qm3h), na.approx(Qm3h), Qm3h),
                           SS_mgL2 = ifelse(is.na(SS_mgL2), na.approx(SS_mgL2), SS_mgL2), 
                           SRP_mgL2 = ifelse(is.na(SRP_mgL2), na.approx(SRP_mgL2), SRP_mgL2), 
                           TP_mgL2 = ifelse(is.na(TP_mgL2), na.approx(TP_mgL2), TP_mgL2)) %>% 
                    mutate(Date = as.Date(dateTime)) %>% 
                    # sum to daily loads
                    group_by(Date) %>% 
                    summarise(across(SS_mgL2:TP_mgL2, mean, na.rm = T)) %>% 
                    rename(SS_mgL = SS_mgL2, SRP_mgL = SRP_mgL2, TP_mgL = TP_mgL2) 
                  
    ### Join loads & conc ----                
    wqMaum_Daily <-  wqMaum_DailyLoads %>% 
                    # join conc and loads
                    full_join(wqMaum_DailyConc, by = "Date") %>% 
                    # join back USGS daily Q, in m3/s
                    # USGS daily Q and summed Q (from above): R2 = 0.98
                    full_join(qMaum, by = "Date") %>% 
                    # remove dates between gap
                    filter(Date <= as.POSIXct("1978-09-30",format = "%Y-%m-%d", tz = "EST") |
                           Date >= as.POSIXct("1981-10-13 19:00:00",format = "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
                    # high flow
                    mutate(HighFlow = ifelse(Qm3s_Daily >= MaumTargetFlow_Daily, "HighFlow", "LowFlow"),
                           Y = as.numeric(strftime(Date, format = "%Y")),
                           M = as.numeric(strftime(Date, format = "%m")),
                           DOY = as.numeric(strftime(Date, format = "%j")),
                           # going to use USGS daily Q instead of summed
                           Qm3day = Qm3s_Daily*60*60*24) %>% 
                    select(Date, Qm3day, gSSday = HrLoad_gSSday, gSRPday = HrLoad_gSRPday, gTPday = HrLoad_gTPday,
                           SS_gm3 = SS_mgL, SRP_gm3 = SRP_mgL, TP_gm3 = TP_mgL, HighFlow:DOY)
    

# export data ----
  # write.csv(wqMaum_Daily, file.path(here::here("04a_generatedDataOnGit"), "02d_MaumeeWatervilleWaterQual.csv"))

# save.image ----
  # save.image(file.path(here::here("03_Rdata"), "02d_WQandQ_MaumeeWatervilleTake3_Rdat"))
  # load(file.path(here::here("03_Rdata"), "02d_WQandQ_MaumeeWatervilleTake3_Rdat"))
