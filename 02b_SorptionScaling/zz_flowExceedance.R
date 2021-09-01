#############################################################################################
# Required packages
#############################################################################################

library(mgcv)
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(readr)
library(readxl)
library(zoo)
library(openxlsx)
library(discharge)

#############################################################################################
# Discharge DAta
#############################################################################################
# USGS hydro network linked data index
# https://waterdata.usgs.gov/blog/nldi-intro/
# Little function for getting Q data I need
StreamQGrabFun <- function(siteNo, pCode, start.date, end.date){
  Qdat<- readNWISuv(siteNumbers = siteNo,
                    parameterCd = pCode,
                    startDate = start.date,
                    endDate = end.date,
                    tz = "America/New_York")
  Qdat <-  renameNWISColumns(Qdat)
  Qdat
}

# SOUTH TURKEY FOOt
# https://waterdata.usgs.gov/monitoring-location/04185440/#parameterCode=00060
qTurk <- StreamQGrabFun(siteNo = "04192599",
                        pCode = "00060",
                        start.date = "2014-12-17",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qTurk, aes(y = Flow_Inst, x = dateTime)) +geom_line()

# FLOW EXCEEDANCE
qTurkfe <- sort(qTurk$Flow_Inst, decreasing = T)
qTurkfe2 <- data.frame(x = 100/length(qTurkfe) * 1:length(qTurkfe), y = qTurkfe)

qTurkfe2[qTurkfe2$x > 24.9 & qTurkfe2$x <25.1,]
#125

TurkTargetFlow <- quantile(qTurk$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; TurkTargetFlow

ggplot(qTurkfe2, aes(y = log(y), x = x)) +
  geom_line() +
  geom_hline(yintercept = log(TurkTargetFlow))

ggplot(qTurkfe2, aes(y = y, x = x)) +
  geom_line() +
  geom_hline(yintercept = TurkTargetFlow)





# UNNAMED TRIB TO LOST CREEK -> MR_def -> MR_wat
# https://waterdata.usgs.gov/monitoring-location/04185440/#parameterCode=00060
qUTLC <- StreamQGrabFun(siteNo = "04185440",
                        pCode = "00060",
                        start.date = "1990-10-01",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# this is missing two Years
ggplot(qUTLC, aes(y = Flow_Inst, x = dateTime)) +geom_line()


# FLOW EXCEEDANCE
qUTLCfe <- sort(qUTLC$Flow_Inst, decreasing = T)
qUTLCfe2 <- data.frame(x = 100/length(qUTLCfe) * 1:length(qUTLCfe), y = qUTLCfe)

qUTLCfe2[qUTLCfe2$x > 24.99 & qUTLCfe2$x <25.01,]
# 2.31

qUTLCTargetFlow <- quantile(qUTLC$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; qUTLCTargetFlow

ggplot(qUTLCfe2, aes(y = log(y), x = x)) +
  geom_line() +
  geom_hline(yintercept = log(qUTLCTargetFlow))


# WEST
# https://waterdata.usgs.gov/monitoring-location/04185440/#parameterCode=00060
qWest <- StreamQGrabFun(siteNo = "04192574",
                        pCode = "00060",
                        start.date = "2014-12-30",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qWest, aes(y = Flow_Inst, x = dateTime)) +geom_point()

# FLOW EXCEEDANCE
qWestfe <- sort(qWest$Flow_Inst, decreasing = T)
qWestfe2 <- data.frame(x = 100/length(qWestfe) * 1:length(qWestfe), y = qWestfe)

qWestfe2[qWestfe2$x > 24.99 & qWestfe2$x <25.01,]
# 19.8

qWestTargetFlow <- quantile(qWest$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; qWestTargetFlow

ggplot(qWestfe2, aes(y = log(y), x = x)) +
  geom_line() +
  geom_hline(yintercept = log(qWestTargetFlow))


# Maumee waterville
# this only goes back to 88
qMaum <- StreamQGrabFun(siteNo = "04193500",
                        pCode = "00060",
                        start.date = "1975-10-01",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qMaum, aes(y = Flow_Inst, x = dateTime)) +geom_point()

# FLOW EXCEEDANCE
qMaumfe <- sort(qMaum$Flow_Inst, decreasing = T)
qMaumfe2 <- data.frame(x = 100/length(qMaumfe) * 1:length(qMaumfe), y = qMaumfe)

qMaumfe2[qMaumfe2$x > 24.99 & qMaumfe2$x <25.01,]
# 7870

qMaumTargetFlow <- quantile(qMaum$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; qMaumTargetFlow

ggplot(qMaumfe2, aes(y = log(y), x = x)) +
  geom_line() +
  geom_hline(yintercept = log(qMaumTargetFlow))
