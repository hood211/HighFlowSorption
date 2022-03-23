# JMH; Fall 2020
# Downloads Q data for Unnamed Trib to Lost Creek (UTLC) as well as downstream monitoring stations
# Maumee @ Defiance and Maumee @ Waterville
# Combines and processes Q and solute chem data
# estimates 75% flow
# Builds models to predict missing SRP and SS data, then predicts missing values
# Generates final stream Q, SRP, SS csv for STF

# Libraries ----
library(mgcv)
library(tidyverse)
library(lubridate)
library(dataRetrieval)
# library(discharge)


# Discharge DAta ----
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


# UNNAMED TRIB TO LOST CREEK -> MR_def -> MR_wat
# https://waterdata.usgs.gov/monitoring-location/04185440/#parameterCode=00060
qUTLC <- StreamQGrabFun(siteNo = "04185440",
                        pCode = "00060",
                        start.date = "1990-10-01",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# this is missing two Years
ggplot(qUTLC, aes(y = Flow_Inst, x = dateTime)) +geom_point()

# MAUMEE NEAR DEFIANCE
qMRdef <- StreamQGrabFun(siteNo = "04192500",
                         pCode = "00060",
                         start.date = "1990-10-01",
                         end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# MAUMEE NEAR WATERVILLE
qMRwat <- StreamQGrabFun(siteNo = "04193500",
                         pCode = "00060",
                         start.date = "1987-10-01",
                         end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))


# Water quality data ----
wqUTLC <- read.csv("01_RawData/00wq_lostcreekdata.csv") %>% 
  mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M"))%>% 
    select(-Days741001, -Future, -Month)



# just take a look
ggplot(wqUTLC, aes(y = QCFS, x = dateTime)) +geom_line()


# Combine water quality and Q data ----
#round data to the nearest 15 min interval
# TZ message doesn't seem to impact the results
wqUTLC.15min <-  wqUTLC %>% 
  mutate(dateTime.15mins = round_date(dateTime, unit = "15 minutes")) %>% 
  arrange(dateTime.15mins)

UTLC <- qUTLC %>% 
  full_join(wqUTLC.15min, by = c("dateTime" = "dateTime.15mins")) %>% 
  select(-agency_cd, -dateTime.y, -TP_mgL, -c(NO23_mgL:Cond_umho)) %>% 
  filter(dateTime >= as.POSIXct("2008-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>% #big gap in Q data, starting after that
  filter(dateTime < as.POSIXct("2019-07-10 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(dateTime) %>% 
  mutate(SS_mgL = ifelse(SS_mgL <= 0, 0, SS_mgL),
         SRP_mgL = ifelse(SRP_mgL <= 0, 0, SRP_mgL)) %>% 
  mutate(Y = as.numeric(as.character(strftime(dateTime, format = "%Y"))),
         M = as.numeric(as.character(strftime(dateTime, format = "%m"))),
         DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))))

# Check alignment of samples
ggplot() +
  geom_point(data = UTLC, aes(y = Flow_Inst, x = DOY)) +
  geom_point(data = UTLC, aes(y = QCFS, x = DOY), color = "red") +
  facet_grid(Y ~., scales = "free_y")


ggplot(UTLC, aes(y = Flow_Inst, x = DOY))+
  geom_point(size = 0.5) +
  facet_grid(Y~.)

# Determine the 75% flow ----
# missing some early season data in 2009-2011
UTLCTargetFlow <- quantile(UTLC$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; UTLCTargetFlow

ggplot(UTLC %>% 
         mutate(AboveTargetFlow = ifelse(Flow_Inst > UTLCTargetFlow, "above", "below")), aes(y = log(Flow_Inst), x = DOY, color = AboveTargetFlow)) +
  geom_point(size = 0.5) +
  facet_grid(Y~.)

# subset to flows above 75%
UTLCgTar <- UTLC %>% 
  filter(Flow_Inst > UTLCTargetFlow) %>% 
  mutate(lSRP_mgL = log(SRP_mgL + 0.0001),
         lSS_mgL = log(SS_mgL + 0.0001),
         lFlow_Inst = log(Flow_Inst)) %>% 
  filter(SS_mgL > 0.01) %>% 
  filter(SRP_mgL > 0.01)
  

# predict missing SRP data ----
ggplot(UTLCgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(Y))) +
  geom_point()

ggplot(UTLCgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(M))) +
  geom_point()

## models ----
# models to predict missing data
# this is the best model, high R2, yet simple
gm_SRP_UTLC <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc"), data = UTLCgTar); summary(gm_SRP_UTLC$gam)
gm_SRP_UTLC1 <- gamm(lSRP_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc"), data = UTLCgTar); summary(gm_SRP_UTLC1$gam)
gm_SRP_UTLC2 <- gamm(lSRP_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = UTLCgTar); summary(gm_SRP_UTLC2$gam)
gm_SRP_UTLC3 <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc") + s(Y, bs = "re"), data = UTLCgTar); summary(gm_SRP_UTLC3$gam)

# use model to predict values
UTLCgTar$lSRP_mgLp <- predict.gam(gm_SRP_UTLC$gam, UTLCgTar, type = "response")

# predicted v. response
ggplot(UTLCgTar, aes(y = lSRP_mgLp, x = lSRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(UTLCgTar, aes(y = exp(lSRP_mgLp), x = SRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# predict missing SS data ----
ggplot(UTLCgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(Y))) +
  geom_point()

ggplot(UTLCgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(M))) +
  geom_point()

# models ----
# models to predict missing data
gm_SSUTLC0 <- gamm(lSS_mgL ~ s(lFlow_Inst), data = UTLCgTar); summary(gm_SSUTLC0$gam)
# best by R2 - also like the Flow spline
gm_SSUTLC <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc"), data = UTLCgTar); summary(gm_SSUTLC$gam)
gm_SSUTLCa <- gamm(lSS_mgL ~ s(Flow_Inst) + s(M, bs = "cc"), data = UTLCgTar); summary(gm_SSUTLCa$gam)
gm_SSUTLC1 <- gamm(lSS_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc"), data = UTLCgTar); summary(gm_SSUTLC1$gam)
gm_SSUTLC2 <- gamm(lSS_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = UTLCgTar); summary(gm_SSUTLC2$gam)
gm_SSUTLC3 <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc") + s(Y, bs = "re"), data = UTLCgTar); summary(gm_SSUTLC3$gam)


# use model to predict values
UTLCgTar$lSS_mgLp <- predict.gam(gm_SSUTLC$gam, UTLCgTar, type = "response")

# predicted v. response
ggplot(UTLCgTar, aes(y = lSS_mgLp, x = lSS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(UTLCgTar, aes(y = exp(lSS_mgLp), x = SS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


# back transform
UTLCgTar2 <- UTLCgTar %>% 
  mutate(SRP_mgLp = exp(lSRP_mgL),
         SS_mgLp = exp(lSS_mgLp))

# Recombine ----
# predict SRP and SS at flows greater than target
# 18065 zeros in flow!
UTLC2 <- UTLC %>% 
  mutate(lFlow_Inst = ifelse(Flow_Inst > 0, log(Flow_Inst), as.numeric("NA")))

# this makes calcs for low flows, those values have not been validated
# values have not been validated and are not used, but numbers retained to maintain code continuity
UTLC2$lSRP_mgLp = predict.gam(gm_SRP_UTLC$gam, UTLC2, type = "response")
UTLC2$lSS_mgLp = predict.gam(gm_SSUTLC$gam, UTLC2, type = "response")


UTLCfin <- UTLC2 %>% 
  select(dateTime, M, Flow_Inst, QCFS, SS_mgL, SRP_mgL, lSRP_mgLp, lSS_mgLp) %>% 
  mutate(SRP_mgLp = exp(lSRP_mgLp),
         SS_mgLp = exp(lSS_mgLp)) %>% 
  select(dateTime, QCFS_utlc = Flow_Inst, SS_mgL, SRP_mgL, SS_mgLp, SRP_mgLp) %>% 
  mutate(gTarPerFlow = ifelse(QCFS_utlc >UTLCTargetFlow, "Gtarget", "Ltarget")) %>% 
  mutate(Qm3s_utlc = QCFS_utlc /35.31467,
         Qm3.15min_utlc = (Qm3s_utlc *60*15),# m3/s to m3/m to m3/15, which is the time b/w samples
         SS_mg.15minP_utlc = SS_mgLp/1000 * Qm3.15min_utlc,#convert mg/L to mg/m3; 
         SRP_mg.15minP_utlc = SRP_mgLp/1000 * Qm3.15min_utlc) 

  
ggplot(UTLCfin, aes(y = SS_mgLp, SS_mgL)) + geom_point()
ggplot(UTLCfin, aes(y = SRP_mgLp, SRP_mgL)) + geom_point()

# include downstream sites ----
# First prepare MR data
qMRdef2 <-  qMRdef %>% 
  select(dateTime, QCFS_mrdef = Flow_Inst) %>% 
  mutate(Qm3s_mrdef = QCFS_mrdef /35.31467) %>% 
  select(dateTime, Qm3s_mrdef)

qMRwat2 <- qMRwat  %>% 
  select(dateTime, QCFS_mrwat = Flow_Inst) %>% 
  mutate(Qm3s_mrwat = QCFS_mrwat /35.31467) %>% 
  select(dateTime, Qm3s_mrwat)


UTLCfin2 <-  UTLCfin %>% 
  left_join(qMRdef2, by = "dateTime") %>% 
  left_join(qMRwat2, by = "dateTime")


# Export ----
# write.csv(UTLCfin2, "04a_generatedDataOnGit/02d_UTLC.csv")
# save.image("03_Rdata/02a_CombUTLC_Rdat")
# load("03_Rdata/02a_CombUTLC_Rdat")



