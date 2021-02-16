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
# load data
#############################################################################################
# load(".02_West_Rdat")

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

# UNNAMED TRIB TO LOST CREEK -> MR_def -> MR_wat
# https://waterdata.usgs.gov/monitoring-location/04185440/#parameterCode=00060
qWest <- StreamQGrabFun(siteNo = "04192574",
                        pCode = "00060",
                        start.date = "2014-12-30",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qWest, aes(y = Flow_Inst, x = dateTime)) +geom_point()

#############################################################################################
#Water quality data
#############################################################################################

wqWest <- read.csv("01_RawData/00wq_Westdata.csv") %>% 
  mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M"))%>% 
    select(-Days741001, -Future, -Month)



# just take a look
ggplot(wqWest, aes(y = QCFS, x = dateTime)) +geom_line()


#############################################################################################
# Combine water quality and Q data
#############################################################################################

#round data to the nearest 15 min interval
wqWest.15min <-  wqWest %>% 
  mutate(dateTime.15mins = round_date(dateTime, unit = "15 minutes")) %>% 
  arrange(dateTime.15mins)

West <- qWest %>% 
  full_join(wqWest.15min, by = c("dateTime" = "dateTime.15mins")) %>% 
  select(-agency_cd, -dateTime.y, -TP_mgL, -c(NO23_mgL:Cond_umho)) %>% 
  filter(dateTime >= as.POSIXct("2015-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>% #big gap in Q data, starting after that
  filter(dateTime < as.POSIXct("2019-07-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(dateTime) %>% 
  mutate(SS_mgL = ifelse(SS_mgL <= 0, 0, SS_mgL),
         SRP_mgL = ifelse(SRP_mgL <= 0, 0, SRP_mgL)) %>% 
  mutate(Y = as.numeric(as.character(strftime(dateTime, format = "%Y"))),
         M = as.numeric(as.character(strftime(dateTime, format = "%m"))),
         DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))))

# Check alignment of samples
ggplot() +
  geom_point(data = West, aes(y = Flow_Inst, x = DOY)) +
  geom_point(data = West, aes(y = QCFS, x = DOY), color = "red") +
  facet_grid(Y ~., scales = "free_y")

# There's some missing chunks - don't know what I can do about that
ggplot(West, aes(y = Flow_Inst, x = DOY))+
  geom_point(size = 0.5) +
  facet_grid(Y~.)

#############################################################################################
# Determine the 75% flow
#############################################################################################

WestTargetFlow <- quantile(West$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; WestTargetFlow

ggplot(West %>% 
         mutate(AboveTargetFlow = ifelse(Flow_Inst > WestTargetFlow, "above", "below")), aes(y = log(Flow_Inst), x = DOY, color = AboveTargetFlow)) +
  geom_point(size = 0.5) +
  facet_grid(Y~.)

WestgTar <- West %>% 
  filter(Flow_Inst > WestTargetFlow) %>% 
  mutate(lSRP_mgL = log(SRP_mgL + 0.0001),
         lSS_mgL = log(SS_mgL + 0.0001),
         lFlow_Inst = log(Flow_Inst)) %>% 
  filter(SS_mgL > 0.01) %>% 
  filter(SRP_mgL > 0.01)
  
#############################################################################################
# predict missing SRP data
#############################################################################################
ggplot(WestgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(Y))) +
  geom_point()

ggplot(WestgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(M))) +
  geom_point()

# models to predict missing data
gm_SRP_West0 <- gamm(lSRP_mgL ~ s(lFlow_Inst), data = WestgTar); summary(gm_SRP_West0$gam)
# best by R2
gm_SRP_West <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc"), data = WestgTar); summary(gm_SRP_West$gam)
# gm_SRP_West1 <- gamm(lSRP_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc", k = 5), data = WestgTar); summary(gm_SRP_West1$gam) # not worth evaluating year
gm_SRP_West2 <- gamm(lSRP_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = WestgTar); summary(gm_SRP_West2$gam)
# gm_SRP_West3 <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc", k = 5) + s(Y, bs = "re"), data = WestgTar); summary(gm_SRP_West3$gam)

# use model to predict values
WestgTar$lSRP_mgLp <- predict.gam(gm_SRP_West$gam, WestgTar, type = "response")

# predicted v. response
ggplot(WestgTar, aes(y = lSRP_mgLp, x = lSRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(WestgTar, aes(y = exp(lSRP_mgLp), x = SRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

#############################################################################################
# predict missing SS data
#############################################################################################
ggplot(WestgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(Y))) +
  geom_point()

ggplot(WestgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(M))) +
  geom_point()

# models to predict missing data
gm_SSWest0 <- gamm(lSS_mgL ~ s(lFlow_Inst), data = WestgTar); summary(gm_SSWest0$gam)
# best by R2 going to go with this one
gm_SSWest <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc"), data = WestgTar); summary(gm_SSWest$gam)
# gm_SSWest1 <- gamm(lSS_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc"), data = WestgTar); summary(gm_SSWest1$gam)
gm_SSWest2 <- gamm(lSS_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = WestgTar); summary(gm_SSWest2$gam)
# gm_SSWest3 <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc") + s(Y, bs = "re"), data = WestgTar); summary(gm_SSWest3$gam)


# use model to predict values
WestgTar$lSS_mgLp <- predict.gam(gm_SSWest$gam, WestgTar, type = "response")

# predicted v. response
ggplot(WestgTar, aes(y = lSS_mgLp, x = lSS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(WestgTar, aes(y = exp(lSS_mgLp), x = SS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


# back transform
WestgTar2 <- WestgTar %>% 
  mutate(SRP_mgLp = exp(lSRP_mgL),
         SS_mgLp = exp(lSS_mgLp))

#############################################################################################
# Recombine
#############################################################################################
West2 <- West %>% 
  mutate(lFlow_Inst = ifelse(Flow_Inst > 0, log(Flow_Inst), as.numeric("NA")))


# this makes calcs for low flows, those values have not been validated
# MIGHT THINK ABOUT KILLING THIS FOR < TARGET BELOW
West2$lSRP_mgLp = predict.gam(gm_SRP_West$gam, West2, type = "response")
West2$lSS_mgLp = predict.gam(gm_SSWest$gam, West2, type = "response")


Westfin <- West2 %>% 
  select(dateTime, M, Flow_Inst, QCFS, SS_mgL, SRP_mgL, lSRP_mgLp, lSS_mgLp) %>% 
  mutate(SRP_mgLp = exp(lSRP_mgLp),
         SS_mgLp = exp(lSS_mgLp)) %>% 
  select(dateTime, QCFS_West = Flow_Inst, SS_mgL, SRP_mgL, SS_mgLp, SRP_mgLp) %>% 
  mutate(gTarPerFlow = ifelse(QCFS_West >WestTargetFlow, "Gtarget", "Ltarget")) %>% 
  mutate(Qm3s_West = QCFS_West /35.31467,
         Qm3.15min_West = (Qm3s_West *60*15),# m3/s to m3/m to m3/15, which is the time b/w samples
         SS_mg.15minP_West = SS_mgLp/1000 * Qm3.15min_West,#convert mg/L to mg/m3; 
         SRP_mg.15minP_West = SRP_mgLp/1000 * Qm3.15min_West) 

ggplot(Westfin, aes(y = SS_mgLp, SS_mgL)) + geom_point()
ggplot(Westfin, aes(y = SRP_mgLp, SRP_mgL)) + geom_point()

#############################################################################################
# Export
#############################################################################################
write.csv(Westfin, "04_generatedData/02f_West.csv")
save.image("03_Rdata/02_West_Rdat")




