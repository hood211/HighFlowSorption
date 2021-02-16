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
# load data
#############################################################################################
# load(".02_Turk_Rdat")

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
qTurk <- StreamQGrabFun(siteNo = "04192599",
                        pCode = "00060",
                        start.date = "2014-12-17",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qTurk, aes(y = Flow_Inst, x = dateTime)) +geom_line()

#############################################################################################
#Water quality data
#############################################################################################

wqTurk <- read.csv("01_RawData/00wq_turkeydata.csv") %>% 
  mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M", tz = "America/New_York"))%>% 
    select(-Days741001, -Future, -Month)



# just take a look
ggplot(wqTurk, aes(y = QCFS, x = dateTime)) +geom_line()


#############################################################################################
# Combine water quality and Q data
#############################################################################################

#round data to the nearest 15 min interval
wqTurk.15min <-  wqTurk %>% 
  mutate(dateTime.15mins = round_date(dateTime, unit = "15 minutes")) %>% 
  arrange(dateTime.15mins)

Turk <- qTurk %>% 
  full_join(wqTurk.15min, by = c("dateTime" = "dateTime.15mins")) %>% 
  select(-agency_cd, -dateTime.y, -TP_mgL, -c(NO23_mgL:Cond_umho)) %>% 
  filter(dateTime >= as.POSIXct("2008-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>% #big gap in Q data, starting after that
  filter(dateTime < as.POSIXct("2019-07-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(dateTime) %>% 
  mutate(SS_mgL = ifelse(SS_mgL <= 0, 0, SS_mgL),
         SRP_mgL = ifelse(SRP_mgL <= 0, 0, SRP_mgL)) %>% 
  mutate(Y = as.numeric(as.character(strftime(dateTime, format = "%Y"))),
         M = as.numeric(as.character(strftime(dateTime, format = "%m"))),
         DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))))

# Check alignment of samples
ggplot() +
  geom_point(data = Turk, aes(y = Flow_Inst, x = DOY)) +
  geom_point(data = Turk, aes(y = QCFS, x = DOY), color = "red") +
  facet_grid(Y ~., scales = "free_y")

# There's some missing chunks - don't know what I can do about that
ggplot(Turk, aes(y = Flow_Inst, x = DOY))+
  geom_line(size = 0.5) +
  facet_grid(Y~.)

#############################################################################################
# Determine the 75% flow
#############################################################################################

TurkTargetFlow <- quantile(Turk$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; TurkTargetFlow

ggplot(Turk %>% 
         mutate(AboveTargetFlow = ifelse(Flow_Inst > TurkTargetFlow, "above", "below")), aes(y = log(Flow_Inst), x = DOY, color = AboveTargetFlow)) +
  geom_point(size = 0.5) +
  facet_grid(Y~.)

TurkgTar <- Turk %>% 
  filter(Flow_Inst > TurkTargetFlow) %>% 
  mutate(lSRP_mgL = log(SRP_mgL + 0.0001),
         lSS_mgL = log(SS_mgL + 0.0001),
         lFlow_Inst = log(Flow_Inst)) %>% 
  filter(SS_mgL > 0.01) %>% 
  filter(SRP_mgL > 0.01)
  
#############################################################################################
# predict missing SRP data
#############################################################################################
ggplot(TurkgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(Y))) +
  geom_point()

ggplot(TurkgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(M))) +
  geom_point()

# models to predict missing data
# gonna go with this one - just simplier
gm_SRP_Turk0 <- gamm(lSRP_mgL ~ s(lFlow_Inst), data = TurkgTar); summary(gm_SRP_Turk0$gam)
gm_SRP_Turk <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc", k = 5), data = TurkgTar); summary(gm_SRP_Turk$gam)
# gm_SRP_Turk1 <- gamm(lSRP_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc", k = 5), data = TurkgTar); summary(gm_SRP_Turk1$gam) # not worth evaluating year
# best by R2
gm_SRP_Turk2 <- gamm(lSRP_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = TurkgTar); summary(gm_SRP_Turk2$gam)
# gm_SRP_Turk3 <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc", k = 5) + s(Y, bs = "re"), data = TurkgTar); summary(gm_SRP_Turk3$gam)

# use model to predict values
TurkgTar$lSRP_mgLp <- predict.gam(gm_SRP_Turk$gam, TurkgTar, type = "response")

# predicted v. response
ggplot(TurkgTar, aes(y = lSRP_mgLp, x = lSRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(TurkgTar, aes(y = exp(lSRP_mgLp), x = SRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

#############################################################################################
# predict missing SS data
#############################################################################################
ggplot(TurkgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(Y))) +
  geom_point()

ggplot(TurkgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(M))) +
  geom_point()

# models to predict missing data
gm_SSTurk0 <- gamm(lSS_mgL ~ s(lFlow_Inst), data = TurkgTar); summary(gm_SSTurk0$gam)
# going to go with this one
gm_SSTurk <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc", k = 5), data = TurkgTar); summary(gm_SSTurk$gam)
# gm_SSTurk1 <- gamm(lSS_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc"), data = TurkgTar); summary(gm_SSTurk1$gam)
# best by R2
gm_SSTurk2 <- gamm(lSS_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = TurkgTar); summary(gm_SSTurk2$gam)
# gm_SSTurk3 <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc") + s(Y, bs = "re"), data = TurkgTar); summary(gm_SSTurk3$gam)


# use model to predict values
TurkgTar$lSS_mgLp <- predict.gam(gm_SSTurk$gam, TurkgTar, type = "response")

# predicted v. response
ggplot(TurkgTar, aes(y = lSS_mgLp, x = lSS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(TurkgTar, aes(y = exp(lSS_mgLp), x = SS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


# back transform
TurkgTar2 <- TurkgTar %>% 
  mutate(SRP_mgLp = exp(lSRP_mgL),
         SS_mgLp = exp(lSS_mgLp))

#############################################################################################
# Recombine
#############################################################################################
Turk2 <- Turk %>% 
  mutate(lFlow_Inst = ifelse(Flow_Inst > 0, log(Flow_Inst), as.numeric("NA")))

# this makes calcs for low flows, those values have not been validated
# MIGHT THINK ABOUT KILLING THIS FOR < TARGET BELOW
Turk2$lSRP_mgLp = predict.gam(gm_SRP_Turk$gam, Turk2, type = "response")
Turk2$lSS_mgLp = predict.gam(gm_SSTurk$gam, Turk2, type = "response")

Turkfin <- Turk2 %>% 
  select(dateTime, M, Flow_Inst, QCFS, SS_mgL, SRP_mgL, lSRP_mgLp, lSS_mgLp) %>% 
  mutate(SRP_mgLp = exp(lSRP_mgLp),
         SS_mgLp = exp(lSS_mgLp)) %>% 
  select(dateTime, QCFS_Turk = Flow_Inst, SS_mgL, SRP_mgL, SS_mgLp, SRP_mgLp) %>% 
  mutate(gTarPerFlow = ifelse(QCFS_Turk >TurkTargetFlow, "Gtarget", "Ltarget")) %>% 
  mutate(Qm3s_Turk = QCFS_Turk /35.31467,
         Qm3.15min_Turk = (Qm3s_Turk *60*15),# m3/s to m3/m to m3/15, which is the time b/w samples
         SS_mg.15minP_Turk = SS_mgLp/1000 * Qm3.15min_Turk,#convert mg/L to mg/m3; 
         SRP_mg.15minP_Turk = SRP_mgLp/1000 * Qm3.15min_Turk) 

ggplot(Turkfin, aes(y = SS_mgLp, SS_mgL, color = gTarPerFlow)) + geom_point()
ggplot(Turkfin, aes(y = SRP_mgLp, SRP_mgL, color = gTarPerFlow)) + geom_point()

#############################################################################################
# Export
#############################################################################################
write.csv(Turkfin, "04_generatedData/02f_Turk.csv")
save.image("03_Rdata/02_Turk_Rdat")




