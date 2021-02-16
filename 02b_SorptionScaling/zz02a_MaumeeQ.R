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
# load("")

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
qMaum <- StreamQGrabFun(siteNo = "04193500",
                        pCode = "00060",
                        start.date = "1975-10-01",
                        end.date = "2020-09-01") %>% 
  mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 

# looks good
ggplot(qMaum, aes(y = Flow_Inst, x = dateTime)) +geom_point()

#############################################################################################
#Water quality data
#############################################################################################

wqMaum <- read.csv("01_data/00wq_maumeedata.csv")[,1:15] %>% 
  mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M")) %>% 
    select(-Days741001, -Future, -Month) %>% 
  mutate(QCFS = as.numeric(QCFS))

wqMaum[wqMaum$dateTime == as.POSIXct("0016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S"),]$dateTime = as.POSIXct("2016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S")
 
# dateTime = ifelse(dateTime == as.POSIXct("0016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S"),
#                   as.POSIXct("2016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S"),
#                   dateTime

# just take a look
ggplot(wqMaum, aes(y = QCFS, x = dateTime)) +geom_line()


#############################################################################################
# Combine water quality and Q data
#############################################################################################

#round data to the nearest 15 min interval
wqMaum.15min <-  wqMaum %>% 
  mutate(dateTime.15mins = round_date(dateTime, unit = "15 minutes")) %>% 
  arrange(dateTime.15mins)

Maum <- qMaum %>% 
  full_join(wqMaum.15min, by = c("dateTime" = "dateTime.15mins")) %>% 
  select(-agency_cd, -dateTime.y, -c(NO23_mgL:Cond_umho)) %>% 
  filter(dateTime >= as.POSIXct("1987-10-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>% #big gap in Q data, starting after that
  filter(dateTime < as.POSIXct("2019-07-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(dateTime) %>% 
  mutate(SS_mgL = ifelse(SS_mgL <= 0, 0, SS_mgL),
         SRP_mgL = ifelse(SRP_mgL <= 0, 0, SRP_mgL),
         TP_mgL = ifelse(TP_mgL <= 0, 0, TP_mgL)) %>% 
  mutate(Y = as.numeric(as.character(strftime(dateTime, format = "%Y"))),
         M = as.numeric(as.character(strftime(dateTime, format = "%m"))),
         DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))))

# Check alignment of samples
ggplot() +
  geom_point(data = Maum, aes(y = Flow_Inst, x = DOY)) +
  geom_point(data = Maum, aes(y = QCFS, x = DOY), color = "red") +
  facet_grid(Y ~., scales = "free_y")



#############################################################################################
# Determine the 75% flow
#############################################################################################

MaumTargetFlow <- quantile(Maum$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]; MaumTargetFlow

ggplot(Maum %>% 
         mutate(AboveTargetFlow = ifelse(Flow_Inst > MaumTargetFlow, "above", "below")), aes(y = log(Flow_Inst), x = DOY, color = AboveTargetFlow)) +
  geom_point(size = 0.5) +
  facet_grid(Y~.)

StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")

  # subset to stumpf window
  
MaumgTar <- Maum %>% 
  filter(M >= 3 & M <= 7) %>% 
  filter(Flow_Inst > MaumTargetFlow) %>% 
  mutate(lSRP_mgL = log(SRP_mgL + 0.0001),
         lSS_mgL = log(SS_mgL + 0.0001),
         lFlow_Inst = log(Flow_Inst)) %>% 
  filter(SS_mgL > 0.01) %>% 
  filter(SRP_mgL > 0.01) %>% 
  mutate(Mf = as.factor(as.character(M)))
  
#############################################################################################
# predict missing SRP data
#############################################################################################
ggplot(MaumgTar, aes(y = SRP_mgL, x = Flow_Inst, color= Mf)) +
  geom_point() +
  facet_wrap(vars(Y),ncol = 5)

ggplot(MaumgTar, aes(y = log(Flow_Inst), x = log(QCFS), color= Mf)) +
  geom_point() 

ggplot(MaumgTar, aes(y = log(SRP_mgL), x = log(Flow_Inst), color=as.factor(M))) +
  geom_point()

# models to predict missing data
gm_SRP_Maum0 <- gamm(lSRP_mgL ~ s(lFlow_Inst), data = MaumgTar); summary(gm_SRP_Maum0$gam)
# best by R2
gm_SRP_Maum <- gamm(lSRP_mgL ~ s(lFlow_Inst, by = Mf), data = MaumgTar); summary(gm_SRP_Maum$gam)

# best model by R2
gm_SRP_Maum1 <- gamm(lSRP_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc", k = 5), data = MaumgTar); summary(gm_SRP_Maum1$gam) 
gm_SRP_Maum2 <- gamm(lSRP_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = MaumgTar); summary(gm_SRP_Maum2$gam)
gm_SRP_Maum3 <- gamm(lSRP_mgL ~ s(lFlow_Inst) + s(M, bs = "cc", k = 5) + s(Y, bs = "re"), data = MaumgTar); summary(gm_SRP_Maum3$gam)

# use model to predict values
MaumgTar$lSRP_mgLp <- predict.gam(gm_SRP_Maum$gam, MaumgTar, type = "response")

# predicted v. response
ggplot(MaumgTar, aes(y = lSRP_mgLp, x = lSRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(MaumgTar, aes(y = exp(lSRP_mgLp), x = SRP_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

#############################################################################################
# predict missing SS data
#############################################################################################
ggplot(MaumgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(Y))) +
  geom_point()

ggplot(MaumgTar, aes(y = lSS_mgL, x = lFlow_Inst, color=as.factor(M))) +
  geom_point()

# models to predict missing data
gm_SSMaum0 <- gamm(lSS_mgL ~ s(lFlow_Inst), data = MaumgTar); summary(gm_SSMaum0$gam)
# going to go with this one
gm_SSMaum <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc"), data = MaumgTar); summary(gm_SSMaum$gam)
gm_SSMaum1 <- gamm(lSS_mgL ~ s(Flow_Inst, by = Y) + s(M, bs = "cc"), data = MaumgTar); summary(gm_SSMaum1$gam)
gm_SSMaum2 <- gamm(lSS_mgL ~ te(Flow_Inst, M, bs = c("cr","cc")), data = MaumgTar); summary(gm_SSMaum2$gam)
gm_SSMaum3 <- gamm(lSS_mgL ~ s(lFlow_Inst) + s(M, bs = "cc") + s(Y, bs = "re"), data = MaumgTar); summary(gm_SSMaum3$gam)


# use model to predict values
MaumgTar$lSS_mgLp <- predict.gam(gm_SSMaum$gam, MaumgTar, type = "response")

# predicted v. response
ggplot(MaumgTar, aes(y = lSS_mgLp, x = lSS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot(MaumgTar, aes(y = exp(lSS_mgLp), x = SS_mgL)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


# back transform
MaumgTar2 <- MaumgTar %>% 
  mutate(SRP_mgLp = exp(lSRP_mgL),
         SS_mgLp = exp(lSS_mgLp))

#############################################################################################
# Recombine
#############################################################################################
Maumfin <- Maum %>% 
  select(dateTime, Flow_Inst, QCFS, SS_mgL, SRP_mgL) %>% 
  left_join(MaumgTar2 %>% 
              select(dateTime, SRP_mgLp, SS_mgLp), by = "dateTime") %>% 
  select(dateTime, QCFS_Maum = Flow_Inst, SS_mgL, SRP_mgL, SS_mgLp, SRP_mgLp) %>% 
  mutate(gTarPerFlow = ifelse(QCFS_Maum >MaumTargetFlow, "Gtarget", "Ltarget")) %>% 
  mutate(Qm3s_Maum = QCFS_Maum /35.31467,
         Qm3.15min_Maum = (Qm3s_Maum *60*15),# m3/s to m3/m to m3/15, which is the time b/w samples
         SS_mg.15minP_Maum = SS_mgLp/1000 * Qm3.15min_Maum,#convert mg/L to mg/m3; 
         SRP_mg.15minP_Maum = SRP_mgLp/1000 * Qm3.15min_Maum) 


#############################################################################################
# Export
#############################################################################################
write.csv(Maumfin, "02f_Maum.csv")
save.image(".02_Maum_Rdat")




