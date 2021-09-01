library(mgcv)
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(readr)
library(readxl)
library(zoo)
library(openxlsx)
library(discharge)
library(MuMIn)
library(GGally)
library(ggpubr)
library(egg)
library(grid)
#load functions
source("02b_SorptionScaling/zz_TravelTimeEstFuns.R")

#############################################################################################
#Water quality data
#############################################################################################

wqMaum <- read.csv("01_RawData/00wq_maumeedata.csv")[,1:15] %>% 
  mutate(dateTime= as.POSIXct(dateTime, format = "%m/%d/%Y %H:%M")) %>% 
  select(-Days741001, -Future, -Month) %>% 
  mutate(Qm3s = as.numeric(QCFS)/35.31467) %>% 
  mutate(SS_mgL = ifelse(SS_mgL == -9, as.numeric("NA"), SS_mgL),
         TP_mgL = ifelse(TP_mgL == -9, as.numeric("NA"), TP_mgL),
         SRP_mgL = ifelse(SRP_mgL == -9, as.numeric("NA"), SRP_mgL),
         SS_mgL = ifelse(SS_mgL <= 0, 0, SS_mgL),
         TP_mgL = ifelse(TP_mgL <= 0, 0, TP_mgL),
         SRP_mgL = ifelse(SRP_mgL <= 0, 0, SRP_mgL))

wqMaum[wqMaum$dateTime == as.POSIXct("0016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S"),]$dateTime = as.POSIXct("2016-06-27 10:10:00", format= "%Y-%m-%d %H:%M:%S")

# just take a look
library(ggthemes)
ggplot(wqMaum, aes(y = TP_mgL, x = dateTime)) +geom_point() +
  theme_few() +
  ylab("Total phosphorus (mg P/L") +
  xlab("Time")
ggplot(wqMaum, aes(y = QCFS, x = dateTime)) +geom_line()

# cummulative sum
ggplot(wqMaum %>% 
         mutate(Y = as.character(strftime(dateTime, format = "%Y"))), aes(Qm3s, color = Y )) + 
  stat_ecdf(geom = "step") +
  geom_vline(xintercept = MaumTargetFlow)


# how many zeros
dim(wqMaum %>% filter(SS_mgL == 0)) #5
dim(wqMaum %>% filter(TP_mgL == 0)) #0
dim(wqMaum %>% filter(SRP_mgL == 0)) #926

ggplot(wqMaum %>% filter(SRP_mgL == 0), aes(y = SRP_mgL, x = dateTime)) +geom_point()

#############################################################################################
# Discharge DAta
# Need regular sampling periods to determine percentile flow
# THIS STARTS IN 1987
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

# better data, but this is just from 1 OCt 1987 to 1 Sep 2020
# using this because the data is more regular than Heidelberg data used later
# this is similar to what is on the guaging station site which has 87 water years = 7500 cfs
MaumTargetFlow <- quantile(qMaum$Flow_Inst, probs = 0.75, na.rm = TRUE)[[1]]/35.31467 ; MaumTargetFlow
hist(log10(qMaum$Flow_Inst))


#############################################################################################
# ESTIMATE IMPACT OF SORPTION
#############################################################################################

# get P sorbed data
# get P sorbed/gDM above flow target
StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")

WC <- read.csv("04_generatedData/05d_WCboot.csv") 
UTLC <- read.csv("04_generatedData/05d_UTLCboot.csv")
SFT <- read.csv("04_generatedData/05d_SFTboot.csv") 



WC2 <- WC %>%   
  filter(gTarPerFlow == "Gtarget") %>% 
  select(dateTime, gPsorb_gSS) %>% 
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
  filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
  mutate(stream = "WC")



UTLC2 <- UTLC %>%   
  filter(gTarPerFlow == "Gtarget") %>% 
  select(dateTime, gPsorb_gSS) %>% 
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
  filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
  mutate(stream = "UTLC")



SFT2 <- SFT %>%   
  filter(gTarPerFlow == "Gtarget") %>% 
  select(dateTime, gPsorb_gSS) %>% 
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
  filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
  mutate(stream = "SFT")

gPsrbGdmData <- rbind(WC2, UTLC2, SFT2)


hist(gPsrbGdmData$gPsorb_gSS)

#####################
# Visualize
#####################
# WC3 <- WC %>% 
#   select(dateTime, gTarPerFlow, DRPloadUS.gPwin, DRPloadUS.gPwin, Svol_gPwindow,Qm3m = Qm3s_West, bootN) %>% 
#   mutate(stream = "WC")
# 
# SFT3 <- SFT %>% 
#   select(dateTime, gTarPerFlow,  DRPloadUS.gPwin, DRPloadUS.gPwin, Svol_gPwindow,Qm3m,  bootN) %>% 
#   mutate(stream = "STF")
# 
# UTLC3 <- UTLC %>% 
#   select(dateTime, gTarPerFlow,  DRPloadUS.gPwin, DRPloadUS.gPwin, Svol_gPwindow,Qm3m = Qm3s_utlc,   bootN) %>% 
#   mutate(stream = "UTLC")
# 
# sorbDat <- rbind(WC3, SFT3, UTLC3)
# 
# 
# StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
# EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")
# sorbDat.Shf <- sorbDat %>% 
#   # subset to stumpf window
#   filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
#   mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
#          Svol_gPwindow.hf = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0)) %>% 
#   select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Svol_gPwindow.hf, Qm3m, bootN, stream) %>% 
#   group_by(bootN, stream) %>% 
#   summarise_at(vars(DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Svol_gPwindow.hf), list(sum = sum), na.rm = TRUE) %>% 
#   mutate(DRPsorbed.hf = Svol_gPwindow.hf_sum/DRPloadUS.gPwin.hf_sum,
#          DRPsorbed.allF = Svol_gPwindow.hf_sum/DRPloadUS.gPwin_sum) 
# 
# ggplot(sorbDat.Shf) +
#   geom_density(aes(x = DRPsorbed.hf, fill = stream))+
#   geom_rug(aes(x = DRPsorbed.hf, y = 0), position = position_jitter(height = 0)) +
#   xlab("Portion DRP sorbed at Waterville (> 75% flows)") +
#   xlim(0,1)
#####################
# BOOTSTRAP
#####################
  BootNum <- 500
  streamWQ_ibs <- as.data.frame(matrix(nrow = 1, ncol = 14)) %>% 
    mutate(bootN = as.numeric("NA"))
  names(streamWQ_ibs) <- c("dateTime", "SmpTimeWindowDay", "SS_mgL", "TP_mgL", "SRP_mgL", "Qm3s",
                           "gPsorbedgDM", "Qm3window", "gDRPwindow", "gSSwindow", "gTPwindow", "gPsorbWindowP", "gPsorbWindow", "gDRPwindowWOsorp", "bootN")
  
  # bootstrap uses a random value of EPC, Scap, and msS > target flow
  # This allows for selecting a different value for each (ignoring sample data)
  # Maximizes variation
  for(i in 1:BootNum) {
    SorpSmpN <- seq(1,length(gPsrbGdmData$gPsorb_gSS))
    gPsorbedgDMe_i = gPsrbGdmData$gPsorb_gSS[sample(SorpSmpN, size = dim(wqMaum)[1], replace = TRUE)]
    streamWQ_i <-  MaumFun(WQDat = wqMaum,
                           gPsorbedgDMe = gPsorbedgDMe_i) %>%
      mutate(bootN = i)
    streamWQ_ibs <- rbind(streamWQ_ibs, streamWQ_i)
  }
  
  
  wqMaum2 <- streamWQ_ibs %>% 
    mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC"),
           DOY = strftime(dateTime, format = "%j"),
           M = as.numeric(as.character(strftime(dateTime, format = "%m"))),
           Y = as.factor(as.character(strftime(dateTime, format = "%Y")))) %>% 
    filter(M >= 3 & M < 7) %>% 
    group_by(Y, DOY, bootN) %>% 
    summarize_at(vars(gPsorbedgDM:gDRPwindowWOsorp), list("s" = sum), na.rm = TRUE) %>% 
    group_by(Y, bootN) %>% 
    summarize_at(vars(gPsorbedgDM_s:gDRPwindowWOsorp_s), list("s" = sum), na.rm = TRUE) %>% 
    rename(gPsorbedgDM = gPsorbedgDM_s_s, # THIS NEEDS TO BE DELETED - DOESN'T MAKE SENSE TO SUM
           Qm3marJun = Qm3window_s_s, gDRPmarJun= gDRPwindow_s_s, gSSmarJun = gSSwindow_s_s,
           gTPmarJun = gTPwindow_s_s,
           gPsorbMarJun = gPsorbWindow_s_s, 
           gPsorbMarJunP = gPsorbWindowP_s_s, gDRMarJunWOsorp= gDRPwindowWOsorp_s_s) %>% 
    group_by(Y) %>% 
    mutate(perSorb1 = gPsorbMarJun/gDRPmarJun,
           perSorb2 = gPsorbMarJunP/gDRPmarJun,
           DRPMarJunWvWOsorp = (gPsorbMarJunP + gDRPmarJun)/gDRPmarJun) %>% 
    summarize_at(vars(gPsorbedgDM:DRPMarJunWvWOsorp), funs(!!!p_funs)) %>% 
    mutate(Yn = as.numeric(as.character(Y)),
           Yg = ifelse(Yn < 1980, "70s",
                  ifelse(Yn >= 1980 & Yn < 1990, "80s",
                     ifelse(Yn >= 1990 & Yn < 2000, "90s",
                       ifelse(Yn >= 2000 & Yn < 2010, "00s",
                          ifelse(Yn >= 2010, "10s", "blah")))))) %>% 
    mutate(Yg = as.factor(Yg),
           Yg = fct_relevel(Yg, c("70s","80s", "90s", "00s", "10s")))
  

  
  # cyano index
  wqMaum03 <- wqMaum2 %>% 
    mutate(mtDRPmarJun_50per = gDRPmarJun_50per/1e6,
           mtPsorbMarJunP_50per  = gPsorbMarJunP_50per/1e6,
           mtDRPmarJun_50perWoS = mtDRPmarJun_50per + mtPsorbMarJunP_50per,
           cyanoIndex = 0.48 * 10^(mtDRPmarJun_50per * 0.00387), # says a * 10^-3 = 0.001?
           cyanoIndexWoS = 0.48 * 10^(mtDRPmarJun_50perWoS * 0.00387))%>% 
    filter(Yn >= 2003) %>% 
    select(Yn, mtDRPmarJun_50per:cyanoIndexWoS) %>% 
    pivot_longer(cols = c(mtDRPmarJun_50per:cyanoIndexWoS), names_to = "var", values_to = "values") 
  
  # write.csv(wqMaum2, "04_generatedData/06d_Scale2WVall.csv")
  # write.csv(wqMaum03, "04_generatedData/06d_Scale2WVhabs.csv")
  
  
  ################  
  # % P loads during high flows
  ################  
  MaumLoads <- MaumFun(WQDat = wqMaum,
          gPsorbedgDMe = gPsorbedgDMe_i)
  
  ggplot(MaumLoads %>% 
           mutate(Y = as.character(strftime(dateTime, format = "%Y")),
                  DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))), 
                  aes(y = gDRPwindow, x = DOY)) +
    geom_line() +
    facet_wrap(vars(Y))
  
  ggplot(wqMaum %>% 
           mutate(Y = as.character(strftime(dateTime, format = "%Y")),
                  DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))), 
         aes(y = Qm3s, x = DOY)) +
    geom_line() +
    facet_wrap(vars(Y))

 
  MaumLoadsAn0 <- MaumLoads %>% 
    mutate(TargetFlow = as.factor(ifelse(Qm3s > MaumTargetFlow, "G75th", "L75th"))) %>% 
    filter(!is.na(TargetFlow)) %>% 
    mutate(Y = as.character(strftime(dateTime, format = "%Y")),
           M = as.factor(as.character(strftime(dateTime, format = "%m"))),
           TargetMonths = as.factor(ifelse(M %in% c("03", "04", "05", "06"), "MarJun",
                                           ifelse(M == "07", "Jul","AugFeb")))) %>% 
    # remove 1978, 79, 80,81; could leave 78 for spring stuff
    filter(Y != 1978 & Y != 1979 & Y != 1980 & Y != 1981) %>% 
    select(dateTime, SmpTimeWindowDay, gDRPwindow, gSSwindow, TargetFlow,Y, M, TargetMonths) %>% 
    # select(-c(SS_mgL:Qm3window), -c(gSSwindow:gDRPwindowWOsorp), -SmpTimeWindowDay) %>% 
    group_by(Y, TargetMonths, TargetFlow) %>% 
    summarise(gDRPwindow= sum(gDRPwindow, na.rm = TRUE),
              gSSwindow = sum(gSSwindow, na.rm = TRUE))  
  
  MaumLoadsAn0b <- MaumLoadsAn0 %>% 
    mutate(key = as.factor(paste0(TargetMonths,"_",TargetFlow))) %>% 
    mutate(key2 = fct_rev(fct_relevel(key, c("AugFeb_L75th","MarJun_L75th","Jul_L75th","AugFeb_G75th", "MarJun_G75th", "Jul_G75th"))))
  
# sum of loads for HABsForecast talk
  MaumLoadsAn0bS <- MaumLoadsAn0b %>% 
    mutate(TargetMonths2 = fct_recode(TargetMonths, MarJul = "MarJun", MarJul = "Jul")) %>%
    filter(Y >= 2003) %>%
    group_by(Y, TargetFlow, TargetMonths2) %>% 
    summarize(across(gDRPwindow:gSSwindow, sum)) %>%
    group_by(TargetFlow, TargetMonths2) %>% 
    summarize(across(gDRPwindow:gSSwindow, mean)) %>%
    pivot_wider(id_cols = "TargetMonths2", names_from = "TargetFlow", values_from = c("gDRPwindow", "gSSwindow")) %>% 
    mutate(gDRPwin = gDRPwindow_G75th + gDRPwindow_L75th,
           gSSwin = gSSwindow_G75th + gSSwindow_L75th,
           DRP_perTotG75th = gDRPwindow_G75th/gDRPwin,
           SS_perTot75th = gSSwindow_G75th/gSSwin)
  

  
FigS1a <-   ggplot(MaumLoadsAn0b, 
         aes(y = gDRPwindow/1e6, x = as.numeric(Y), fill = key2)) +
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("pink", "firebrick1", "firebrick4","lightblue", "steelblue3", "steelblue4"), name = "Loading window",
                      labels = c("<25% flow ex.: Jul",
                                 "<25% flow ex.: Mar-Jun",
                                 "<25% flow ex.: Aug-Feb",
                                 ">25% flow ex.: Jul",
                                 ">25% flow ex.: Mar-Jun",
                                 ">25% flow ex.: Aug-Feb"))+
    theme_bw() +
    ylab("DRP load (ton P)") +
    xlab("Year") +
    theme(legend.position = c(0.12,0.79),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 12, face = "bold"),
          legend.key.height = unit(0.45,"cm"),
          legend.spacing.y = unit(0.05,"cm"),
          legend.margin = margin(0.1,0,0,0, unit="cm"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) 
  

  
FigS1b <-   ggplot(MaumLoadsAn0b, 
         aes(y = gSSwindow/1e12, x = as.numeric(Y), fill = key2)) +
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("pink", "firebrick1", "firebrick4","lightblue", "steelblue3", "steelblue4"), name = "Loading window",
                      labels = c("<25% flow ex.: Jul",
                                 "<25% flow ex.: Mar-Jun",
                                 "<25% flow ex.: Aug-Feb",
                                 ">25% flow ex.: Jul",
                                 ">25% flow ex.: Mar-Jun",
                                 ">25% flow ex.: Aug-Feb"))+
    theme_bw() +
    ylab(expression(paste("SS load (ton dry mass x ",10^-6,")"))) +
    xlab("Year") +
    theme(legend.position = "none",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          legend.key.height = unit(0.5,"cm"),
          legend.spacing.y = unit(0.05,"cm"),
          legend.margin = margin(0.1,0,0,0, unit="cm"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) 
  
FigS1a.g <- ggplotGrob(FigS1a)
FigS1b.g <- ggplotGrob(FigS1b)

FigS1a.gtf <- gtable_frame(FigS1a.g, width = unit(1, "null"), height = unit(0.75, "null"))
FigS1b.gtf <- gtable_frame(FigS1b.g, width = unit(1, "null"), height = unit(0.75, "null"))
FigS1.gtf <- gtable_frame(gtable_rbind(FigS1a.gtf, FigS1b.gtf),
                          width = unit(1, "null"), height = unit(1.5, "null"))

png("05_Figures/06aFigS1.png", units="in", width=8, height=8, res=300)
grid.newpage()
grid.draw(FigS1.gtf)
grid.text("a", x = unit(0.02,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.02,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()






# Calc percent loads during different windows  
  MaumLoadsAn_DRP <- MaumLoadsAn0 %>% 
      pivot_wider(id_cols = Y, names_from = c(TargetFlow, TargetMonths), values_from = gDRPwindow) %>% 
    mutate(G75th_Jul = ifelse(is.na(G75th_Jul), 0, G75th_Jul)) %>% 
    mutate(gDRPannual = G75th_AugFeb + L75th_AugFeb + G75th_MarJun + L75th_MarJun + G75th_Jul + L75th_Jul,
           gDRP_MarJuly = G75th_MarJun + L75th_MarJun + G75th_Jul + L75th_Jul,
           gDRP_MarJulyG75th = G75th_MarJun + G75th_Jul,
           perMarJulyG75th_gDRP_MarJuly = gDRP_MarJulyG75th/gDRP_MarJuly,
           perG75thJul_MarJulyG75th = G75th_Jul/gDRP_MarJulyG75th)

  MaumLoadsAn_DRP %>% 
    ungroup() %>% 
    summarise(mean = mean(perG75thJul_MarJulyG75th, na.rm = T)*100,
              sd = sd(perG75thJul_MarJulyG75th, na.rm = T)*100)
  
  MaumLoadsAn_DRP %>% 
    ungroup() %>% 
    summarise(mean = mean(perMarJulyG75th_gDRP_MarJuly, na.rm = T)*100,
              sd = sd(perMarJulyG75th_gDRP_MarJuly, na.rm = T)*100)
    
  MaumLoadsAn_SS <- MaumLoadsAn0 %>% 
    pivot_wider(id_cols = Y, names_from = c(TargetFlow, TargetMonths), values_from = gSSwindow) %>% 
    mutate(G75th_Jul = ifelse(is.na(G75th_Jul), 0, G75th_Jul)) %>% 
    mutate(gDRPannual = G75th_AugFeb + L75th_AugFeb + G75th_MarJun + L75th_MarJun + G75th_Jul + L75th_Jul,
           gDRP_MarJuly = G75th_MarJun + L75th_MarJun + G75th_Jul + L75th_Jul,
           gDRP_MarJulyG75th = G75th_MarJun + G75th_Jul,
           perMarJulyG75th_gDRP_MarJuly = gDRP_MarJulyG75th/gDRP_MarJuly,
           perG75thJul_MarJulyG75th = G75th_Jul/gDRP_MarJulyG75th)
  
  MaumLoadsAn_SS %>% 
    ungroup() %>% 
    summarise(mean = mean(perG75thJul_MarJulyG75th, na.rm = T)*100,
              sd = sd(perG75thJul_MarJulyG75th, na.rm = T)*100)
  
  MaumLoadsAn_SS %>% 
    ungroup() %>% 
    summarise(mean = mean(perMarJulyG75th_gDRP_MarJuly, na.rm = T)*100,
              sd = sd(perMarJulyG75th_gDRP_MarJuly, na.rm = T)*100)

  
  
    

  

  
  
  ################  
  save.image("03_Rdata/06a_ScaleAtWaterville_Rdat")
  # load("03_Rdata/06a_ScaleAtWaterville_Rdat")
  ################
  
