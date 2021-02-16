library(tidyverse)
# load("05z_scaleUTLC_rdat")

#############################################################################################
# Get data
#############################################################################################
#load functions
source("02b_SorptionScaling/zz_TravelTimeEstFuns.R")

# water quality and Q for target stream and downstream
streamWQ <- read.csv("04_generatedData/02d_UTLC.csv")

# slope and segment length info
segInfo <- read.csv("04_generatedData/03d_slopes.csv")

# sorption data
sorp <- read.csv("04_generatedData/04d_SorpDat.csv") %>% 
  filter(Stream == "UTLC")

# Watershed areas
gsWA <- read.csv("01_RawData/ScalingValTable.csv")


#############################################################################################
# Get data
#############################################################################################

# Here's what we need for scaling

# sorption stuff
  # mass-specific sorption rate, all values for bootstraping
  # EPC - use equation to estimate EPC?????
  # Scap


# for each station
  # mean annual Q at each gauging station along path
  # WA for each station
  # Q at 15 min intervals
  # DRP at 15 min intervals
  # SS at 15 min intervals

# for each reach
  # slope for each reach
  # lg of each reach


#############################################################################################
# Determine mean annual Q
#############################################################################################
# plot to investigate gaps
ggplot(streamWQ %>% 
         select(dateTime, Qm3s_utlc, Qm3s_mrdef, Qm3s_mrwat) %>% 
         pivot_longer(!dateTime, names_to = "stream", values_to = "Q_m3s") %>% 
         mutate(Y = as.factor(as.character(strftime(dateTime, format = "%Y"))),
                DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))), aes(y = Q_m3s, x = DOY, color = stream)) +
  geom_point(size = 0.75) +
  facet_grid(Y ~stream)

MeanAnnQ <- streamWQ %>% 
  select(dateTime, Qm3s_utlc, Qm3s_mrdef, Qm3s_mrwat) %>% 
  pivot_longer(!dateTime, names_to = "stream", values_to = "Q_m3s") %>% 
  mutate(Y = as.factor(as.character(strftime(dateTime, format = "%Y"))),
         DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))) %>% 
  group_by(stream, Y) %>% 
  summarise(Q_m3s = mean(Q_m3s, na.rm = TRUE))

MeanAnnQS <- MeanAnnQ %>% 
  mutate(Y = as.numeric(as.character(Y))) %>% 
  filter(Y >= 2014 & Y < 2019) %>% 
  group_by(stream) %>% 
  summarise(Q_m3s = mean(Q_m3s, na.rm = TRUE))


#############################################################################################
# Scale
#############################################################################################
# Scaling up sorption rate using the median sorption rate
# Input descriptions
# Qdat = data frame with discharge, SS, and DRP data for 2019, see below for units
# SorDat = data frame with sorption stuff
# StreamName = one of Whitney's six focal streams
# StreamQ = name of column in Qdat with Q from gs at top of reach, m3/min
# StreamWA = Watershed area above gs, m2
# ReachLg = legnth of stream from US gs to DS gs, m
# ReachSlope = mean slope of reach from US gs to DS gs, m/m
# ManningsCoef = Mannings Coefficient, UNITS???
# TimeWin = time between gs samples (in Qdat), min
# SS_gDMm3 = SS at US gs, gDM/m3
# tt_pref = one of tt_job, tt_mor, tt_du

# ScalingFun <- function(Qdat, SorDat, StreamName, StreamQ, StreamWA, StreamQan, ReachLg, ReachSlope, ManningsCoef,
#                        TimeWin, SS_gDMm3, tt_pref, StreamDRP){
  # subset sorption data
  SorDatSN <- sorp %>% 
    filter(Stream == "UTLC") %>% 
    filter(Percentile_flow >= 75)
  
  # sorption parameters
  
  Sms_gPgDMmin.Med <- quantile(SorDatSN$Sms_gPgDMmin, probs = 0.5)
  StreamEPC <- quantile(SorDatSN$EPC_gPm3, probs = 0.5)
  Scap_gPgDM <- quantile(SorDatSN$Scap_gPgDM, probs = 0.5)
  
  #no relationship with Q or % flow within stream, so boot strap best approach
  summary(lm(Sms_gPgDMmin ~ Q_m3min, SorDatSN))
  ggplot(SorDatSN, aes(y = log(Sms_gPgDMmin), x = log(Q_m3min))) +geom_point()

  ###################### 
  #scale for reach one between UTLC to  Maumee @ Defiance
  #####################
  # reach parameters
  StreamWA1 <- as.numeric(gsWA[gsWA$USGS.gauge == "Unnamed Trib to Lost Ck",]$WA_m2)
  slope_mm1 <- segInfo[segInfo$Gauge_start == "UTLC",]$SlopeMm
  seglg_m1 <- segInfo[segInfo$Gauge_start == "UTLC",]$SegLgM
  #3-1215: 0.026- Indian Fork below Atwood Dam, near New Cumberland, OH
  # 5-Misc 0.037 - Middle fork of Vermilion River near Danville, Ill. Bed of gravel and small cobbles, banks lined with trees small underbrush
  ManningsCoef <- 0.037/60 # from Barnes report (pretty sure units are metric s/m^(1/3) so converting to min)
  StreamQan1 <- MeanAnnQS[MeanAnnQS$stream == "Qm3s_utlc",]$Q_m3s*60
  
  # other stuff
  TimeWin <- 15 #mins
  
  # Prep data
  streamWQ2 <- streamWQ %>% 
    mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d %H:%M:%S"),
           Yf = as.factor(as.character(strftime(dateTime, format = "%Y"))),
           Y = as.numeric(as.character(Yf)),
           DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))), 
           gTarPerFlow = as.factor(gTarPerFlow)) %>% 
    filter(Y == 2019) %>% 
    mutate(Qm3m = Qm3s_utlc*60, #Q in m3/min
           SSgDMm3 = SS_mgLp*1000/1000,
           DRPgPm3 = SRP_mgLp*1000/1000)
  
# run function
  streamWQatMRDEF <- ScalingFun.step1(WQdat = streamWQ2,
                           Qm3m = Qm3m, 
                           SSgDMm3 = SSgDMm3, 
                           DRPgPm3 = DRPgPm3, 
                           slope_mm = segInfo[segInfo$Gauge_start == "UTLC",]$SlopeMm, 
                           ManningsCoef = 0.037/60,
                           seglg_m = segInfo[segInfo$Gauge_start == "UTLC",]$SegLgM, 
                           Sms_gPgDMmin.Med = quantile(SorDatSN$Sms_gPgDMmin, probs = 0.5), 
                           Scap_gPgDM = quantile(SorDatSN$Scap_gPgDM, probs = 0.5), 
                           StreamEPC = quantile(SorDatSN$EPC_gPm3, probs = 0.5), 
                           TimeWin = 15)
  

  

  # do the travel times seem reasonable?
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = (tt_du)/seglg_m1*100, x = Qm3m)) +
    geom_point() +
    ylab("Mins to travel a 100 m reach")
  
  # How much does sorption influence DRP
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = porDRPloadSorbed, x = DOY)) +
    geom_point() 
  
  # Does sorb get saturated?
  # Always reaching capacity! No reason to go on.
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = Sms_porCap, x = dateTime)) +
    geom_point() +
    ylim(0,5)
  
  ##########################
  # How far does it travel? Seg lg - 51221.81
  # This only goes to next station
  ##########################
  
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = SorbDistM, x = dateTime)) +
    geom_point() 
  
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget")) +
    geom_segment(aes(x = 0, xend = SorbDistM, y = 1, yend = 0.8), position = "jitter", arrow = arrow(), size = 0.25, color = "blue") +
    xlim(0,seglg_m1) +
    ylim(0,1.5)
  
  
  # Why does the portion of DRP sorbed increase?
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = porDRPloadSorbed, x = SSgDMm3, color = DOY)) +
    geom_point() +
    scale_x_log10() +
    scale_x_log10()
  
  ggplot(streamWQatMRDEF %>% 
           filter(gTarPerFlow == "Gtarget"), aes(y = gPsorb_gSS, x = DOY)) +
    geom_point()


  ##########################
  # Portion sorbed during stumpf window - high flows
  ##########################
  StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
  EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")
  streamWQatMRDEF.stfHF <- streamWQatMRDEF %>% 
    # subset to stumpf window
    filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
    mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
           DRPloadDS.gPwin.wS.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadDS.gPwin.wS, 0)) %>% 
    summarise_at(vars(DRPloadUS.gPwin, DRPloadUS.gPwin.hf, DRPloadDS.gPwin.wS.hf), list(sum = sum), na.rm = TRUE)

  UTLC_PercentSorbStmp = streamWQatMRDEF.stfHF$DRPloadDS.gPwin.wS.hf_sum/ streamWQatMRDEF.stfHF$DRPloadUS.gPwin.hf_sum

  #####################
  # BOOTSTRAP
  #####################
  BootNum <- 500
  streamWQ_ibs <- as.data.frame(matrix(nrow = 1, ncol = dim(streamWQatMRDEF)[2])) %>% 
                    mutate(bootN = as.numeric("NA"))
  names(streamWQ_ibs) <- c(names(streamWQatMRDEF), "bootN")
  
  # bootstrap uses a random value of EPC, Scap, and msS > target flow
  # This allows for selecting a different value for each (ignoring sample data)
  # Maximizes variation
  for(i in 1:BootNum) {
    SorpSmpN <- seq(1,length(SorDatSN$Sms_gPgDMmin))
    Sms_gPgDMmin.Med_i = SorDatSN$Sms_gPgDMmin[sample(SorpSmpN, size = dim(streamWQatMRDEF)[1], replace = TRUE)]
    Scap_gPgDM_i = SorDatSN$Scap_gPgDM[sample(SorpSmpN, size = dim(streamWQatMRDEF)[1], replace = TRUE)] 
    StreamEPC_i = SorDatSN$EPC_gPm3[sample(SorpSmpN, size =dim(streamWQatMRDEF)[1], replace = TRUE)]
    streamWQ2legnth <- dim(streamWQ2)[1]
    streamWQ_i <- ScalingFun.step1(WQdat = streamWQ2,
                                   Qm3m = Qm3m, 
                                   SSgDMm3 = SSgDMm3, 
                                   DRPgPm3 = DRPgPm3, 
                                   slope_mm = segInfo[segInfo$Gauge_start == "UTLC",]$SlopeMm, 
                                   ManningsCoef = 0.037/60,
                                   seglg_m = segInfo[segInfo$Gauge_start == "UTLC",]$SegLgM, 
                                   Sms_gPgDMmin.Med = Sms_gPgDMmin.Med_i, 
                                   Scap_gPgDM = Scap_gPgDM_i, 
                                   StreamEPC = StreamEPC_i, 
                                   TimeWin = 15)%>% 
                      mutate(bootN = i)
    streamWQ_ibs <- rbind(streamWQ_ibs,streamWQ_i)
  }

  #####################
  # Little functions for calculating probs
  #####################
  p <- c(0.025, 0.5, 0.975)
  
  p_names <- map_chr(p, ~paste0(.x*100, "per"))
  
  p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = p_names)
  
  p_funs
  
  #####################
  # only make calcuate for intervals with Q over target % flow
  #####################
  
  streamWQatMRDEFbs.Sdoy <- streamWQ_ibs %>% 
    mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
           Svol_gPwindow.hf = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0)) %>% 
    select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Svol_gPwindow.hf, Qm3m) %>% 
    group_by(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Qm3m) %>% 
    summarize_at(vars(Svol_gPwindow.hf), funs(!!!p_funs)) %>% 
    mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC"),
           DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))) %>% 
    group_by(DOY) %>% 
    summarize_at(vars(DRPloadUS.gPwin:'97.5per'), list(".s" = sum)) %>% 
    select(DOY, "DRPloadUS.gPwin" = "DRPloadUS.gPwin_.s","DRPloadUS.gPwin.hf" = "DRPloadUS.gPwin.hf_.s", "Svol_gPwindow.hf_2.5per" = "2.5per_.s", 
           "Svol_gPwindow.hf_50per" = "50per_.s", "Svol_gPwindow.hf_97.5per" = "97.5per_.s") %>% 
    mutate(DRPdif_2.5per = Svol_gPwindow.hf_2.5per/DRPloadUS.gPwin.hf,
           DRPdif_50per = Svol_gPwindow.hf_50per/DRPloadUS.gPwin.hf,
           DRPdif_97.5per = Svol_gPwindow.hf_97.5per/DRPloadUS.gPwin.hf) %>% 
    mutate(DRPdif_50per = ifelse(DRPdif_50per == 0,1,DRPdif_50per),
           DRPdif_2.5per = ifelse(DRPdif_2.5per == 0,1,DRPdif_2.5per),
           DRPdif_97.5per = ifelse(DRPdif_97.5per == 0,1, DRPdif_97.5per))
  

  
  ggplot(streamWQatMRDEFbs.Sdoy) + 
    geom_point(aes(y = DRPdif_50per, x = DOY)) +
    geom_line(aes(y = DRPdif_50per, x = DOY), color = "grey") +
    geom_ribbon(aes(ymin = DRPdif_2.5per, ymax = DRPdif_97.5per, x = DOY), fill = "green", alpha = 0.25) +
    geom_line(aes(y = DRPdif_2.5per, x = DOY), color = "grey") +
    geom_line(aes(y = DRPdif_97.5per, x = DOY), color = "grey") +
    ylab("Difference in DRP load at > 75% flows")
  
  # Portion sorbed during stumpf window - high flows
  StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
  EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")
  streamWQatMRDEF.stfHF <- streamWQ_ibs %>% 
    # subset to stumpf window
    filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
    mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
           Svol_gPwindow = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0)) %>% 
    select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Svol_gPwindow, Qm3m, bootN) %>% 
    group_by(bootN) %>% 
    summarise_at(vars(DRPloadUS.gPwin, DRPloadUS.gPwin.hf, Svol_gPwindow), list(sum = sum), na.rm = TRUE) %>% 
    mutate(DRPsorbed.hf = Svol_gPwindow_sum/DRPloadUS.gPwin.hf_sum,
           DRPsorbed.allF = Svol_gPwindow_sum/DRPloadUS.gPwin_sum) 
  
  ggplot(streamWQatMRDEF.stfHF) +
    geom_density(aes(x = DRPsorbed.hf))+
    geom_rug(aes(x = DRPsorbed.hf, y = 0), position = position_jitter(height = 0)) +
    xlab("Annual difference in DRP load at > 75% flows")
  
  
  # export bootstrapped data
  write.csv(streamWQ_ibs, "04_generatedData/05d_UTLCboot.csv")
  # save image
  save.image("03_Rdata/05z_scaleUTLC_rdat")
  