

TTfunJob <- function(Q, WA, Q_annual, ReachLg, slope){
  # INPUTS
  # Q = discharge at time of measurements -m3/s
  # Q = 10000
  # # WA = drainage area - m2
  # WA = 12240000
  # # Qannual = mean annual discharge - m3/s
  # Q_annual = 0.106
  # # ga = acceleration of gravity - m/s2
  # # Reach length - m
  # ReachLg = 100
  # # slope - m/m
  # slope = 0.126
  
  # JOBSON etal 1997 eq 9  
  # v_pd = dimentionless peak velocity
  # v_p = !!!!!!!!!!
  # v_pd <- v_p * WA * 1/Q 

  # JOBSON et al. 1997 eq 10  
  # WA_d = dimmentionless peak velocity
  ga <-  9.80665 #m/s2
  WA_d <- (WA^1.25 * ga^0.5 )/Q_annual
 
  # JOBSON et al. 1997 eq 11 
  # dimensionless relative discharge
  Qannual_d <- Q/Q_annual
  
  # velocity at concentration peak - m/s
  V_p <-  0.094 + 0.0143 * (WA_d^-0.919) * (Qannual_d^-0.469) * (slope^0.159) * Q/WA # 
  
  # travel time - s
  tt <- ReachLg/V_p
  
  tt
}

TTfunMor <- function(Q, WA, Q_annual, ReachLg, slope){
  # Morales et al. 2007, Journal of hydrology
  # WA_d = dimmentionless peak velocity
  ga <-  9.80665 #m/s2
  WA_d <- (WA^1.25 * ga^0.5 )/Q_annual
  
  # JOBSON et al. 1997 eq 11 
  # dimensionless relative discharge
  Qannual_d <- Q/Q_annual
  
  # velocity at concentration peak - m/s
  # V_p <-  6171.3 * (WA_d^0.2667) * (Qannual_d^-0.5956) * (slope^0.4548) * Q/WA
  V_p <-  2900207 * (Qannual_d^-0.5956) * (slope^0.5149) * Q/WA
  # travel time - s
  tt <- ReachLg/V_p
  
  tt
}

TTfunDu <- function(Q,slope,n, ReachLg){
  V_p <- slope^(3/8) * Q^(1/4) * n^(-3/4)
  tt <- ReachLg/V_p
}

# Function for scaling up sorption rate using the median sorption rate
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
# needs to have 
# streamWQ2 <- streamWQ %>% 
#   mutate(dateTime = as.POSIXct(dateTime, format = "%Y-%m-%d"),
#          Yf = as.factor(as.character(strftime(dateTime, format = "%Y"))),
#          Y = as.numeric(as.character(Yf)),
#          DOY = as.numeric(as.character(strftime(dateTime, format = "%j"))), 
#          gTarPerFlow = as.factor(gTarPerFlow)) %>% 
#   filter(Y == 2019) %>% 
#   mutate(Qm3m = Qm3s_utlc*60, #Q in m3/min #NOTE THAT THESE NEED TO BE FOR THE FOCAL STREAM
#          SSgDMm3 = SS_mgLp*1000/1000,
#          DRPgPm3 = SRP_mgLp*1000/1000)

ScalingFun.step1 <- function(WQdat, Qm3m, SSgDMm3, DRPgPm3, slope_mm, ManningsCoef,
                             seglg_m, Sms_gPgDMmin.Med, Scap_gPgDM, StreamEPC, TimeWin){
  WQdat2 <- WQdat %>% 
    mutate(tt_du = TTfunDu(Q = Qm3m, slope_mm, ManningsCoef, seglg_m)) %>%  # units = minutes
    # potential ms sorption rate gP/gDM for the entire tt - needs MINS
    mutate(
      # potential ms sorption rate gP/gDM for the entire tt - needs MINS
          Sms_gPgDMtt.P = Sms_gPgDMmin.Med * tt_du, #using tt_du here seems most reliable
           # final ms sorption rate gP/gDM for entire tt
           Sms_gPgDMtt.F = ifelse(Sms_gPgDMtt.P > Scap_gPgDM, Scap_gPgDM, Sms_gPgDMtt.P),
           # what percent of capacity?
           Sms_porCap = Sms_gPgDMtt.P/Scap_gPgDM,
           # potential volumetric sorption rate gP/m3/tt
           Svol_gPm3TT.P = Sms_gPgDMtt.F * SSgDMm3,
           # difference between DRP and EPC gP/m3
           DifDRP.EPC = ifelse((DRPgPm3 - StreamEPC)>0, DRPgPm3 - StreamEPC,0),
           # final volumetric sorption rate gP/m3/tt
           Svol_gPm3TT.F = ifelse(Svol_gPm3TT.P < DifDRP.EPC, Svol_gPm3TT.P, DifDRP.EPC),
           # potential volumetric sorption rate (w/o constraints on capacity or EPC)
           Svol_gPm3TT.P2 = Sms_gPgDMtt.P * SSgDMm3,
           # how far did the particle travel while sorbing
           SorbDistM = Svol_gPm3TT.F/Svol_gPm3TT.P2 * seglg_m) %>% 
    # NOW THE SCALING
    # DRP moving DS during between samples gP/window
    mutate(DRPloadUS.gPwin = DRPgPm3 * Qm3m * TimeWin,
           # SS moving DS during window g DM/window
           SSloadUS.gDMwin = SSgDMm3 * Qm3m * TimeWin,
           # Vol sorption for all water in window gDM/window
           Svol_gPwindow = Svol_gPm3TT.F * Qm3m * TimeWin,
           # DRP load at DS site without sorption; gP/window
           DRPloadDS.gPwin.wS = DRPloadUS.gPwin + Svol_gPwindow,
           # portion of DRP sorbed
           porDRPloadSorbed = Svol_gPwindow/DRPloadUS.gPwin,
           gPsorb_gSS = Svol_gPwindow/SSloadUS.gDMwin,
           Sms_gPgDMmin.Med = Sms_gPgDMmin.Med,
           Scap_gPgDM = Scap_gPgDM,
           StreamEPC = StreamEPC)
  
  WQdat2
}


ScalingFun.2leg <- function(WQdat.l1, WQdat.l2, Qm3m, SSgDMm3, DRPgPm3, slope_mm.l1, slope_mm.l2, ManningsCoef,
                             seglg_m.l1, seglg_m.l2, Sms_gPgDMmin.Med, Scap_gPgDM, StreamEPC, TimeWin){
  WQdatFin <- WQdat.l1  %>% 
      rename(Qm3m.l1 = Qm3m) %>% 
      mutate(tt_du.l1 = TTfunDu(Q = Qm3m.l1, slope_mm.l1, ManningsCoef, seglg_m.l1)) %>%  # units = minutes
      mutate(dateTimeAtl2 = round_date(dateTime + (tt_du.l1*60), unit = "15 minutes")) %>% 
      left_join(WQdat.l2 %>% 
                  select(dateTime, Qm3m2.l2 = Qm3s_Turk), by = c("dateTimeAtl2" = "dateTime")) %>% 
      mutate(tt_du.l2 = TTfunDu(Q = Qm3m2.l2, slope_mm.l2, ManningsCoef, seglg_m.l2)) %>% 
      # LEG NUMBER 1

      mutate(
             # potential ms sorption rate gP/gDM for the entire tt - needs MINS
             Sms_gPgDMtt.l1.P = Sms_gPgDMmin.Med * tt_du.l1, #using tt_du here seems most reliable
             # final ms sorption rate gP/gDM for entire tt
             Sms_gPgDMtt.l1.F = ifelse(Sms_gPgDMtt.l1.P > Scap_gPgDM, Scap_gPgDM, Sms_gPgDMtt.l1.P),
             # what percent of capacity?
             Sms_porCap.l1 = Sms_gPgDMtt.l1.P/Scap_gPgDM,
             # potential volumetric sorption rate gP/m3/tt
             Svol_gPm3TT.l1.P = Sms_gPgDMtt.l1.F * SSgDMm3,
             # difference between DRP and EPC gP/m3
             DifDRP.EPC.l1 = ifelse((DRPgPm3 - StreamEPC)>0, DRPgPm3 - StreamEPC,0),
             # final volumetric sorption rate gP/m3/tt
             Svol_gPm3TT.l1.F = ifelse(Svol_gPm3TT.l1.P < DifDRP.EPC.l1, Svol_gPm3TT.l1.P, DifDRP.EPC.l1),
             # potential volumetric sorption rate (w/o constraints on capacity or EPC)
             Svol_gPm3TT.l1.P2 = Sms_gPgDMtt.l1.P * SSgDMm3,
             # how far did the particle travel while sorbing
             SorbDistM.l1 = Svol_gPm3TT.l1.F/Svol_gPm3TT.l1.P2 * seglg_m.l1) %>% 
      # LEG NUMBER 2
      # potential ms sorption rate gP/gDM for the entire tt - needs MINS
      mutate(Sms_gPgDMtt.l2.P = Sms_gPgDMmin.Med * tt_du.l2, #using tt_du here seems most reliable
             # how much of the capacity is remaining for this leg?
             Scap_gPgDMremaning.l2 = ifelse(Sms_porCap.l1 >1, Scap_gPgDM*0, Scap_gPgDM*(1-Sms_porCap.l1)),
             # final ms sorption rate gP/gDM for entire tt
             Sms_gPgDMtt.l2.F = ifelse(Sms_gPgDMtt.l2.P > Scap_gPgDMremaning.l2, Scap_gPgDMremaning.l2, Sms_gPgDMtt.l2.P),
             # what percent of capacity?
             Sms_porCap.l2 = Sms_gPgDMtt.l2.F/Scap_gPgDMremaning.l2,
             # potential volumetric sorption rate gP/m3/tt
             Svol_gPm3TT.l2.P = Sms_gPgDMtt.l2.F * SSgDMm3, #Just for sediments from WC
             # difference between DRP and EPC gP/m3
             DifDRP.EPC.l2 = ifelse(((DRPgPm3- Svol_gPm3TT.l1.F) - StreamEPC)>0, (DRPgPm3- Svol_gPm3TT.l1.F) - StreamEPC, 0),
             # final volumetric sorption rate gP/m3/tt
             Svol_gPm3TT.l2.F = ifelse(Svol_gPm3TT.l2.P < DifDRP.EPC.l2, Svol_gPm3TT.l2.P, DifDRP.EPC.l2),
             # potential volumetric sorption rate (w/o constraints on capacity or EPC)
             Svol_gPm3TT.l2.P2 = Sms_gPgDMtt.l2.P * SSgDMm3,
             # how far did the particle travel while sorbing
             SorbDistM.l2 = Svol_gPm3TT.l2.F/Svol_gPm3TT.l2.P2 * seglg_m.l2) %>% 
      # Total volumetric sorption
      mutate(Svol_gPm3TT.F = Svol_gPm3TT.l1.F + Svol_gPm3TT.l2.F,
             # DRP moving DS during between samples gP/window
             DRPloadUS.gPwin = DRPgPm3 * Qm3m.l1 * TimeWin,
             # SS moving DS during window g DM/window
             SSloadUS.gDMwin = SSgDMm3 * Qm3m.l1 * TimeWin,
             # Vol sorption for all water in window gDM/window
             Svol_gPwindow = Svol_gPm3TT.F * Qm3m.l1 * TimeWin,
             # DRP load at DS site without sorption; gP/window
             DRPloadDS.gPwin.wS = DRPloadUS.gPwin + Svol_gPwindow,
             # portion of DRP sorbed
             porDRPloadSorbed = Svol_gPwindow/DRPloadUS.gPwin,
             # sorption rate per g SS
             gPsorb_gSS = Svol_gPwindow/SSloadUS.gDMwin,
             SorbDistM = SorbDistM.l1 + SorbDistM.l2,
             Sms_gPgDMmin.Med = Sms_gPgDMmin.Med,
             Scap_gPgDM = Scap_gPgDM,
             StreamEPC = StreamEPC)
  
  WQdatFin
}

MaumFun <- function(WQDat, gPsorbedgDMe){
  WQDat2 <- WQDat %>% 
    select(dateTime, SmpTimeWindowDay, SS_mgL, TP_mgL, SRP_mgL, Qm3s) %>% 
    # estimate P sorbed
    mutate(gPsorbedgDM = gPsorbedgDMe,
           # flow in window
           Qm3window = (Qm3s*60*60*24) * SmpTimeWindowDay, #convert Q to 1/day - units gP/window (window varies)
           # estimate DRP load b/w samples - observed
           gDRPwindow = SRP_mgL * Qm3window,
           # estimate SS load b/w samples - observed
           gSSwindow = SS_mgL * Qm3window,
           gTPwindow = TP_mgL * Qm3window,
           # estiamte DRP load b/w - sorbed
           # this is potential, which can exceed DRP if not bounded
           gPsorbWindowP = gPsorbedgDM * SS_mgL *  Qm3window,
           # bound sorption to available P
           gPsorbWindow = ifelse(gPsorbWindowP > gDRPwindow, gDRPwindow, gPsorbWindowP),
           # estiamte DRP load - observed + sorbed)
           gDRPwindowWOsorp = gDRPwindow + gPsorbWindow)
  
  WQDat2
}


#####################
# Little fuctions for calculating probs
#####################
p <- c(0.025, 0.5, 0.975)

p_names <- map_chr(p, ~paste0(.x*100, "per"))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

p_funs