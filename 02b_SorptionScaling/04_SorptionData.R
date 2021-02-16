library(tidyverse)

#############################################################################################
# load data
#############################################################################################
sorp <- read.csv("01_RawData/RawSorptionData.csv", row.names = 1)[,-1] %>% 
  # remove outliers - excluded from other analyses
  filter(sorp_capacity != 1544.88) %>%
  filter(NAP != -810.38) %>%
  filter(NAP != 2680.00) %>%
  filter(NAP != 1401.03) %>% 
  filter(slope_ambP > 0) %>% 
  mutate(Stream = as.factor(Stream),
         Stream = fct_recode(Stream, PC = "Platter Creek",
                                     STF = "South Turkey Foot",
                                     WC = "West Creek",
                                     UTLC = "Unnamed Trib to Lost Creek",
                                     LFR = "Little Flat Rock",
                                     PR = "Potato Run")) %>% 
          # mass-specific sorption rate (g P/gDM/min)
  mutate(Sms_gPgDMmin = mass_sorp_P_DM_m/1000/1000/60, #convert from ug P to g P and from h to min
         # volumetric sorption rate (g P/m3/min)
         Svol_gPm3Min = Sms_gPgDMmin*(gDM_L_avg*1000), #convert gDM/L to gDM/m3
         # Native P (g P/gDM)
         NAP_gPgDM = NAP/1000/1000, # convert from ugP/gDM to gP/gDM
         # Sorption capacity (Smax - NAP; gP/gDM)
         Scap_gPgDM = sorp_capacity/1000/1000, #convert from µgP/gDM to gP/gDM
         # Smax (gP/gDM)
         Smax_gPgDM = (Scap_gPgDM + NAP_gPgDM)) %>% 
  #converting some other parameters
  mutate(AmbP_gPm3 = AmbP/1000/1000*1000, #converting from ugP/L to gPm3
         EPC_gPm3 = EPC/1000/1000*1000, #converting from ugP/L to gP/m3
         SS_gDMm3 = gDM_L_avg*1000, #converting from gDM/L to gDM/m3
         PP_gPm3 = ugP_L_avg/1000/1000*1000, # converting from ugP/L to gP/m3
         Q_m3min = Discharge /35.31467*60) %>%  #convert from CFS to m3/min
  mutate(SampleDate = as.POSIXct(SampleDate, format = "%m/%d/%y")) %>% 
  select(Stream, SampleDate, Percentile_flow, percentPP, percentOM, Sms_gPgDMmin:Q_m3min)

ggplot(sorp, aes(y = Sms_gPgDMmin, x = Q_m3min)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(method = "lm")

summary(lm(log(Sms_gPgDMmin) ~ log(Q_m3min), sorp))



write.csv(sorp,"04_generatedData/04d_SorpDat.csv")


         
         
