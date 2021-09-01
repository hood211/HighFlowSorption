library(tidyverse)
library(MuMIn)
library(ggpubr)
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
  select(Stream, SampleDate, Percentile_flow, percentPP, percentOM, Sms_gPgDMmin:Q_m3min) %>% 
  mutate(Percentile_flow26th = Percentile_flow^6,
         Svol_gPm3Min.log = log(Svol_gPm3Min),
         SS_gDMm3.log = log(SS_gDMm3))

# range of sorption - mgP/m3/h
summary(sorp$Svol_gPm3Min)*1000*60


# SORPTION ~ Q + STREAM
# what is correlation between SS and Q
cor(sorp$SS_gDMm3.log, sorp$Percentile_flow26th, method = "pearson") #0.64
SSQbS.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th * Stream, sorp)
SSQpS.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th + Stream, sorp)
SSQ.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th, sorp)
model.sel(SSQbS.lm, SSQpS.lm, SSQ.lm)
summary(SSQpS.lm)

ggplot(sorp, aes(y = SS_gDMm3.log, x = Percentile_flow26th, color = Stream)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)

f1 <- ggplot(sorp, aes(y = SS_gDMm3.log, x = Q_m3min, color = Stream)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_color_manual(name="Stream", 
                                       values = c("LFR"="light blue", 
                                                  "PC"="light pink", 
                                                  "PR"="indianred1", 
                                                  "STF"="light salmon", 
                                                  "UTLC"="medium sea green",
                                                  "WC" = "plum2")) +
  theme_classic() + theme(legend.position = "none")  + 
  theme(axis.title.y = element_text(color ="black", size=15), 
        axis.title.x = element_text(color="black", size = 15), 
        axis.text.y= element_text(size=13),
        axis.text.x = element_text(size = 13)) 



# What is relationship between P sorp and SS
SSSbS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log * Stream, sorp)
SSSpS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log + Stream, sorp)
SSS.lm <- lm(Svol_gPm3Min.log  ~ SS_gDMm3.log, sorp)
model.sel(SSSbS.lm, SSSpS.lm, SSS.lm)

summary(SSSpS.lm)


# WHITNEY MAKE THIS LOOK LIKE A BETTER FIGURE
## WMK edited 3/31/21
ggplot(sorp, aes(y = Svol_gPm3Min.log, x = SS_gDMm3.log, color = Stream)) + 
  geom_point(size=2) +
  scale_color_manual(name="Stream", 
                                    values = c("LFR"="light blue", 
                                               "PC"="light pink", 
                                               "PR"="indianred1", 
                                               "STF"="light salmon", 
                                               "UTLC"="medium sea green",
                                               "WC" = "plum2")) +
  theme_classic() + theme(legend.position = "none")  + 
  theme(axis.title.y = element_text(color ="black", size=15), 
        axis.title.x = element_text(color="black", size = 15), 
        axis.text.y= element_text(size=13),
        axis.text.x = element_text(size = 13)) + ylab("log (Volumetric P Sorption (g P "~m^3*" Min))") +
  xlab("log(Suspended Sediment (g DM "~m^3*"))") + 
  theme(legend.position = c(0.9,0.2))

ggsave("SorpVol_SS.tif", plot = last_plot(), device = "tiff",
       width = 5, height = 6, units = "in", dpi = 300)                                        

# What effect would a decline in DRP have on sorption desorption
qqnorm(sorp$AmbP_gPm3); qqline(sorp$AmbP_gPm3)
qqnorm(sorp$EPC_gPm3); qqline(sorp$EPC_gPm3)

DRP_pQ_EPCpS.lm <- lm(AmbP_gPm3 ~ EPC_gPm3 + Percentile_flow^6 +Stream, sorp)
DRP_bQ_EPC.lm <- lm(AmbP_gPm3 ~ EPC_gPm3 * Percentile_flow^6, sorp)
DRP_pQ_EPC.lm <- lm(AmbP_gPm3 ~ EPC_gPm3 + Percentile_flow^6, sorp)
DRP_EPCbS.lm <- lm(AmbP_gPm3 ~ EPC_gPm3 * Stream, sorp)
DRP_EPCpS.lm <- lm(AmbP_gPm3 ~ EPC_gPm3 + Stream, sorp)
DRP_EPC.lm <- lm(AmbP_gPm3 ~ EPC_gPm3, sorp)
model.sel(DRP_pQ_EPCpS.lm, DRP_bQ_EPC.lm, DRP_pQ_EPC.lm, DRP_EPCbS.lm, DRP_EPCpS.lm, DRP_EPC.lm)
summary(DRP_Q_EPCbS.lm)
summary(DRP_EPC.lm)


EPC_pQ_DRPpS.lm <- lm(EPC_gPm3 ~ AmbP_gPm3 + Percentile_flow^6 +Stream, sorp)
EPC_bQ_DRP.lm <- lm(EPC_gPm3 ~ AmbP_gPm3 * Percentile_flow^6, sorp)
EPC_pQ_DRP.lm <- lm(EPC_gPm3 ~ AmbP_gPm3 + Percentile_flow^6, sorp)
EPC_DRPbS.lm <- lm(EPC_gPm3 ~ AmbP_gPm3 * Stream, sorp)
EPC_DRPpS.lm <- lm(EPC_gPm3 ~ AmbP_gPm3 + Stream, sorp)
EPC_DRP.lm <- lm(EPC_gPm3 ~ AmbP_gPm3, sorp)
model.sel(EPC_pQ_DRPpS.lm, EPC_bQ_DRP.lm, EPC_pQ_DRP.lm, EPC_DRPbS.lm, EPC_DRPpS.lm, EPC_DRP.lm)
summary(EPC_DRP.lm)

summary(sorp %>% 
          mutate(DRPmEPC = AmbP_gPm3/EPC_gPm3))


f2 <- ggplot() +
  geom_point(data = sorp, aes(x = EPC_gPm3, y = AmbP_gPm3, color = Stream)) +
  stat_smooth(data = sorp, aes(x = EPC_gPm3, y = AmbP_gPm3), method = "lm") +
  xlim(0,0.06) +
  ylim(0,0.6) +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_abline(intercept = 0, slope = 5.4*.5, color = "grey") +
  geom_segment()

ggplot(sorp, aes(y = AmbP_gPm3, x = Percentile_flow^6, size = EPC_gPm3)) +
  geom_point()

summary(lm(AmbP_gPm3 ~ Percentile_flow^6 + EPC_gPm3, sorp))

plot(residuals(SSSpS.lm) ~ sorp$Percentile_flow26th)

write.csv(sorp,"04_generatedData/04d_SorpDat.csv")

png("05_Figures/04_TribSorbFigs.png", units = "in", height = 12, width = 6, res = 300)
ggarrange(f1,f2, nrow = 2, ncol =1)
dev.off()

         
         
