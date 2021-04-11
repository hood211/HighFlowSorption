library(tidyverse)
library(MuMIn)
library(ggpubr)
library(grid)
library(egg)
library(smatr)

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


#############################################################################################
# SS load ~ Q + STREAM
#############################################################################################
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
  scale_x_log10() 



#############################################################################################
# P sorption
#############################################################################################

# What is relationship between P sorp and SS
SSSbS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log * Stream, sorp)
SSSpS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log + Stream, sorp)
SSS.lm <- lm(Svol_gPm3Min.log  ~ SS_gDMm3.log, sorp)
SSQbS.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th * Stream, sorp)
SSQpS.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th + Stream, sorp)
SSQ.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th, sorp)

model.sel(SSSbS.lm, SSSpS.lm, SSS.lm, SSQbS.lm, SSQpS.lm, SSQ.lm)

summary(SSSpS.lm)




sorp2 <- sorp %>% 
  mutate(Svol_gPm3Min.log.pred = predict(SSSpS.lm))


# Fig 2
## WMK edited 3/31/21 - JMH 4 apr 21
Fig2 <- ggplot() + 
  geom_point(data = sorp2, aes(y = Svol_gPm3Min.log, x = SS_gDMm3.log, color = Stream), size=2) +
  geom_line(data = sorp2, aes(y = Svol_gPm3Min.log.pred, x = SS_gDMm3.log, color = Stream), size=1) +
  scale_color_manual(name="Stream", 
                     values = c("LFR"="light blue", 
                                "PC"="light pink", 
                                "PR"="indianred1", 
                                "STF"="light salmon", 
                                "UTLC"="medium sea green",
                                "WC" = "plum2")) +
  theme_bw() +
  theme(legend.position = c(0.9,0.25),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank()) + 
  ylab(expression(paste("log P sorption (g P ", m^-3," ", min^-1,")" ))) +
  xlab(expression(paste("log Suspended sediment (g dry mass ",m^-3,")"))) +
  annotate("text", x = 5.5, y = -16.5, label = "R^2 == 0.79", parse = T, size = 5)+
  annotate("text", x = 5.5, y = -17.5, label = "P < 0.001", size = 5)
  


png("05_Figures/Fig2.png", units = "in", height = 5, width = 6, res = 300)
Fig2
dev.off()


#############################################################################################
# EPC
#############################################################################################
qqnorm(sorp$AmbP_gPm3); qqline(sorp$AmbP_gPm3)
qqnorm(sorp$EPC_gPm3); qqline(sorp$EPC_gPm3)


# Predicting EPC with DRP - SMA model
EPC_DRPbS.sma <- sma(EPC_gPm3 ~ AmbP_gPm3 * Stream, sorp)
summary(EPC_DRPbS.sma) #slopes are equal

EPC_DRPpS.sma <- sma(EPC_gPm3 ~ AmbP_gPm3 + Stream, sorp)
summary(EPC_DRPpS.sma) #no difference in elevation

EPC_DRP.sma <- sma(EPC_gPm3 ~ AmbP_gPm3, sorp)
summary(EPC_DRP.sma)

sorp3 <- sorp2 %>% 
  mutate(EPCpred = predict(EPC_DRP.lm))

#Fig S3
FigS3a <- ggplot() +
  geom_point(data = sorp, aes(y = EPC_gPm3, x = AmbP_gPm3, color = Stream)) +
  geom_line(data = sorp3, aes(y = EPCpred, x = AmbP_gPm3), color = "black") +
  scale_color_manual(name="Stream", 
                     values = c("LFR"="light blue", 
                                "PC"="light pink", 
                                "PR"="indianred1", 
                                "STF"="light salmon", 
                                "UTLC"="medium sea green",
                                "WC" = "plum2")) +
  theme_bw() +
  theme(legend.position = c(0.9,0.25),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank()) + 
  ylab(expression(atop("Equilibrium P concentration",paste("(",EPC[0],"; mg P ",L^-1,")")))) +
  xlab(expression(atop("Dissolved reactive P",paste( "(DRP; mg P ", L^-1,")")))) +
  annotate("text", x = 0.28, y = 0.001, label = "R^2 == 0.48", parse = T, size = 5, hjust = 0.5) +
  annotate("text", x = 0.4, y = 0.0008, label = "P < 0.001", size = 5, hjust = 0.5) +
  annotate("text", x= 0.30, y = 0.005, label = "EPC = 0.008 + 0.09 x DRP", hjust = 0.5, size = 5)

FigS3b <- ggplot() +
  geom_point(data = sorp, aes(y = EPC_gPm3, x = AmbP_gPm3, color = Stream)) +
  geom_line(data = sorp3, aes(y = EPCpred, x = AmbP_gPm3), color = "black") +
  xlim(0,0.600) +
  ylim(0,0.600) +
  scale_color_manual(name="Stream", 
                     values = c("LFR"="light blue", 
                                "PC"="light pink", 
                                "PR"="indianred1", 
                                "STF"="light salmon", 
                                "UTLC"="medium sea green",
                                "WC" = "plum2")) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank()) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  ylab(expression(atop("Equilibrium P concentration",paste("(",EPC[0],"; mg P ",L^-1,")")))) +
  xlab(expression(atop("Dissolved reactive P",paste( "(DRP; mg P ", L^-1,")")))) +
  annotate("text", x = 0.3, y = 0.34, label = "Sediments are P source", hjust = 0, size = 5, angle = 41) +
  annotate("text", x = 0.325, y = 0.28, label = "Sediments are P sink", hjust = 0, size = 5, angle = 40)


FigS3a.g <- ggplotGrob(FigS3a)
FigS3b.g <- ggplotGrob(FigS3b)

FigS3a.gtf <- gtable_frame(FigS3a.g, width = unit(0.9, "null"), height = unit(0.9, "null"))
FigS3b.gtf <- gtable_frame(FigS3b.g, width = unit(0.9, "null"), height = unit(0.9, "null"))

FigS3.gtf <- gtable_frame(gtable_cbind(FigS3a.gtf, FigS3b.gtf),
                          width = unit(18, "null"), height = unit(0.9, "null"))


png("05_Figures/FigS3.png", units = "in", height = 5, width = 12, res = 300)
grid.newpage()
grid.draw(FigS3.gtf)
grid.text("a", x = unit(0.02,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.52,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()


#####################
# SAVE/LOAD
#####################

# save.image("03_Rdata/04_SorptionData_rdat")
# load("03_Rdata/04_SorptionData_rdat")

         
         
