# JMH Jan 2022
# Cleans up sorption measurement data
# Generates Fig. 2 and Fig. S3 
# Addresses some questions about the relationship between P sorption & SS
# Examines variation in sorption at median SS
# Statistics for EPC

library(tidyverse)
library(MuMIn)
# library(ggpubr)
library(grid)
library(egg)
library(smatr)


# load data ----
sorp <-  read.csv("01_RawData/00_WK_IsoAmb_JMH.csv", row.names = 1) %>% 
  # remove outliers - excluded from other analyses
  # JAN AND FEB ALREADY EXCLUDED
  filter(sorp_capacity != 1544.88) %>%
  filter(NAP != -810.38) %>%
  filter(NAP != 2680.00) %>%
  filter(NAP != 1401.03) %>% 
  filter(slope_ambP > 0) %>% 
  mutate(Stream = as.factor(Stream),
         Stream = fct_recode(Stream, PC = "Platter",
                                     STF = "S_Turkey",
                                     WC = "West",
                                     UTLC = "Lost_Trib",
                                     LFR = "LFR",
                                     PR = "Potato_Run")) %>% 
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
  mutate(Percentile_flow26th = Percentile_flow^6,
         Svol_gPm3Min.log = log(Svol_gPm3Min),
         SS_gDMm3.log = log(SS_gDMm3)) 

# range of sorption - mgP/m3/h
summary(sorp$Svol_gPm3Min)*1000*60


# Predicting SS concentration ----
    
    # what is correlation between SS and Q
    cor(sorp$SS_gDMm3.log, sorp$Percentile_flow26th, method = "pearson") #0.70
    
    # best model of SS concentration
    SSQbS.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th * Stream, sorp)
    SSQpS.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th + Stream, sorp)
    SSQ.lm <- lm(SS_gDMm3.log ~ Percentile_flow26th, sorp)
    model.sel(SSQbS.lm, SSQpS.lm, SSQ.lm)
    summary(SSQbS.lm)
    
    # quick plots
    ggplot(sorp, aes(y = SS_gDMm3.log, x = Percentile_flow26th, color = Stream)) +
      geom_point() +
      stat_smooth(method = "lm", se = FALSE)
    
    ggplot(sorp, aes(y = SS_gDMm3.log, x = Q_m3min, color = Stream)) +
      geom_point() +
      stat_smooth(method = "lm", se = FALSE) +
      scale_x_log10() 

# Predicting P sorption ----
    # Model predicting P srob
    PSSSbS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log * Stream, sorp)
    PSSSpS.lm <- lm(Svol_gPm3Min.log ~ SS_gDMm3.log + Stream, sorp)
    PSSS.lm <- lm(Svol_gPm3Min.log  ~ SS_gDMm3.log, sorp)
    PSSQbS.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th * Stream, sorp)
    PSSQpS.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th + Stream, sorp)
    PSSQ.lm <- lm(Svol_gPm3Min.log ~ Percentile_flow26th, sorp)
    
    model.sel(PSSSbS.lm, PSSSpS.lm, PSSS.lm, PSSQbS.lm, PSSQpS.lm, PSSQ.lm)
    
    summary(PSSSbS.lm)
    
    sorp$Svol_gPm3Min_pred <- exp(predict(PSSSbS.lm, sorp))
    plot(sorp$Svol_gPm3Min_pred ~ sorp$Svol_gPm3Min)

    # Calcs supporting Results ----
# At the median sediment concentration across all streams (233.89 g dry mass m-3), 
# P sorption was ?-fold higher between the stream with the lowest (???) and highest (???) P sorption.

    SSmedian <- median(sorp$SS_gDMm3)
    sorpAtMedianSS <- as.data.frame(cbind(SS_gDMm3.log = log(rep(SSmedian, each = 6)), 
                                          Stream = unique(as.character(sorp$Stream))))
    sorpAtMedianSS$SS_gDMm3.log = as.numeric(sorpAtMedianSS$SS_gDMm3.log)
    # predicted sorption in mg P/m3/h
    sorpAtMedianSS$Svol_gPm3Min_pred <- exp(predict(PSSSbS.lm, sorpAtMedianSS))  *60*1000

    # percent difference b/w min and max sorption at median SS
    max(sorpAtMedianSS$Svol_gPm3Min_pred)/min(sorpAtMedianSS$Svol_gPm3Min_pred)-1

# Add to dataframe
    sorp2 <- sorp %>% 
      mutate(Svol_gPm3Min.log.pred = predict(PSSSbS.lm),
             Svol_mgPm3hr.log.pred = log(exp(Svol_gPm3Min.log.pred)*1000*60),
             Svol_mgPm3hr.log = log(exp(Svol_gPm3Min.log)*1000*60))


# Predicting EPC ----
# check normality
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
  # coefs from EPC_DRP.sma
  mutate(EPCpred = 0.002926867 + 0.1251420*(AmbP_gPm3))

# Fig S4 ----
FigS4 <- ggplot() +
  geom_point(data = sorp, aes(y = EPC_gPm3, x = AmbP_gPm3*0.5, color = Stream)) +
  geom_line(data = sorp3, aes(y = EPCpred, x = AmbP_gPm3*0.5), color = "black") +
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
  theme(legend.position = c(0.125,0.75),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank()) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  ylab(expression(atop("Equilibrium P concentration",paste("(",EPC[0],"; mg P ",L^-1,")")))) +
  xlab(expression(atop("Dissolved reactive P",paste( "(DRP; mg P ", L^-1,")")))) +
  annotate("text", x = 0.3, y = 0.34, label = "Sediment desorbs P", hjust = 0, size = 5, angle = 41) +
  annotate("text", x = 0.325, y = 0.28, label = "Sediment sorbs P", hjust = 0, size = 5, angle = 40) +
  annotate("text", x = 0.41, y = 0.18, label = "R^2 == 0.56", parse = T, size = 5, hjust = 0.5) +
  annotate("text", x = 0.53, y = 0.175, label = "; P < 0.001", size = 5, hjust = 0.5) +
  annotate("text", x= 0.46, y = 0.13, label = "EPC = 0.003 + 0.13 x DRP", hjust = 0.5, size = 5)

png("05_Figures/01_FigS4AgTalk.png", units = "in", height = 5, width = 6, res = 300)
FigS4
dev.off()

# Fig S5 ----
## JMH 27 Feb 2022
FigS5 <- ggplot() + 
  geom_point(data = sorp2, aes(y = Svol_mgPm3hr.log, x = SS_gDMm3.log, color = Stream), size=2) +
  geom_line(data = sorp2, aes(y = Svol_mgPm3hr.log.pred, x = SS_gDMm3.log, color = Stream), size=1) +
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
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank()) + 
  ylab(expression(paste("log P sorption (mg P ", m^-3," ", h^-1,")" ))) +
  xlab(expression(paste("log Suspended sediment (g dry mass ",m^-3,")"))) +
  annotate("text", x = 5.5, y = -3.0, label = "R^2 == 0.90", parse = T, size = 5)+
  annotate("text", x = 5.5, y = -4.0, label = "P < 0.001", size = 5)



png("05_Figures/01_FigS5.png", units = "in", height = 5, width = 6, res = 300)
FigS5
dev.off()


# EXPORT DATA ----
# write.csv(sorp2 %>%
#             # only greater than 75 percentile flow
#             filter(Percentile_flow >= 75), "04_generatedData/01_CleanedSorpIsoDat.csv")

# SAVE/LOAD ----
# save.image("03_Rdata/01_SorpDataMungeFig_rdat")
# load("03_Rdata/01_SorpDataMungeFig_rdat")

         
         
