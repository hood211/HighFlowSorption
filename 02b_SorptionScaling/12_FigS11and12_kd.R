# JMH Feb 2020
#  Partitioning analyses (kd) 
# for both tribs and Maumee River
# generates Fig. S11 and S12

# library ----
library(tidyverse)
library(mgcv)
library(MuMIn)
library(mgcViz)
library(tidymv)
library(ggpubr)
library(rsq)

# Get data ----
# isotherm data 
isos <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "01_CleanedSorpIsoDat.csv"), row.names = 1) %>% 
            select(Date, Stream, AmbP, gDM_L_avg, ugP_L_avg, est_K, est_Qmax, est_NAP, P_change_ugP_gMODatAmb) %>% 
            # partitioning of isotherms
            # this is from Pu et al. 2021, Chemosphere, 263 128334
            mutate(Kd_iso = P_change_ugP_gMODatAmb/AmbP*1000, # convert from L/g to mL/g
                   Kd_Lin = ugP_L_avg/(AmbP * gDM_L_avg)*1000, # convert from L/g to mL/g
                   logKd_iso = log(Kd_iso),
                   logKd_Lin = log(Kd_Lin))

# Maumee at Waterville nutrient data
wqMaum <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "02d_MaumeeWatervilleWaterQual.csv"), row.names = 1) %>% 
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d"),
         HighFlow = ifelse(HighFlow == "HighFlow", "G75th", "L75th"),
         # Y = as.character(strftime(dateTime, format = "%Y")),
         Mc = as.factor(as.character(strftime(Date, format = "%m"))),
         # create month groups: Mar-June, Jul-Sep, Oct-Feb
         TargetMonths = as.factor(ifelse(Mc %in% c("03", "04", "05", "06"), "MarJun",
                                         ifelse(Mc %in% c("07", "08", "09"), "JulSep","OctFeb")))) %>% 
  # remove 1978, 79, 80,81; could leave 78 for spring stuff
  filter(Y != 1978 & Y != 1979 & Y != 1980 & Y != 1981) %>% 
          # particulate P, 0.15 is the average TDP/SRP ratio from Baker 2014
  mutate(PP_gm3 = TP_gm3 - SRP_gm3 - SRP_gm3*0.15,
         perP = PP_gm3/SS_gm3,
         # equation from Lin et al. 2016 Science of the Total Environment, section 2.5
         # units: mL/g, same as Lin et al. 2016
         Kd = PP_gm3/(SRP_gm3 * SS_gm3)*1e6,# convert from m3/g to mL/g
         logKd = log(Kd),
         logSS_gm3 = log(SS_gm3),
         logQm3day = log(Qm3day),
         Yf = as.factor(as.character(Y))) %>% 
  filter(TargetMonths == "MarJun") %>% 
  filter(HighFlow == "G75th")

# average Kd's
summary(isos$Kd_Lin)
wqMaum %>% 
  filter(Y == 2019) %>% 
  group_by(Y) %>% 
  summarise(mean = mean(Kd, na.rm = T))


# Patterns in MR Kd ----
# not using log link because conceptual model focuses on log-log relationships
# Generally done using linear model, but allowing for nonlinear change in elevation with Y
# used model 4 to evaluate several family options for each response including 
# Gamma-log, Gamma-identity, Gaussian-log, Gaussian-identity, scat-log, and scat-identity
# also tw with sqrt, log, identity, power

qqnorm(wqMaum$logKd); qqline(wqMaum$logKd)

# starting with very complex model
gam_Kd_1 <- gam(logKd ~ te(logSS_gm3, Qm3day) + te(logSS_gm3, Y) + te(logSS_gm3, DOY),
                 data = wqMaum,
                 family = scat(link = "identity"),
                 method = "ML")

gam_Kd_2 <- gam(logKd ~ te(logSS_gm3, Qm3day) + s(Y) + te(logSS_gm3, DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_3 <- gam(logKd ~ te(logSS_gm3, Qm3day) + te(logSS_gm3, Y) + s(DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_4 <- gam(logKd ~ te(logSS_gm3, Qm3day) + s(Y) + s(DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_5 <- gam(logKd ~ te(logSS_gm3, Y) + + s(DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_6 <- gam(logKd ~ te(logSS_gm3, DOY) + s(Y),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_7 <- gam(logKd ~ s(logSS_gm3) + s(Qm3day) + s(Y) + s(DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")

gam_Kd_8 <- gam(logKd ~ s(logSS_gm3) + s(Y) + s(DOY),
                data = wqMaum,
                family = scat(link = "identity"),
                method = "ML")


# going with 8 instead of 4 because it is MUCH simpler and has a similar variation explained
model.sel(gam_Kd_1, gam_Kd_2, gam_Kd_3, gam_Kd_4, gam_Kd_5, gam_Kd_6, gam_Kd_7, gam_Kd_8)

# this is not fantastic
check(getViz(gam_Kd_8))
summary(gam_Kd_8)
plot(gam_Kd_8)

# Kd's in iso's 
ggplot(isos, aes(y = logKd_iso, x = logKd_Lin, color = Stream)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

ggplot() +
  geom_point(data = isos, aes(y = logKd_iso, x = log(gDM_L_avg), color = Stream), shape = 21) +
  geom_point(data = isos, aes(y = logKd_Lin, x = log(gDM_L_avg), color = Stream), shape = 22) +
  geom_abline(intercept = 0, slope = 1)


# combine Kd datasets-
wqMaumKd <- wqMaum %>% 
            mutate(Stream = "MaumeeWV",
                   KdType = "Lin_Maumee") %>% 
            select(Date, Stream,KdType, logKd, logSS_gm3) 

IsoKd <- isos %>% 
            select(Date, Stream, gDM_L_avg, logKd_iso, logKd_Lin) %>% 
            distinct() %>% 
            pivot_longer(cols = logKd_iso:logKd_Lin, names_to = "KdType", values_to = "logKd") %>% 
            separate(KdType, into = c("blah", "KdType"), sep = "_") %>% 
            mutate(KdType = paste0(KdType, "_", Stream),
                   logSS_gm3 = log(gDM_L_avg * 1000)) %>% 
            select(Date, Stream,KdType, logKd, logSS_gm3) 

Kds <- rbind(wqMaumKd, IsoKd) %>% 
          separate(KdType, into = c("KdType", "Stream2"), sep = "_")

## make predictions ----
Fig1Pred <- predict_gam(gam_Kd_8,
                    values = list(logSS_gm3 = seq(min(wqMaum$logSS_gm3), max(wqMaum$logSS_gm3), length = 200),
                                  #median values
                                  Y = 2000,
                                  DOY = 107))

Fig2Pred <- predict_gam(gam_Kd_8,
                        values = list(Y = seq(1975,2019, by =1),
                                      #median values
                                      logSS_gm3 = 4.923169,
                                      DOY = 107)) 

Fig3Pred <- predict_gam(gam_Kd_8,
                        values = list(DOY = seq(60,182, by = 1),
                                      #median values
                                      Y = 2000,
                                      logSS_gm3 = 4.923169)) 

# Patterns in Trib Kd ----
  qqnorm(Kds$logKd); qqline(Kds$logKd)
  
  Kds_Lin <- Kds %>% 
        filter(KdType == "Lin" & Stream != "MaumeeWV")

  glm_KdTrib_1 <- glm(logKd ~ logSS_gm3 * Stream,
                      family = Gamma(link = "identity"),
                      data = Kds_Lin)
  
  glm_KdTrib_2 <- glm(logKd ~ logSS_gm3 + Stream,
                      family = Gamma(link = "identity"),
                      data = Kds_Lin)
  
  glm_KdTrib_3 <- glm(logKd ~ logSS_gm3,
                      family = Gamma(link = "identity"),
                      data = Kds_Lin)

  model.sel(glm_KdTrib_1, glm_KdTrib_2, glm_KdTrib_3)
  plot(glm_KdTrib_2)
  summary(glm_KdTrib_2)
  rsq(glm_KdTrib_2)
  
  Kds_Lin$logKd_pred <- predict(glm_KdTrib_2, Kds_Lin)
  
# Kd ~ Watershed area ----
  PredictTrib = as.data.frame(cbind(Stream = unique(Kds_Lin$Stream), 
                                    logSS_gm3 = rep(log(100), length = 6),
                                    WA = c(39.6, 50.5, 10.8, 300.4, 40.2, 44.5)))
  PredictTrib2 <- PredictTrib %>% 
    mutate(logSS_gm3 = as.numeric(logSS_gm3),
           WA = as.numeric(WA))
  
  PredictTrib2$logKd_pred <- predict(glm_KdTrib_2, PredictTrib2,se = T)$fit
  PredictTrib2$logKd_pred.se <- predict(glm_KdTrib_2, PredictTrib2,se = T)$se.fit
  
  PredictTrib2[7,] <- c("MR", 5, 17000, as.numeric(NA), as.numeric(NA))
  
  MaumKdPred2019 <- predict_gam(gam_Kd_8,
                                values = list(logSS_gm3 = log(100),
                                              #median values
                                              Y = 2019,
                                              DOY = 107))
  
  PredictTrib2[7,4] <- as.numeric(MaumKdPred2019[4])
  PredictTrib2[7,5] <- as.numeric(MaumKdPred2019[5])
  
  PredictTrib3 <- PredictTrib2 %>% 
    mutate_at(vars(logSS_gm3:logKd_pred.se), as.numeric) %>% 
    mutate(log10WA = log10(WA))
  
  summary(PredictTrib3 %>% filter(Stream != "MR"))
  
  summary(lm(logKd_pred ~ log10WA, data = PredictTrib3))
  
 
  
  
  

# plots ----
  ## Fig. S10 ----
  # Fig. S10a ----
  pS10a <-   ggplot() +
              geom_point(data = Kds_Lin, aes(y = logKd, x = logSS_gm3, color = Stream), size = 3) +
              geom_line(data = Kds_Lin, aes(y = logKd_pred, x = logSS_gm3, color = Stream))+ 
              theme_bw() +
              theme(legend.position = c(0.8,0.7),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    legend.key.height = unit(0.5,"cm"),
                    legend.spacing.y = unit(0.05,"cm"),
                    legend.margin = margin(0.1,0,0,0, unit="cm"),
                    legend.background = element_rect(fill = "transparent"),
                    axis.title = element_text(size = 18),
                    axis.text = element_text(size = 12),
                    panel.background = element_rect(fill = "transparent"),
                    panel.grid.minor = element_blank()) +
              ylab("log Kd") +
              xlab(expression(paste("log suspended sediment (g ",m^-3,")"))) +
              ylim(7,16)+
              xlim(2,8) +
              annotate("text", x = 2, y = 16, label = "a", size = 8, fontface = "bold")
            
    
  # Fig. S10b ----
 pS10b <-  ggplot() +
    geom_point(data = wqMaum, aes(y = logKd, x = logSS_gm3), shape = 21, fill = "transparent") +
    geom_ribbon(data = Fig1Pred, aes(ymin = fit - se.fit, ymax = fit + se.fit, x= logSS_gm3), fill = "pink", alpha = 0.5)  +
    geom_line(data = Fig1Pred, aes(y = fit, x= logSS_gm3))  + 
    theme_bw() +
    theme(legend.position = c(0.21,0.70),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 12),
          legend.key.height = unit(0.5,"cm"),
          legend.spacing.y = unit(0.05,"cm"),
          legend.margin = margin(0.1,0,0,0, unit="cm"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) +
    ylab("log Kd") +
    xlab(expression(paste("log suspended sediment (g ",m^-3,")")))+
   ylim(7,16)+
   xlim(2,8) +
   annotate("text", x = 2, y = 16, label = "b", size = 8, fontface = "bold")

  # Fig. S10c ----
  pS10c <- ggplot() +
    geom_ribbon(data = Fig2Pred, aes(ymin = fit - se.fit, ymax = fit + se.fit, x= Y), fill = "pink", alpha = 0.5)  +
    geom_line(data = Fig2Pred, aes(y = fit, x= Y)) +
    theme_bw() +
    theme(legend.position = c(0.21,0.70),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 12),
          legend.key.height = unit(0.5,"cm"),
          legend.spacing.y = unit(0.05,"cm"),
          legend.margin = margin(0.1,0,0,0, unit="cm"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) +
    ylab("log Kd") +
    xlab("Year") +
    ylim(9.5, 11.5) +
    annotate("text", x = 1975, y = 11.5, label = "c", size = 8, fontface = "bold")
  
  # Fig. S10d ----
  pS10d <- ggplot() +
    geom_ribbon(data = Fig3Pred, aes(ymin = fit - se.fit, ymax = fit + se.fit, x= DOY), fill = "pink", alpha = 0.5)  +
    geom_line(data = Fig3Pred, aes(y = fit, x= DOY)) +
    theme_bw() +
    theme(legend.position = c(0.21,0.70),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 12),
          legend.key.height = unit(0.5,"cm"),
          legend.spacing.y = unit(0.05,"cm"),
          legend.margin = margin(0.1,0,0,0, unit="cm"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) +
    ylab("log Kd") +
    xlab("DOY")+
    ylim(9.5, 11.5) +
    annotate("text", x = 60, y = 11.5, label = "d", size = 8, fontface = "bold")
  
  png("05_Figures/12_FigS10.png", units="in", width=9, height=7, res=300)
  ggarrange(pS10a, pS10b, pS10c, pS10d, nrow = 2, ncol = 2)
  dev.off()
  
  
## Fig S11 ----
  pS11 <- ggplot(PredictTrib3, aes(y = logKd_pred, x = WA)) +
    geom_point(size = 3, shape = 21, fill = "lightblue") +
    stat_smooth(method = "lm", color = "black", fill = "lightgrey")+
    scale_x_log10()+
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank()) +
    ylab("log Kd") +
    xlab(expression(paste("Watershed area (",km^2,")"))) +
    annotate("text", x = 1200, y = 9.5, label = "y = 9.0 + 0.4x", size = 6, hjust = 0.5) +
    annotate("text", x = 1200, y = 9.4, label = "R2 = 0.92; P < 0.001", size = 6, hjust = 0.5)
  
  png("05_Figures/12_FigS11.png", units="in", height=5, width=6,res=300)
  pS11
  dev.off()
  
# save image ----
  save.image(file.path(here::here("03_Rdata"), "12_FigS10and11_Kd_Rdat"))
  