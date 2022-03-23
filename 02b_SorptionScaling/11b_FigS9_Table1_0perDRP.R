# JMH Feb 2021
# This makes Fig S9, was Fig 4
# Also does analyses for Table 1

# libraries ----
library(tidyverse)
library(grid)
library(egg)
library(MuMIn)
library(rcompanion)

# Data ----
# original data that went into bootstrapping
# same code for 08 Maumee loads
MR <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "02d_MaumeeWatervilleWaterQual.csv"), row.names = 1) %>% 
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d"),
         HighFlow = ifelse(HighFlow == "HighFlow", "G75th", "L75th"),
         # Y = as.character(strftime(dateTime, format = "%Y")),
         Mc = as.factor(as.character(strftime(Date, format = "%m"))),
         # create month groups: Mar-June, Jul, Aug-Feb
         TargetMonths = as.factor(ifelse(Mc %in% c("03", "04", "05", "06"), "MarJun",
                                         ifelse(Mc %in% c("07", "08", "09"), "JulSep","OctFeb")))) %>% 
        # remove 1978, 79, 80,81; could leave 78 for spring stuff
        filter(Y != 1978 & Y != 1979 & Y != 1980 & Y != 1981) %>% 
        select(Date, Qm3day:gTPday, HighFlow, Y, M, TargetMonths) %>% 
        # sum to month group for low and high flow
        group_by(Y, TargetMonths, HighFlow) %>% 
        summarise(Qm3window = sum(Qm3day, na.rm = T),
                  gDRPwindow= sum(gSRPday, na.rm = TRUE),
                  gSSwindow = sum(gSSday, na.rm = TRUE),
                  gTPwindow = sum(gTPday, na.rm = TRUE))  %>% 
        # easy codes for ploting
        mutate(key = as.factor(paste0(TargetMonths,"_",HighFlow))) %>% 
        # mutate(key2 = fct_rev(fct_relevel(key, c("AugFeb_L75th","MarJun_L75th","Jul_L75th","AugFeb_G75th", "MarJun_G75th", "Jul_G75th")))) %>% 
        mutate(key2 = fct_rev(fct_relevel(key, c("OctFeb_L75th","MarJun_L75th","JulSep_L75th","OctFeb_G75th", "MarJun_G75th", "JulSep_G75th")))) %>% 
        mutate(DRP_TP = gDRPwindow/gTPwindow) %>% 
        filter(key == "MarJun_G75th") %>% 
        mutate(Qm3e6_MjHf = Qm3window/1e6,
               mtDRP_MjHf = gDRPwindow/1e6,
               mtSS_MjHf = gSSwindow/1e6,
               mtTP_MjHf = gTPwindow/1e6,
               Yn1975 = 1975,
               Yg = as.factor(ifelse(Y >=2003, "G2003", "L2003")),
               Yg = fct_relevel(Yg,"L2003", "G2003"),
               Y2 = str_sub(Y,3))

# raw bootstrapped scaled data, ----
# I just need gPsorbedgDM
# this is just the median g P sorbed per g DM, which was drawn from the empirical measurements
# takes a couple mins
MaumBSraw <- read.csv(file.path(here::here("04a_generatedDataTooBigForGit"), "04_RawBootstrapResults_Maumee.csv"))[-1,-1]
gPsorbedgDMmedian <- quantile(MaumBSraw$gPsorbedgDM, probs = 0.5)


# SS by Q and Y mod ----
# this model was created in 09_SStimeseriesFigS6
# SS ~ Q * Y
SSqy.glm <- readRDS(file.path(here::here("04a_generatedDataOnGit"), "09_SS_QQandYglm"))
summary(SSqy.glm)




# HOW MUCH DID DECLINES IN SS AND P SORB CONTRIBUTE TO CHANGES IN DRP ----

# SS PREDICTIONS ----
# SS predictions from model above
MR$mtSS_MjHf_pred = predict(SSqy.glm, MR)

plot(MR$mtSS_MjHf_pred ~ MR$mtSS_MjHf)

# SS for each yearly Q with 1975 SS concentrations (accounts for Q)
MR$mtSS_MjHf_w75seds = predict(SSqy.glm, MR %>% 
                                              mutate(Y = Yn1975))

# SS for each yearly Q with 2019 SS concentrations (accounts for Q)
MR$mtSSmarJun_w2019seds = predict(SSqy.glm, MR %>% 
                                       mutate(Y = 2019))

ggplot() +
         geom_point(data = MR, aes(y = mtSS_MjHf_w75seds, x = Qm3e6_MjHf), color = "red") +
          geom_point(data = MR, aes(y = mtSSmarJun_w2019seds, x = Qm3e6_MjHf), color = "blue") 

# P SORPTION PREDICTIONS ----
MR2 <- MR %>% 
                # P sorbed with predicted seds, using median gPsorb/gDM used in bootstrap
                # note on units: mult mt SS by proportion so final units are mt P
          mutate(mtPsorbP_MjHf_obs = mtSS_MjHf * gPsorbedgDMmedian,
                # P sorbed with 1975 predicted seds, using median gPsorb/gDM used in bootstrap
                 mtPsorbP_MjHf_w75seds = mtSS_MjHf_w75seds * gPsorbedgDMmedian,
                 # P sorbed with 2019 predicted seds, using median gPsorb/gDM used in bootstrap
                 mtPsorbP_MjHf_w2019seds = mtSSmarJun_w2019seds * gPsorbedgDMmedian,
                 # differences in sorption with different SS predictions
                 mtPsorbP_MjHf_obsM75 = mtPsorbP_MjHf_obs - mtPsorbP_MjHf_w75seds,
                 mtPsorbP_MjHf_19m75 = mtPsorbP_MjHf_w75seds - mtPsorbP_MjHf_w2019seds,
                 #  75-02 obs & 03-19 predicted with 75 loads (Psorb75 - PsorbPred)
                 mtDRP_MjHf_w75ss = ifelse(Yg == "L2003", mtDRP_MjHf, mtDRP_MjHf - (mtPsorbP_MjHf_w75seds - mtPsorbP_MjHf_obs)))

ggplot() +
  geom_point(data = MR2, aes(y = mtDRP_MjHf_w75ss, x = Y), color = "red") +
  geom_point(data = MR2, aes(y = mtDRP_MjHf, x = Y))

# Table 1 ----
## DRP~Q with observed loads ----
  # not normal
  qqnorm(MR2$mtDRP_MjHf); qqline(MR2$mtDRP_MjHf)
  # not bad
  qqnorm(MR2$Qm3e6_MjHf); qqline(MR2$Qm3e6_MjHf)
  
  ggplot(MR2, aes(y = mtDRP_MjHf,x = Qm3e6_MjHf, label = Y)) +
      geom_text()
  
  
  # Q X Y - MOST LIKELY MODEL
  DRP_Q_obs_glm_int <- glm(mtDRP_MjHf ~ Qm3e6_MjHf * Yg,
                           data = MR2,
                           family = Gamma(link = "identity")) #
  plot(DRP_Q_obs_glm_int)
  summary(DRP_Q_obs_glm_int)
  
  
  # Q + Y
  DRP_Q_obs_glm_add <- glm(mtDRP_MjHf ~ Qm3e6_MjHf + Yg,
                           data = MR,
                           family = Gamma(link = "identity"))
  plot(DRP_Q_obs_glm_add)
  summary(DRP_Q_obs_glm_add)
  
  # model comparison
  AICc(DRP_Q_obs_glm_int)
  AICc(DRP_Q_obs_glm_add)
  
  nagelkerke(DRP_Q_obs_glm_int)
  nagelkerke(DRP_Q_obs_glm_add)
  
  plot(residuals(DRP_Q_obs_glm_add) ~ MR2$Qm3e6_MjHf)
  
  # predictions for paper
  DRP_Q_obs_pred4000 <- as.data.frame(cbind(Qm3e6_MjHf = c(4000,4000), Yg = c("L2003","G2003")))
  DRP_Q_obs_pred4000$Qm3e6_MjHf <-  as.numeric(DRP_Q_obs_pred4000$Qm3e6_MjHf)
  DRP_Q_obs_pred4000$DRP_obs_pred <- predict(DRP_Q_obs_glm_int, DRP_Q_obs_pred4000)
  351.1412/223.3459
  
  ## DRP ~ Q with 75 loads  ----
  # Q x Y
  DRP_Q_obs_glm_75loads_int <- glm(mtDRP_MjHf_w75ss ~ Qm3e6_MjHf * Yg,
                                   data = MR2,
                                   family = Gamma(link = "identity"))
  summary(DRP_Q_obs_glm_75loads_int)

  
  # Q + Y - BETTER MODEL
  DRP_Q_obs_glm_75loads_add <- glm(mtDRP_MjHf_w75ss ~ Qm3e6_MjHf + Yg,
                                   data = MR2,
                                   family = Gamma(link = "identity"))
  summary(DRP_Q_obs_glm_75loads_add)

  # model comparison
  AICc(DRP_Q_obs_glm_75loads_int)
  AICc(DRP_Q_obs_glm_75loads_add)
  nagelkerke(DRP_Q_obs_glm_75loads_int)
  nagelkerke(DRP_Q_obs_glm_75loads_add)
  
  # predictions for paper
  DRP_Q_obs_pred4000$DRP_hist_pred <- predict(DRP_Q_obs_glm_75loads_add, DRP_Q_obs_pred4000)
  250.3-224.6
  250.3279/224.5751-1
  
  
  ## FIgure 4 ----
  # prepare data
  MR2_L03 <- MR2 %>% 
    filter(Y < 2003)
  
  MR2_L03$mtDRP_MjHf_pred <- (predict(DRP_Q_obs_glm_int, MR2_L03))
  
  MR2_G03 <- MR2 %>% 
    filter(Y >= 2003)
  MR2_G03$mtDRP_MjHf_pred <- ( predict(DRP_Q_obs_glm_int, MR2_G03))
  MR2_G03$mtDRP_MjHf_pred75ss <- (predict(DRP_Q_obs_glm_75loads_add, MR2_G03))
  
  # Figure
  Fig4 <-ggplot() +
    geom_line(data = MR2_L03, aes(y = mtDRP_MjHf_pred, x = Qm3e6_MjHf), color = "dodgerblue", size = 1.3) +
    geom_line(data = MR2_G03, aes(y = mtDRP_MjHf_pred75ss, x = Qm3e6_MjHf), color = "salmon", size = 1.3) +
    geom_line(data = MR2_G03, aes(y = mtDRP_MjHf_pred, x = Qm3e6_MjHf), color = "firebrick", size = 1.3) +
    geom_segment(data = MR2_G03, aes(y = mtDRP_MjHf, 
                                     yend = mtDRP_MjHf_w75ss + 5,
                                     x= Qm3e6_MjHf, xend = Qm3e6_MjHf),
                 arrow = arrow(length = unit(0.1, "inches")), color = "grey30", size = 0.75, alpha = 0.75) +
    geom_point(data = MR2_L03, aes(y = mtDRP_MjHf, x = Qm3e6_MjHf), fill = "dodgerblue", shape = 21, size = 4) +
    geom_point(data = MR2_G03, aes(y =   mtDRP_MjHf_w75ss, x = Qm3e6_MjHf), fill = "salmon", shape = 21, size = 4) +
    geom_point(data = MR2_G03, aes(y = mtDRP_MjHf, x = Qm3e6_MjHf), fill = "firebrick", shape = 21, size = 4)+
    scale_y_continuous("DRP load (tons P)") +
    xlab(expression(paste("Discharge (",m^3, "x ", 10^6,")"))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          panel.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 16),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 1, color = "black"),
          axis.line.y = element_line(size = 1, color = "black")) +
    geom_point(aes(x = 400, y = 400), shape = 21, fill = "firebrick", size = 3) +
    geom_point(aes(x = 400, y = 375), shape = 21, fill = "salmon", size = 3) +
    geom_point(aes(x = 400, y = 350), shape = 21, fill = "dodgerblue", size = 3) +
    annotate("text", x = 500, y = 400, label = "≥ 2003 observed", size = 7, color = "black", hjust = 0) +
    annotate("text", x = 500, y = 375, label = "≥ 2003 historic SS loads", size = 7, color = "black", hjust = 0) +
    annotate("text", x = 500, y = 350, label = "< 2003 observed", size = 7, color = "black", hjust = 0) 
  
  
  Fig4.g <- ggplotGrob(Fig4)
  Fig4.gtf <- gtable_frame(Fig4.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))
  
  png("05_Figures/11b_FigS9.png", units = "in", height = 7, width = 7, res = 300)
  grid.newpage()
  grid.draw(Fig4.gtf)
  dev.off()

# save ----
  # write.csv(MR2_G03, file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_G03.csv"))
  # write.csv(MR2_L03, file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_L03.csv"))
  # save.image(file.path(here::here("03_Rdata"), "11b_FigS9_Table1_0perDRP_Rdat"))
  # load(file.path(here::here("03_Rdata"), "11b_FigS9_Table1_0perDRP_Rdat"))
  