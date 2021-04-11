library(tidyverse)
library(MuMIn)
library(ggpubr)
library(tibble)
library(ggpmisc)
library(grid)
library(egg)
library(png)

# Import data -----------
# subset of data for HABS stuff 2003 to 2019
Mhab <- read.csv("04_generatedData/06d_Scale2WVhabs.csv", row.names = 1) %>% 
  mutate(var = as.factor(var))


Mall <- read.csv("04_generatedData/06d_Scale2WVall.csv", row.names = 1) %>% 
  mutate(Yg = as.factor(Yg)) %>% 
  # these were were summed over time intervals - doesn't make sense
  select(-c(gPsorbedgDM_2.5per, gPsorbedgDM_50per, gPsorbedgDM_97.5per)) %>% 
  mutate(gPsorbedgDM = gPsorbMarJunP_50per/gSSmarJun_50per)


# How much P was sorbed? -------------
summary(Mall$gPsorbMarJunP_50per)/1e6 #tons


# time series of DRP exports w and w/p sorption -----
# - Fig 4a 

Fig4A <- ggplot() +
  geom_ribbon(data = Mall, aes(ymin = gDRPmarJun_50per/1e6, ymax = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, x = Y), fill = "pink", alpha = 40/100) +
  geom_segment(data = Mall, aes(y = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, yend =  gDRPmarJun_50per/1e6+5, x = Y, xend=Y), 
               arrow = arrow(length = unit(0.1, "inches")), color = "black", size = 1.15)+
  geom_point(data = Mall, aes(y = gDRPmarJun_50per/1e6, x = Y), fill = "black", shape = 21, size = 3) +
  geom_point(data = Mall, aes(y = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, x = Y), fill = "grey", shape = 21, size = 3) +
  ylab("DRP load (tons P)") +
  xlab("Year") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black"),
        panel.background = element_rect(fill = "transparent")) +
  geom_point(aes(y =620, x = 1975), color = "grey", size = 5) +
  annotate("text", label = "Without P sorption", y = 620, x = 1976, size = 8, hjust = 0) +
  geom_point(aes(y =585, x = 1975), color = "black", size = 5) +
  annotate("text", label = "Observed", y = 585, x = 1976, size = 8, hjust = 0) 

  
# Cyano Fig - 4c -----


Fig4C <- ggplot(Mhab %>% 
                  filter(var == "cyanoIndex" | var == "cyanoIndexWoS") %>% 
                  mutate( Yn2 = str_sub(Yn,3),
                          var = fct_recode(var, "Observed" ="cyanoIndex", "Without P sorption" = "cyanoIndexWoS")),
                aes(y = values, x = var,label = Yn2))+
    geom_boxplot(fill = "grey80", alpha = 0.7) + #"#7fbf7b"
    geom_text(position = "jitter", size = 4) +
    scale_y_log10() +
    # annotate("text", x = 0.5, y = 90, label = "e", size = 8, fontface = "bold") +
    theme_bw() +
    ylab("Max cyanobacterial index") +
    xlab(NULL) +
    theme(axis.title = element_text(size = 24),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 20, angle = 25, hjust = 1),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(size = 1, color = "black"),
          axis.line.y = element_line(size = 1, color = "black"))

# * CyanoHABS effect size ----
  MhabW <- Mhab %>% 
    pivot_wider(names_from = "var", values_from = "values")
  
  MhabWs <- MhabW %>% 
    summarize(across(c(mtDRPmarJun_50per,mtDRPmarJun_50perWoS, cyanoIndex:cyanoIndexWoS), mean, na.rm = TRUE)) %>% 
    mutate(effectSize = 1-cyanoIndex/cyanoIndexWoS,
           HigherWSorption = mtDRPmarJun_50perWoS/mtDRPmarJun_50per)
  
# relationship between cyano Index and sorption
ggplot(Mhab %>% 
         pivot_wider(names_from = "var", values_from = "values") ,
       aes(y = cyanoIndex,x = mtPsorbMarJunP_50per)) +
  geom_point()

# What the equation is doing
ggplot(Mhab %>% 
         pivot_wider(names_from = "var", values_from = "values") ,
       aes(y = cyanoIndex,x = mtDRPmarJun_50per, label = Yn)) +
  geom_text()



################
#  SS drives P sorption --------
################

# Get data together
Mall2 <- Mall %>% 
  mutate(Yg2 = ifelse(Yn < 2003, "<2003", ">2003"),
         mtSSmarJun_50per = gSSmarJun_50per/1e6, # tonne or 1 billion g
         Qm3e6_marJun_50per = Qm3marJun_50per/1e6,
         mtDRPmarJun_50per = gDRPmarJun_50per/1e6,
         L10.mtSSmarJun_50per = log10(mtSSmarJun_50per),
         L10.Qm3e6_marJun_50per = log10(Qm3e6_marJun_50per),
         L10.mtDRPmarJun_50per = log10(mtDRPmarJun_50per),
         Yn2 = str_sub(Yn,3))


# Differences before and after 2003
Mall2 %>% 
  group_by(Yg2) %>% 
  summarize_at(vars(perSorb2_50per), list(mean = mean, sd = sd), na.rm = TRUE)

# Teasing apart relationships among Q, SS, DRP, Sorp, and time ------
# DISCHARGE v TIME

ggplot(Mall2, aes(y = L10.Qm3e6_marJun_50per, x = Y)) +
  geom_point()

QvY.l <- lm(L10.Qm3e6_marJun_50per ~ Y, Mall2); 
summary(QvY.l) # year is marginal


QvY1985.l <- lm(L10.Qm3e6_marJun_50per ~ Y, Mall2 %>% 
              filter(Yn >= 1985)); 
summary(QvY1985.l) # year is now significant


# SUSPENDED SEDIMENTS v Q --------

# Decline in SS loads over time
ggplot(Mall2, aes(y = mtSSmarJun_50per, x = Yn)) +
  geom_point() 

# best model - log10
SSvQintYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per * Yn, Mall2); 
SSvQpYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yn, Mall2)
SSvQpYg2.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yg2, Mall2)
SSvQ.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall2)
SSvt.l <- lm(L10.mtSSmarJun_50per ~ Yn, Mall2)

# what's the best model - Q+Yn, Q+Yg2
model.sel(SSvQintYn.l, SSvQpYn.l, SSvQ.l, SSvQpYg2.l)
# Best model is SS ~ Q + Y
summary(SSvQpYn.l)
# plot(SSvQpYn.l)

# Look at how much SS has changed over time by focusing on one med-high Q
# declines in SS at 4000
# predicts SS in metric tons
  SSAt4000.1975 <- 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = log10(4000), Yn = 1975))
  SSAt4000.2019 <- 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = log10(4000), Yn = 2019))
  # percent decline in SS
  1-SSAt4000.2019/SSAt4000.1975
  
  # decline in SS per decade - kiloton
  (SSAt4000.1975 - SSAt4000.2019) /(2019-1975)*10/1e3
  
  
  SSAt4000.1975 - SSAt4000.2019

# declines in P sorption at 4000
  # predicts metric tons P
  PsorbAt4000.1975 <- SSAt4000.1975 * quantile(Mall2$gPsorbedgDM, probs = 0.5)
  PsorbAt4000.2019 <- SSAt4000.2019 * quantile(Mall2$gPsorbedgDM, probs = 0.5)
  # percent decline in P sorption
  1-PsorbAt4000.2019/PsorbAt4000.1975
  # per decade decline in P sorption metric tons P
  (PsorbAt4000.1975 - PsorbAt4000.2019)/(2019-1975)*10
  
    
# residuals from SS ~ Q relationship
  # checked with gam  - linear best model for time
Mall2$SS_Qres <- residuals(SSvQ.l)
summary(lm(SS_Qres ~ Yn, Mall2))


ggplot(Mall2, aes(y = SS_Qres, x = Yn, label= Yn)) +
  geom_text() +
  stat_smooth(method = "lm")

# DRP v Q -------
DRPvQintYn.l1 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per * Yn, Mall2)
DRPvQintYn.l2 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per + Yn, Mall2)
DRPvQintYn.l3 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall2)
# just curious if SS did better than Y, but they are equal
DRPvQintYn.l4 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per * SS_Qres, Mall2)
DRPvQintYn.l5 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per + SS_Qres, Mall2)
model.sel(DRPvQintYn.l1, DRPvQintYn.l2, DRPvQintYn.l3, DRPvQintYn.l4, DRPvQintYn.l5) # lowest AICc is additive model

# best time model DRP ~ Q + Y
summary(DRPvQintYn.l2) # effect of year is marginally significant

# best SS model DRP ~ Q + SSres
summary(DRPvQintYn.l5) # effect of SS is marginally significant

# SOME SUPPLEMENTAL PLOTS FIG S5 -------
Qvtime <- ggplot() +
  geom_line(data = Mall2, aes(y = Qm3e6_marJun_50per, x = Y),  color = "black", size = 0.5, linetype = "dashed") +
  geom_point(data = Mall2, aes(y = Qm3e6_marJun_50per, x = Y), shape = 21, fill = "grey80", size = 4) +
  stat_smooth(data = Mall2 %>% 
                filter(Yn >= 1985), aes(y = Qm3e6_marJun_50per, x = Y), method = "lm", color = "black") +
  scale_y_log10(breaks = c(1000, 2000, 3000, 4000, 5000)) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  ylab(expression(paste("Discharge (",m^3," x ",10^6,")"))) +
  xlab("Year") +
  annotate("text", x = 2006, y = 910, label = "R^2 == 0.24;", parse = T, size = 7) +
  annotate("text", x = 2015, y = 900, label = "P = 0.003", size = 7) 

SSvtime <-  ggplot() +
  geom_line(data = Mall2, aes(y = mtSSmarJun_50per/1e6, x = Y),  color = "black", size = 0.5, linetype = "dashed") +
  geom_point(data = Mall2, aes(y = mtSSmarJun_50per/1e6, x = Y), shape = 21, fill = "grey80", size = 4) +
  # stat_smooth(data = Mall2, aes(y = mtSSmarJun_50per, x = Y), method = "lm")+
  scale_y_log10()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  ylab(expression(paste(atop("Suspended sediment load", "(Megaton dry mass)")))) +
  xlab("Year") +
  annotate("text", x = 2015, y = 0.125, label = "R^2 == 0.02", parse = TRUE, size = 7) +
  annotate("text", x = 2015, y = 0.1, label = "P = 0.43", size = 7)


SSvQ <- ggplot() +
  # geom_point(data = Mall2, aes(y = mtSSmarJun_50per/1e6, x = Qm3e6_marJun_50per), shape = 21, size = 4, fill = "green") +
  geom_text(data = Mall2, aes(y = mtSSmarJun_50per/1e6, x = Qm3e6_marJun_50per, label = Yn2), size = 5) +
  stat_smooth(data = Mall2, aes(y = mtSSmarJun_50per/1e6, x = Qm3e6_marJun_50per), method = "lm", color = "black")+
  scale_y_log10(breaks = c(0.1, 0.25, 0.5, 0.75, 1)) +
  scale_x_log10(breaks = c(1000, 2000, 3000, 4000, 5000))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  ylab(expression(paste(atop("Suspended sediment load", "(Megaton dry mass)")))) +
  xlab(expression(paste("Discharge (",m^3," x ",10^6,")"))) +
  annotate("text", x = 3500, y = 0.145, label = "R^2 == 0.69", parse = TRUE, size = 7) +
  annotate("text", x = 3500, y = 0.12, label = "P[Q] < 0.001", parse = TRUE, size = 7) +
  annotate("text", x = 3500, y = 0.1, label = "P[year] < 0.001", parse = TRUE, size = 7)


SS.QresVtime <- ggplot() +
  geom_line(data = Mall2, aes(y = SS_Qres, x = Yn), color = "black", size = 0.5, linetype = "dashed") +
  geom_point(data = Mall2, aes(y = SS_Qres, x = Yn), shape = 21, fill = "grey80", size = 4) +
  stat_smooth(data = Mall2, aes(y = SS_Qres, x = Yn), 
              method = "lm", color = "black") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  ylab(expression(atop("SS ~ Q residuals", paste("(", log[10]," megaton dry mass)")))) +
  xlab("Year") +
  annotate("text", x = 2015, y = -0.29, label = "R^2 == 0.27", parse = TRUE, size = 7) +
  annotate("text",x = 2015, y = -0.35, label = "P < 0.001", size = 7) +
  ylim(-0.4,0.4) 



# Make plot
FigS5A.g <- ggplotGrob(Qvtime)
FigS5B.g <- ggplotGrob(SSvtime)
FigS5C.g <- ggplotGrob(SSvQ)
FigS5D.g <- ggplotGrob(SS.QresVtime)


FigS5A.gtf <- gtable_frame(FigS5A.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
FigS5B.gtf <- gtable_frame(FigS5B.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
FigS5C.gtf <- gtable_frame(FigS5C.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
FigS5D.gtf <- gtable_frame(FigS5D.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))

FigS5AB.gtf <- gtable_frame(gtable_cbind(FigS5A.gtf, FigS5B.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
FigS5CD.gtf <- gtable_frame(gtable_cbind(FigS5C.gtf, FigS5D.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))


FigS5.gtf <- gtable_rbind(FigS5AB.gtf, FigS5CD.gtf)




png("05_Figures/08_FigS5.png", units = "in", height = 12, width = 16, res = 300)
grid.newpage()
grid.draw(FigS5.gtf)
grid.text("a", x = unit(0.03,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.53,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.03,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("d", x = unit(0.53,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()


# HOW MUCH DID DECLINES IN SS AND P SORB CONTRIBUTE TO CHANGES IN DRP

  # make new dataframe
  Mall3 <- Mall2
  Mall3$Yn1975 <-  1975
  
  # # SS PREDICTIONS
  # # SS predictions from model above
  # Mall3$mtSSmarJun_pred = 10^predict(SSvQpYn.l, Mall3)
  # 
  # # SS for each yearly Q with 1975 concentrations (accounts for Q)
  # Mall3$mtSSmarJun_w75seds = 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = Mall3$L10.Qm3e6_marJun_50per, 
  #                                                             Yn = rep(1975, length = 42)))
  # 
  # # SS for each yearly Q with 2019 concentrations (accounts for Q)
  # Mall3$mtSSmarJun_w2019seds = 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = Mall3$L10.Qm3e6_marJun_50per, 
  #                                                               Yn = rep(2019, length = 42)))
  # 
  # # P SORPTION PREDICTIONS
  # # P sorbed with predicted seds, using median from median values from monte carlo
  # Mall3$mtPsorbMarJunP_wPredSeds = Mall3$mtSSmarJun_pred* quantile(Mall2$gPsorbedgDM, probs = 0.5)
  # 
  # # P sorbed with 75 seds
  # Mall3$mtPsorbMarJunP_w75seds = Mall3$mtSSmarJun_w75seds* quantile(Mall2$gPsorbedgDM, probs = 0.5)
  # 
  # # P sorbed with 2019 seds
  # Mall3$mtPsorbMarJunP_w2019seds = Mall3$mtSSmarJun_w2019seds* quantile(Mall2$gPsorbedgDM, probs = 0.5)
  # 
  # # DIFFERENCES BETWEEN DIFFERENT P SORPTION PREDICTIONS
  # # Dif between P sorb with 1975 and year "y" seds
  # Mall3$mtPsorbMarJunP_75_Ydiff = Mall3$mtPsorbMarJunP_w75seds - Mall3$mtPsorbMarJunP_wPredSeds
  # 
  # # Dif between P sorb with 1975 and 2019 seds
  # Mall3$mtPsorbMarJunP_75_19diff = Mall3$mtPsorbMarJunP_w75seds - Mall3$mtPsorbMarJunP_w2019seds

# PREPARE DATA FOR PLOT
  # GET DRP -Q PREDICTIONS FOR < 2003 AND >2003
  #< 2003
  Mall_L03 <- Mall3 %>% 
    filter(Yn < 2003)
  
  DRP_Q_lm_L03 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall_L03)
  summary(DRP_Q_lm_L03)
  Mall_L03$mtDRPmarJun_predict <- 10^predict(DRP_Q_lm_L03, Mall_L03)
  
  # > 2003
  Mall_G03_0 <- Mall3 %>% 
    mutate(mtDRPmarJun_w75SSa = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff, # difference between 75 and year y
           L10.mtDRPmarJun_w75SSa = log10(mtDRPmarJun_w75SSa))
  
  Mall_G03 <- Mall_G03_0%>% 
    filter(Yn >= 2003) 
  
  # DRP ~ Q - observed
  DRP_Q_lm_G03 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall_G03)
  summary(DRP_Q_lm_G03)
  Mall_G03$mtDRPmarJun_predict <- 10^predict(DRP_Q_lm_G03, Mall_G03)
  
  # DRP ~ Q - 1975 SS - option 1 1975 - year "y"
  DRP_Q_lm_L03w75SSa.lm <- lm(L10.mtDRPmarJun_w75SSa ~ L10.Qm3e6_marJun_50per, Mall_G03)
  summary(DRP_Q_lm_L03w75SSa.lm)
  Mall_G03$mtDRPmarJun_predict_w75SSa <- 10^predict(DRP_Q_lm_L03w75SSa.lm, Mall_G03)
  
  # DOES < V. > 2003 DIFFER?
  Mall_G03_2 <- Mall_G03 %>% 
    select(Yn, L10.Qm3e6_marJun_50per, L10.mtDRPmarJun_50per, L10.mtDRPmarJun_w75SSa) %>% 
    # group_by(Yn, L10.Qm3e6_marJun_50per) %>% 
    pivot_longer(c(L10.mtDRPmarJun_50per, L10.mtDRPmarJun_w75SSa), names_to = "DRPval", values_to = "L10.mtDRPmarJun")
  
  
  
  # 1- < V. > 2003 DRP ~ Q different?
  DRP_Q_prepost2003int.lm <-  lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per * Yg2, Mall3)
  summary(DRP_Q_prepost2003int.lm) #NS
  
  DRP_Q_prepost2003plus.lm <-  lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per + Yg2, Mall3)
  summary(DRP_Q_prepost2003plus.lm) # plus 2003 sig
  
  
  # 2 - DRP ~ Q relationship differ based on < 2003 and > 2003 with historic SS
  Mall_HistPred <- Mall_G03_0 %>% 
    select(Yn, L10.Qm3e6_marJun_50per, L10.mtDRPmarJun_50per, L10.mtDRPmarJun_w75SSa, Yg2) %>% 
    mutate(L10.mtDRPmarJunHistPred = ifelse(Yn < 2003, L10.mtDRPmarJun_50per, L10.mtDRPmarJun_w75SSa)) %>% 
    select(Yn, Yg2, L10.Qm3e6_marJun_50per, L10.mtDRPmarJunHistPred)
  
  Mall_HistPred_int.lm <- lm(L10.mtDRPmarJunHistPred ~ L10.Qm3e6_marJun_50per * Yg2, data = Mall_HistPred)
  summary(Mall_HistPred_int.lm) #NS
  
  Mall_HistPred_plus.lm <- lm(L10.mtDRPmarJunHistPred ~ L10.Qm3e6_marJun_50per + Yg2, data = Mall_HistPred)
  summary(Mall_HistPred_plus.lm) #NS
  

  
  
  # # DRP ~ Q - 1975 SS - option 1 1975 - year "y"
  # DRP_Q_lm_L03w75SSb.lm <- lm(L10.mtDRPmarJun_w75SSb ~ L10.Qm3e6_marJun_50per, Mall_G03)
  # summary(DRP_Q_lm_L03w75SSb.lm)
  # Mall_G03$mtDRPmarJun_predict_w75SSb <- 10^predict(DRP_Q_lm_L03w75SSb.lm, Mall_G03)
  
  # difference in DRP at 4000 between before and after 2003
  DRPloadAt4000.L03 <- 10^predict(DRP_Q_lm_L03, data.frame(L10.Qm3e6_marJun_50per = log10(4000)))
  DRPloadAt4000.G03 <- 10^predict(DRP_Q_lm_G03, data.frame(L10.Qm3e6_marJun_50per = log10(4000)))
  CyanoHABsAt4000.L03 <- 0.48 * 10^(DRPloadAt4000.L03 * 0.00387)
  CyanoHABsAt4000.G03 <- 0.48 * 10^(DRPloadAt4000.G03 * 0.00387)
  DRPloadAt4000.G03/DRPloadAt4000.L03
  DRPloadAt4000.G03 - DRPloadAt4000.L03
  CyanoHABsAt4000.G03/CyanoHABsAt4000.L03
  

# impact of historical flows on DRP ~ Q relationship
# NOT SURE WHAT THIS IS DOING
# DRP_Q.lm <- lm(log10(gDRPmarJun_50per) ~ log10(Qm3marJun_50per), Mall3); summary(DRP_Q.lm)
# DRP_QYg2.lm  <- lm(log10(gDRPmarJun_50per) ~ log10(Qm3marJun_50per) + Yg2, Mall3); summary(DRP_QYg2.lm)
# 
# Mall4 <- Mall3 %>% 
#   mutate(mtDRPmarJun_w75SSa = ifelse(Y <2003, mtDRPmarJun_50per, mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff),
#          mtDRPmarJun_w75SSb = ifelse(Y < 2003,mtDRPmarJun_50per ,mtDRPmarJun_50per - mtPsorbMarJunP_75_19diff),
#          L10.mtDRPmarJun_w75SSa = log10(mtDRPmarJun_w75SSa),
#          L10.mtDRPmarJun_w75SSb = log10(mtDRPmarJun_w75SSb))
# 
# DRP_QYg2.lm  <- lm(log10(mtDRPmarJun_w75SSa) ~ log10(Qm3marJun_50per) + Yg2, Mall4); summary(DRP_QYg2.lm)
  
  
# PLOT
 Fig4D <-  ggplot()+
    geom_segment(data = Mall_G03, aes(y = mtDRPmarJun_50per, 
                                      yend = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff,
                                      x= Qm3e6_marJun_50per, xend = Qm3e6_marJun_50per),
                 arrow = arrow(length = unit(0.1, "inches")), color = "grey30", size = 0.75, alpha = 0.75) +
    # geom_text(data = Mall_G03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per, label = Yn2), color = "firebrick", fontface = "bold") +
    geom_point(data = Mall_G03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per), color = "firebrick", fontface = "bold", size = 4) +
    geom_point(data = Mall_G03, aes(y = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff, x = Qm3e6_marJun_50per), 
               fill = "salmon", shape = 21, size = 4) +
    geom_point(data = Mall_L03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per), 
               fill = "dodgerblue", shape = 21, size = 4) +
    geom_line(data = Mall_G03, aes(y = mtDRPmarJun_predict, x = Qm3e6_marJun_50per), color = "firebrick", size = 1.3) +
    geom_line(data = Mall_G03, aes(y = mtDRPmarJun_predict_w75SSa, x = Qm3e6_marJun_50per), color = "salmon", size = 1.3)+
    geom_line(data = Mall_L03, aes(y = mtDRPmarJun_predict, x = Qm3e6_marJun_50per), color = "dodgerblue", size = 1.3) +
    scale_y_continuous("DRP load (tons P)",
                       sec.axis = sec_axis(~ 0.48 * 10^(. * 0.00387), #conver to kT DM then metric ton P
                                           "Max cyanobacteria index", breaks = c(1, 2.5, 5, 10, 20, 30))) +
    xlab(expression(paste("Discharge (",m^3, "x ", 10^6,")"))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          panel.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 16),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 1, color = "black"),
          axis.line.y = element_line(size = 1, color = "black")) +
    geom_point(aes(x = 1000, y = 450), shape = 21, fill = "firebrick", size = 3) +
    geom_point(aes(x = 1000, y = 425), shape = 21, fill = "salmon", size = 3) +
    geom_point(aes(x = 1000, y = 400), shape = 21, fill = "dodgerblue", size = 3) +
    annotate("text", x = 1100, y = 450, label = "≥ 2003 observed", size = 7, color = "black", hjust = 0) +
    annotate("text", x = 1100, y = 425, label = "≥ 2003 historic SS loads", size = 7, color = "black", hjust = 0) +
    annotate("text", x = 1100, y = 400, label = "< 2003 observed", size = 7, color = "black", hjust = 0) 




# check to make sure relationship between res and year is linear
# Man, that's legit!
Mall4 <- Mall3 %>% 
  mutate(SSvQ.lres = 10^residuals(SSvQ.l),
         SSflowWeightedMean_mgL = log10(mtSSmarJun_50per/Qm3e6_marJun_50per*1000))



Fig4Dinset <- ggplot(Mall4, aes(y = SSvQ.lres, x = Yn)) +
  stat_smooth(method = "lm", color = "black", fill = "grey80") +
  geom_point(shape = 21, size = 2, fill = "grey30") +
  ylab(expression(atop("Suspended sediment v.", "Discharge residuals (tons)"))) +
  xlab("Year") +
  # scale_y_log10() +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 9),
        panel.border = element_rect(size = 1, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

################
# DRP sorp v. DRP load
################

sorpRange0 <- aes(ymin = gPsorbMarJunP_2.5per/1e6, ymax = gPsorbMarJunP_97.5per/1e6, x = gDRPmarJun_50per/1e6)
Fig4B <- ggplot(Mall2, aes(y = gPsorbMarJunP_50per/1e6, x = gDRPmarJun_50per/1e6,  fill = Yg2)) +
  geom_errorbar(sorpRange0) +
  geom_point(size = 4, shape = 21) +
  # geom_text() +
  geom_segment(aes(x = 0, xend = 300, y = 0, yend = 300), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 500, y = 0, yend = 500*0.5), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 500, y = 0, yend = 500*0.25), color = "grey", linetype = "dashed") +
  annotate("text", x = 280, y = 300, label = "100%", color = "grey40", size = 6) +
  annotate("text", x = 500, y = 260, label = "50%", color = "grey40", size = 6) +
  annotate("text", x = 500, y = 135, label = "25%", color = "grey40", size = 6) +
  # annotate("text", x = 0, y = 300, label = "d", size = 8, fontface = "bold") +
  scale_fill_manual(values = c("dodgerblue", "firebrick"), name = NULL, labels = c("< 2003","≥ 2003")) +
  # scale_color_manual(values = c("firebrick", "dodgerblue"), name = "Year", labels = c("≤ 2002","≥ 2003")) +
  # geom_abline(intercept = 0, slope =1, color = "grey") +
  # geom_abline(intercept = 0, slope =0.5, color = "grey") +
  # geom_abline(intercept = 0, slope = 0.25, color = "grey") +
  theme_bw() +
  xlab("DRP load (tons P)") +
  ylab("P sorption (tons P)") +
  ylim(0,300) +
  xlim(0,500) +
  theme(
    legend.position = c(0.17,0.9),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 21),
    legend.background = element_rect(fill = "transparent"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_line(size = 1, color = "black")
  )

################
# Change in % DRP sorbed b/w <2003 and >=2003
################
Mall2 %>% 
  group_by(Yg2) %>% 
  summarize(perSorb2_50per = mean(perSorb2_50per, na.rm = T))
0.916/0.332

Mall2 %>% 
  group_by(Yg2) %>% 
  summarize(perSorb2_50per_sd = sd(perSorb2_50per, na.rm = T))

# Q
perSorb2_50per
ggplot(Mall2, aes(y = perSorb2_50per, x = log(Qm3marJun_50per), color = Yg2)) +
  geom_point()
perSvQpY <- lm(perSorb2_50per*100 ~ log(Qm3marJun_50per) + Yg2, data = Mall2)
summary(perSvQpY)
cor(Mall2$perSorb2_50per,Mall2$Qm3marJun_50per)


# SS
cor(Mall2$perSorb2_50per,Mall2$gSSmarJun_50per)
ggplot(Mall2, aes(y = perSorb2_50per, x = log(gSSmarJun_50per), color = Yg2)) +
  geom_point()

perSvSSpY <- lm(perSorb2_50per*100 ~ log(gSSmarJun_50per) + Yg2, data = Mall2)
summary(perSvSSpY)

# DRP
ggplot(Mall2, aes(y = perSorb2_50per, x = log(gDRPmarJun_50per), color = Yg2)) +
  geom_point()
cor(Mall2$perSorb2_50per,Mall2$gDRPmarJun_50per)
perSvDRPpY <- lm(perSorb2_50per*100 ~ log(gDRPmarJun_50per) + Yg2, data = Mall2)
summary(perSvDRPpY)


###############
# BLAH
###############
# change in Q over time
# NEED THIS FIGURE WITH 85-PRES LINE
ggplot(Mall2, aes(y = log10(Qm3marJun_50per), x = Yn)) +
  geom_point() 

summary(lm(log10(Qm3marJun_50per) ~ Yn, Mall2))
summary(lm(log10(Qm3marJun_50per) ~ Yn, Mall2 %>% 
                                            filter(Mall2 >= 1985)))


################
# Assemble PLOT
################



Fig4A.g <- ggplotGrob(Fig4A)
Fig4B.g <- ggplotGrob(Fig4B)
Fig4C.g <- ggplotGrob(Fig4C)
Fig4D.g <- ggplotGrob(Fig4D)


Fig4A.gtf <- gtable_frame(Fig4A.g,  width = unit(1.8, "null"), height = unit(0.9, "null"))
Fig4B.gtf <- gtable_frame(Fig4B.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))
Fig4C.gtf <- gtable_frame(Fig4C.g,  width = unit(0.4, "null"), height = unit(0.9, "null"))
Fig4D.gtf <- gtable_frame(Fig4D.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))

# Fig4ABC.gtf <- gtable_frame(gtable_cbind(Fig4AB.gtf, Fig4C.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
# Fig4AB.gtf <- gtable_frame(gtable_cbind(Fig4A.gtf, Fig4B.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4BC.gtf <- gtable_frame(gtable_cbind(Fig4B.gtf, Fig4C.gtf, Fig4D.gtf),  width = unit(1.8, "null"), height = unit(0.9, "null"))

Fig4.gtf <- gtable_rbind(Fig4A.gtf, Fig4BC.gtf)




png("05_Figures/08_Fig4_2.png", units = "in", height = 12, width = 16, res = 300)
grid.newpage()
grid.draw(Fig4.gtf)
grid.text("a", x = unit(0.01,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.01,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.385,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("d", x = unit(0.615,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()


#####################
# SAVE/LOAD
#####################
#didn't get through everything before saving
# save.image("03_Rdata/08_MaumeeFigs3_rdat")
# load("03_Rdata/08_MaumeeFigs3_rdat")


  