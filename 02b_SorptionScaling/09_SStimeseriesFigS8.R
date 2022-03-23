# JMH Feb 2022
# Generates Fig S8 (changes in discharge and suspended sediment load during 1975-2019). 
# Also makes some calculations for the Results.
# this was Fig. S6, hence the names

# library ----
library(tidyverse)
library(MuMIn)
library(grid)
library(egg)
library(jtools)

# data ----
# prepared in 02d_WQandQ_MaumeeWaterville
# took Heid data, which is generally sub daily, then combined with USGS daily discharge
# (USGS 15 min discharge only goes back to 1987)
# built GAMS to predict missing days from Q
# No Heid monitoring in late 70s, early 80s which I've removed
# This is my SmpTimeWindowDay, not Heid's which I didn't understand
maum <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "02d_MaumeeWatervilleWaterQual.csv"), row.names = 1) %>% 
            # filter to Mar-June and high flow
            filter(M >= 3 & M < 7) %>% 
            filter(HighFlow == "HighFlow")%>% 
            # sum loads to Mar-June high flows
            group_by(Y) %>% 
            summarise(Qm3window = sum(Qm3day, na.rm = T),
                      gDRPwindow= sum(gSRPday, na.rm = TRUE),
                      gSSwindow = sum(gSSday, na.rm = TRUE),
                      gTPwindow = sum(gTPday, na.rm = TRUE)) %>% 
            # some conversions
            mutate(Qm3e6_MjHf = Qm3window/1e6,
                   mtDRP_MjHf = gDRPwindow/1e6,
                   mtSS_MjHf = gSSwindow/1e6, # tonne or 1 billion g
                   mtTP_MjHf = gTPwindow/1e6,
                   L10.Qm3e6_MjHf = log10(Qm3e6_MjHf),
                   L10.mtDRP_MjHf = log10(mtDRP_MjHf),
                   L10.mtSS_MjHf = log10(mtSS_MjHf),
                   L10.mtTP_MjHf = log10(mtTP_MjHf),
                   Yn2 = str_sub(Y,3))


# raw bootstrapped scaled data, I just need gPsorbedgDM
  # this is just the median g P sorbed per g DM, which was drawn from the empirical measurements
  # takes a couple mins
  MaumBSraw <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "04_RawBootstrapResults_Maumee.csv"))[-1,-1]
  gPsorbedgDMmedian <- quantile(MaumBSraw$gPsorbedgDM, probs = 0.5)


# SS time series analysis ----
  # Have to account for Q
  # evaluate normality and dispersion
  # not bad, log makes this worse
  qqnorm(maum$mtSS_MjHf); qqline(maum$mtSS_MjHf)
  qqnorm(maum$Qm3e6_MjHf); qqline(maum$Qm3e6_MjHf)
  
  ggplot(maum, aes(y = L10.mtSS_MjHf, x = L10.Qm3e6_MjHf)) +
    geom_point() 
  
  ggplot(maum, aes(y = mtSS_MjHf, x = Qm3e6_MjHf)) +
    geom_point() +
    stat_smooth(method = "lm")

  # model selection
  # checked with gam,  converges on linear model
  SSmod1 <- glm(mtSS_MjHf ~ Qm3e6_MjHf * Y, 
                data = maum,
                family = Gamma(link = "identity"),
                )
  SSmod2 <- glm(mtSS_MjHf ~ Qm3e6_MjHf + Y, 
                data = maum,
                family = Gamma(link = "identity"))
  SSmod3 <- glm(mtSS_MjHf ~ Qm3e6_MjHf , 
                data = maum,
                family = Gamma(link = "identity"))
  SSmod4 <- glm(mtSS_MjHf ~ Y, 
                data = maum,
                family = Gamma(link = "identity"))
  

  # what's the best model?
  model.sel(SSmod1, SSmod2, SSmod3, SSmod4)

  # Best model is SS ~ Q * Y
  summary(SSmod1)
  plot(SSmod1)
  summ(SSmod1)
  
  
  # Q relationship
  effect_plot(SSmod1, pred = Qm3e6_MjHf, interval = T, plot.points = TRUE)
  # raw relationship
  effect_plot(SSmod1, pred = Y, interval = T, plot.points = TRUE)
  # after controling for Q
  effect_plot(SSmod1, pred = Y, interval = T, plot.points = TRUE, partial.residuals = T)
  
# How much has SS changed over time? ----
  # predicts SS in metric tons
  maum$mtSSAt500Qm3e6pred <- predict(SSmod1, data.frame(Qm3e6_MjHf = 500, Y = maum$Y))
  maum$mtSSAt1000Qm3e6pred <- predict(SSmod1, data.frame(Qm3e6_MjHf = 1000, Y = maum$Y))
  maum$mtSSAt2000Qm3e6pred <- predict(SSmod1, data.frame(Qm3e6_MjHf = 2000, Y = maum$Y))
  maum$mtSSAt3000Qm3e6pred <- predict(SSmod1, data.frame(Qm3e6_MjHf = 3000, Y = maum$Y))
  maum$mtSSAt4000Qm3e6pred <- predict(SSmod1, data.frame(Qm3e6_MjHf = 4000, Y = maum$Y))
  
  ggplot() +
    geom_point(data= maum, aes(y = mtSSAt500Qm3e6pred, x = Y), color = "red") +
    geom_point(data= maum, aes(y = mtSSAt1000Qm3e6pred, x = Y), color = "pink") +
    geom_point(data= maum, aes(y = mtSSAt2000Qm3e6pred, x = Y), color = "orange") +
    geom_point(data= maum, aes(y = mtSSAt3000Qm3e6pred, x = Y), color = "lightblue") +
    geom_point(data= maum, aes(y = mtSSAt4000Qm3e6pred, x = Y), color = "blue")
  
  SSAt4000.1975 <- predict(SSmod1, data.frame(Qm3e6_MjHf = 4000, Y = 1975))
  SSAt4000.2019 <- predict(SSmod1, data.frame(Qm3e6_MjHf = 4000, Y = 2019))
  
  # total decline in SS (mT SS)
  SSAt4000.1975 - SSAt4000.2019
  
  # percent decline in SS
  (1-SSAt4000.2019/SSAt4000.1975)*100
  
  # decline in SS per decade - kiloton
  (SSAt4000.1975 - SSAt4000.2019) /(2019-1975)*10/1e3
  
## Associated decline in P sorption? ----
  # change at Q=4000
  # predicts metric tons P
  PsorbAt4000.1975 <- SSAt4000.1975 * gPsorbedgDMmedian
  PsorbAt4000.2019 <- SSAt4000.2019 * gPsorbedgDMmedian
  
  # percent decline in P sorption
  1-PsorbAt4000.2019/PsorbAt4000.1975
  
  # per decade decline in P sorption metric tons P
  (PsorbAt4000.1975 - PsorbAt4000.2019)/(2019-1975)*10 
  

# Fig S8 ----
  
  ## Fig S8a -----
  summary(lm(Qm3e6_MjHf ~ Y, maum[maum$Y >= 1985,]))
  
  FigS6a <- ggplot() +
    geom_line(data = maum, aes(y = Qm3e6_MjHf, x = Y),  color = "black", size = 0.5, linetype = "dashed") +
    geom_point(data = maum, aes(y = Qm3e6_MjHf, x = Y), shape = 21, fill = "grey80", size = 4) +
    stat_smooth(data = maum %>% 
                  filter(Y >= 1985), aes(y = Qm3e6_MjHf, x = Y), method = "lm", color = "black") +
    scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000), limits = c(0,5000)) +
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
    annotate("text", x = 2006, y = 0, label = "R^2 == 0.30;", parse = T, size = 7) +
    annotate("text", x = 2015, y = 0, label = "P < 0.001", size = 7) 
  
  ## Fig S8b ----
  summary(lm(mtSS_MjHf ~ Y, maum))
  
  FigS6b <-  ggplot() +
    geom_line(data = maum, aes(y = mtSS_MjHf/1e6, x = Y),  color = "black", size = 0.5, linetype = "dashed") +
    geom_point(data = maum, aes(y = mtSS_MjHf/1e6, x = Y), shape = 21, fill = "grey80", size = 4) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), lim = c(0,1.5)) +
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
    annotate("text", x = 1992, y = 0.0175, label = "R^2 == 0", parse = TRUE, size = 7, hjust = 0) +
    annotate("text", x = 1997.5, y = 0, label = "; P = 0.81", size = 7, hjust = 0)
  
  
    SSQpred <- as.data.frame(cbind(Qm3e6_MjHf = seq(360, 4600, by = 10),
                          Y = rep(c(1975,2019), each = 425))) 
    SSQpred$mtSS_MjHf_pred = predict(SSmod1, SSQpred)

  
  # Fig S8c ----
    summ(SSmod1)
    
  
  FigS6c <- ggplot() +
      # geom_text(data = maum, aes(y = mtSS_MjHf/1e6, x = Qm3e6_MjHf, label = Yn2), size = 5) +
      geom_line(data = SSQpred, aes(y = mtSS_MjHf_pred/1e6, x = Qm3e6_MjHf, linetype = as.factor(Y))) +
      # geom_point(data = maum, aes(y = mtSS_MjHf/1e6, x = Qm3e6_MjHf), shape = 21, fill = "grey80", size = 4) +
      geom_text(data = maum, aes(y = mtSS_MjHf/1e6, x = Qm3e6_MjHf, label = Yn2), size = 6) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2), lim = c(0,2)) +
      scale_x_continuous(breaks = c(0, 2000, 4000, 6000), lim = c(0,6000))+
      # scale_color_manual(values = )
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey90"),
            panel.background = element_rect(fill = "transparent"),
            legend.position = c(0.15,0.8),
            legend.title = element_blank(),
            legend.text = element_text(size = 22),
            legend.key = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 16),
            panel.border = element_blank(),
            axis.line.x = element_line(size = 1, color = "black"),
            axis.line.y = element_line(size = 1, color = "black")) +
      ylab(expression(paste(atop("Suspended sediment load", "(Megaton dry mass)")))) +
      xlab(expression(paste("Discharge (",m^3," x ",10^6,")"))) +
      annotate("text", x = 3500, y = 0.375, label = "Pseudo-R^2 == 0.75", parse = TRUE, size = 7, hjust = 0) +
      annotate("text", x = 3500, y = 0.23, label = "P[Q] < 0.001", parse = TRUE, size = 7, hjust = 0) +
      annotate("text", x = 3500, y = 0.11, label = "P[year] < 0.028", parse = TRUE, size = 7, hjust = 0) +
      annotate("text", x = 3500, y = 0, label = "P[Q * year] < 0.001", parse = TRUE, size = 7, hjust = 0)
  
  
  
maumSSyrPred <- maum %>% 
    select(Y, mtSSAt500Qm3e6pred:mtSSAt4000Qm3e6pred) %>% 
    pivot_longer(cols = mtSSAt500Qm3e6pred:mtSSAt4000Qm3e6pred, names_to = "Flow_m3e6", values_to = "mtSS_pred") %>% 
    mutate(Flow_m3e6 = fct_recode(Flow_m3e6, 
                                  "500" = "mtSSAt500Qm3e6pred",
                                  "1000" = "mtSSAt1000Qm3e6pred",
                                  "2000" = "mtSSAt2000Qm3e6pred",
                                  "3000" = "mtSSAt3000Qm3e6pred",
                                  "4000" = "mtSSAt4000Qm3e6pred"),
           Flow_m3e6 = fct_relevel(Flow_m3e6, c("500", "1000", "2000", "3000", "4000")))
  
  # Fig S8d ----
FigS6d <- ggplot(maumSSyrPred, aes(y = mtSS_pred/1e6, x = Y, color = Flow_m3e6)) +
    geom_line(size = 1.5)+
    scale_color_manual(values = c("#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d"), name = expression(paste("Discharge (",m^3," x ",10^6,")"))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          panel.background = element_rect(fill = "transparent"),
          legend.position = c(0.8,0.84),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 16),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 1, color = "black"),
          axis.line.y = element_line(size = 1, color = "black")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), lim = c(0,1.5)) +
    ylab(expression(paste(atop("Suspended sediment load", "(Megaton dry mass)")))) +
    xlab("Year")
  
  
  
  
  ## Make plot ---- 
  FigS6a.g <- ggplotGrob(FigS6a)
  FigS6b.g <- ggplotGrob(FigS6b)
  FigS6c.g <- ggplotGrob(FigS6c)
  FigS6d.g <- ggplotGrob(FigS6d)
  
  FigS6a.gtf <- gtable_frame(FigS6a.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
  FigS6b.gtf <- gtable_frame(FigS6b.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
  FigS6c.gtf <- gtable_frame(FigS6c.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
  FigS6d.gtf <- gtable_frame(FigS6d.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
  
  FigS6ab.gtf <- gtable_frame(gtable_cbind(FigS6a.gtf, FigS6b.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
  FigS6cd.gtf <- gtable_frame(gtable_cbind(FigS6c.gtf, FigS6d.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
  
  FigS6.gtf <- gtable_rbind(FigS6ab.gtf, FigS6cd.gtf)
  
  png("05_Figures/09_FigS8.png", units = "in", height = 12, width = 16, res = 300)
  grid.newpage()
  grid.draw(FigS6.gtf)
  grid.text("a", x = unit(0.03,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
  grid.text("b", x = unit(0.53,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
  grid.text("c", x = unit(0.03,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
  grid.text("d", x = unit(0.53,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
  dev.off()

# save data ----
  # write.csv(maum, file.path(here::here("04a_generatedDataOnGit"), "09_MaumeeLoadsMarJuneHighFlow.csv"))

# save best SS, Q, and Y models
  # saveRDS(SSmod1, file.path(here::here("04a_generatedDataOnGit"), "09_SS_QQandYglm"))
  
# save image ----
  # save.image(file.path(here::here("03_Rdata"), "09_SStimeseriesFigS6_Rdata"))
  # load(file.path(here::here("03_Rdata"), "09_SStimeseriesFigS6_Rdata"))
