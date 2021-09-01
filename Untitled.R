################
#  CHANGES OVER TIME
################

# DISCHARGE v TIME
ggplot(Mall2, aes(y = L10.Qm3e6_marJun_50per, x = Y)) +
  geom_point()

QvY.l <- lm(L10.Qm3e6_marJun_50per ~ Y, Mall2); 
summary(QvY.l) # year is marginal

# SUSPENDED SEDIMENTS v Q
  # Decline in SS loads over time
  ggplot(Mall2, aes(y = mtSSmarJun_50per, x = Yn)) +
    geom_point() 
  
  #best model - log10
  SSvQintYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per * Yn, Mall2); 
  SSvQpYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yn, Mall2)
  SSvQpYg2.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yg2, Mall2)
  SSvQ.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall2)
  
  # what's the best model - Q+Yn, Q+Yg2
  model.sel(SSvQintYn.l, SSvQpYn.l, SSvQ.l, SSvQpYg2.l)
  summary(SSvQpYn.l)
  plot(SSvQpYn.l)

  # residuals from SS ~ relationship
  Mall2$SS_Qres <- residuals(SSvQ.l)
  
  #checked with gam still linear
  ggplot(Mall2, aes(y = SS_Qres, x = Yn, label= Yn)) +
    geom_text() +
    stat_smooth(method = "lm")
  
# DRP v Q
  DRPvQintYn.l1 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per * Yn, Mall2)
  DRPvQintYn.l2 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per + Yn, Mall2)
  DRPvQintYn.l3 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall2)
  # just curious if SS did better than Y, but they are equal
  DRPvQintYn.l4 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per * SS_Qres, Mall2)
  DRPvQintYn.l5 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per + SS_Qres, Mall2)
  model.sel(DRPvQintYn.l1, DRPvQintYn.l2, DRPvQintYn.l3, DRPvQintYn.l4, DRPvQintYn.l5) # lowest AICc is additive model
  
  # best time model
  summary(DRPvQintYn.l2) # effect of year is marginally significant
  
  # best SS model
  summary(DRPvQintYn.l5) # effect of SS is marginally significant



# PUT PREDICTIONS IN A DATAFRAME
  # Q to predict at
      # for now, need to update to 25% exceedance
      TargetQ <- (Mall2 %>% 
                    ungroup() %>% 
                    # group_by(Yg2) %>% 
                    summarize(Qm3e6_marJun_50per50per = quantile(Qm3e6_marJun_50per, probs = 0.5)))[[1]]

#   # Make dataframe
#       QpredDF <- data.frame(Yn = seq(from = 1975, to = 2019, by = 1),
#                             L10.Qm3e6_marJun_50per = rep(log10(TargetQ), length = 45)) 
#       
#       QpredDF$SSloadAtTarQ = 10^predict(SSvQpYn.l, QpredDF)
#       QpredDF$PsorbAtTarQ = QpredDF$SSloadAtTarQ * quantile(Mall2$gPsorbedgDM, probs = 0.5)
#       
#       ggplot() +
#         geom_point(data = QpredDF, aes(y = SSloadAtTarQ/10000, x = Yn), color = "brown") +
#         geom_point(data = QpredDF, aes(y = PsorbAtTarQ, x = Yn), color = "blue")
# 
# # DRP EXPORTS
#    ## - INCREASE IN DRP OVER TIME
#       ggplot(Mall2, aes(y = L10.mtDRPmarJun_50per, x= L10.Qm3e6_marJun_50per)) +
#         geom_point()
#       
#       
#  
#       
#   # Put prediction in dataframe
#       QpredDF$DRPloadAtTarQ = 10^predict(DRPvQintYn.l2,QpredDF)
#       QpredDF$DRPloadAtTarQwoPsorbloss = QpredDF$DRPloadAtTarQ + QpredDF$PsorbAtTarQ
#       
#       ggplot() +
#         geom_point(data = QpredDF, aes(y = SSloadAtTarQ/10000, x = Yn), color = "brown") +
#         geom_point(data = QpredDF, aes(y = PsorbAtTarQ, x = Yn), color = "blue")+
#         geom_ribbon(data = QpredDF, aes(ymin = DRPloadAtTarQ, ymax = DRPloadAtTarQwoPsorbloss, x = Yn), fill = "green")
      
      
# new try
      Mall3 <- Mall2
      Mall3$Yn1975 <-  1975
      
      # SS PREDICTIONS
      # SS predictions from model above
      Mall3$mtSSmarJun_pred = 10^predict(SSvQpYn.l, Mall3)
      
      # SS for each yearly Q with 1975 concentrations (accounts for Q)
      Mall3$mtSSmarJun_w75seds = 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = Mall3$L10.Qm3e6_marJun_50per, 
                                                                  Yn = rep(1975, length = 42)))
      
      # SS for each yearly Q with 2019 concentrations (accounts for Q)
      Mall3$mtSSmarJun_w2019seds = 10^predict(SSvQpYn.l, data.frame(L10.Qm3e6_marJun_50per = Mall3$L10.Qm3e6_marJun_50per, 
                                                                    Yn = rep(2019, length = 42)))
      
      # P SORPTION PREDICTIONS
      # P sorbed with predicted seds, using median from median values from monte carlo
      Mall3$mtPsorbMarJunP_wPredSeds = Mall3$mtSSmarJun_pred* quantile(Mall2$gPsorbedgDM, probs = 0.5)
      
      # P sorbed with 75 seds
      Mall3$mtPsorbMarJunP_w75seds = Mall3$mtSSmarJun_w75seds* quantile(Mall2$gPsorbedgDM, probs = 0.5)
      
      # P sorbed with 2019 seds
      Mall3$mtPsorbMarJunP_w2019seds = Mall3$mtSSmarJun_w2019seds* quantile(Mall2$gPsorbedgDM, probs = 0.5)
      
      # DIFFERENCES BETWEEN DIFFERENT P SORPTION PREDICTIONS
      # Dif between P sorb with 1975 and year "y" seds
      Mall3$mtPsorbMarJunP_75_Ydiff = Mall3$mtPsorbMarJunP_w75seds - Mall3$mtPsorbMarJunP_wPredSeds
      
      # Dif between P sorb with 1975 and 2019 seds
      Mall3$mtPsorbMarJunP_75_19diff = Mall3$mtPsorbMarJunP_w75seds - Mall3$mtPsorbMarJunP_w2019seds

# PREPARE DATA FOR PLOT
    # GET DRP -Q PREDICTIONS FOR < 2003 AND >2003
      #< 2003
        Mall_L03 <- Mall3 %>% 
          filter(Yn < 2003)
      
        DRP_Q_lm_L03 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall_L03)
        summary(DRP_Q_lm_L03)
        Mall_L03$mtDRPmarJun_predict <- 10^predict(DRP_Q_lm_L03, Mall_L03)
        
    # > 2003
      Mall_G03 <- Mall3 %>% 
        filter(Yn >= 2003) %>% 
        mutate(mtDRPmarJun_w75SSa = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff,
               mtDRPmarJun_w75SSb = mtDRPmarJun_50per - mtPsorbMarJunP_75_19diff,
               L10.mtDRPmarJun_w75SSa = log10(mtDRPmarJun_w75SSa),
               L10.mtDRPmarJun_w75SSb = log10(mtDRPmarJun_w75SSb))
      
      # DRP ~ Q - observed
        DRP_Q_lm_G03 <- lm(L10.mtDRPmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall_G03)
        summary(DRP_Q_lm_G03)
        Mall_G03$mtDRPmarJun_predict <- 10^predict(DRP_Q_lm_G03, Mall_G03)
      
      # DRP ~ Q - 1975 SS - option 1 1975 - year "y"
        DRP_Q_lm_L03w75SSa.lm <- lm(L10.mtDRPmarJun_w75SSa ~ L10.Qm3e6_marJun_50per, Mall_G03)
        summary(DRP_Q_lm_L03w75SSa.lm)
        Mall_G03$mtDRPmarJun_predict_w75SSa <- 10^predict(DRP_Q_lm_L03w75SSa.lm, Mall_G03)
      
      # DRP ~ Q - 1975 SS - option 1 1975 - year "y"
        DRP_Q_lm_L03w75SSb.lm <- lm(L10.mtDRPmarJun_w75SSb ~ L10.Qm3e6_marJun_50per, Mall_G03)
        summary(DRP_Q_lm_L03w75SSb.lm)
        Mall_G03$mtDRPmarJun_predict_w75SSb <- 10^predict(DRP_Q_lm_L03w75SSb.lm, Mall_G03)
        
    
# PLOT
    ggplot()+
      geom_segment(data = Mall_G03, aes(y = mtDRPmarJun_50per, 
                                     yend = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff,
                                     x= Qm3e6_marJun_50per, xend = Qm3e6_marJun_50per),
                   arrow = arrow(length = unit(0.1, "inches")), color = "grey30", size = 0.75, alpha = 0.75) +
      # geom_text(data = Mall_G03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per, label = Yn2), color = "firebrick", fontface = "bold") +
      geom_point(data = Mall_G03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per), color = "firebrick", fontface = "bold", size = 3) +
      geom_point(data = Mall_G03, aes(y = mtDRPmarJun_50per - mtPsorbMarJunP_75_Ydiff, x = Qm3e6_marJun_50per), 
                 fill = "salmon", shape = 21, size = 3) +
      geom_point(data = Mall_L03, aes(y = mtDRPmarJun_50per, x = Qm3e6_marJun_50per), 
                 fill = "dodgerblue", shape = 21, size = 3) +
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
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 12),
            panel.border = element_blank(),
            axis.line.x = element_line(size = 1, color = "black"),
            axis.line.y = element_line(size = 1, color = "black")) +
      geom_point(aes(x = 1000, y = 450), shape = 21, fill = "firebrick", size = 3) +
      geom_point(aes(x = 1000, y = 425), shape = 21, fill = "salmon", size = 3) +
      geom_point(aes(x = 1000, y = 400), shape = 21, fill = "dodgerblue", size = 3) +
      annotate("text", x = 1100, y = 450, label = "≥ 2003 observed", size = 7, color = "black", hjust = 0) +
      annotate("text", x = 1100, y = 425, label = "≥ 2003 predicted with 1975 SS", size = 7, color = "black", hjust = 0) +
      annotate("text", x = 1100, y = 400, label = "< 2003 observed", size = 7, color = "black", hjust = 0) 


      
   