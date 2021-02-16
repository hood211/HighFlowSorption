library(tidyverse)
library(MuMIn)
library(ggpubr)
library(tibble)
library(ggpmisc)
library(grid)
library(egg)
library(png)

# subset of data for HABS stuff 2003 to 2019
Mhab <- read.csv("04_generatedData/06d_Scale2WVhabs.csv", row.names = 1) %>% 
  mutate(var = as.factor(var))


Mall <- read.csv("04_generatedData/06d_Scale2WVall.csv", row.names = 1) %>% 
  mutate(Yg = as.factor(Yg))


################
# Fig 4A - % sorption v. Y
################

Fig4A <- ggplot() +
  # geom_ribbon(data = Mall, aes(ymin = 0, ymax = gDRPmarJun_50per, x = Y), fill = "lightblue") +
  geom_ribbon(data = Mall, aes(ymin = gDRPmarJun_50per/1e6, ymax = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, x = Y), fill = "pink", alpha = 40/100) +
  geom_point(data = Mall, aes(y = gDRPmarJun_50per/1e6, x = Y), fill = "steelblue", shape = 21, size = 3) +
  geom_line(data = Mall, aes(y = gDRPmarJun_50per/1e6, x = Y), color = "steelblue", alpha = 60/100) +
  geom_point(data = Mall, aes(y = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, x = Y), fill = "red", shape = 21, size = 3) +
  geom_line(data = Mall, aes(y = gDRPmarJun_50per/1e6 + gPsorbMarJunP_50per/1e6, x = Y), color = "red", alpha = 60/100) +
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
  # annotate("text", x = 1975, y = 590, label = "c", size = 8, fontface = "bold") +
  geom_point(aes(y =600, x = 1975), color = "red", size = 5) +
  annotate("text", label = "Without P sorption", y = 600, x = 1977, size = 8, hjust = 0) +
  geom_point(aes(y =565, x = 1975), color = "steelblue", size = 5) +
  annotate("text", label = "Observed", y = 565, x = 1977, size = 8, hjust = 0) 
  


################
# Fig 4A - % sorption v. Y
################
sorpRange2 <- aes(ymin = perSorb2_2.5per*100, ymax = perSorb2_97.5per*100, x = Yn)
# Fig4A <- 
  ggplot(Mall) + 
  stat_smooth(data = Mall, aes(y = perSorb2_50per*100, x = Yn), span = 1, color = "red", fill = "grey80", alpha = 0.5) +
  # geom_point(data = Mall, aes(y = gDRPmarJun_50per/1e6, x = Yn), color = "blue", alpha = 0.4) +
  # geom_line(data = Mall, aes(y = gDRPmarJun_50per/1e6, x = Yn), color = "blue", alpha = 0.4) +
  # geom_point(data = Mall, aes(y = Qm3marJun_50per/1e7, x = Yn), color = "red", alpha = 0.4) +
  # geom_line(data = Mall, aes(y = Qm3marJun_50per/1e7, x = Yn), color = "red", alpha = 0.4) +
  # geom_point(data = Mall, aes(y = gSSmarJun_50per/1e9, x = Yn), color = "brown", alpha = 0.4) +
  # geom_line(data = Mall, aes(y = gSSmarJun_50per/1e9, x = Yn), color = "brown", alpha = 0.4) +
  geom_errorbar(sorpRange2)+
  geom_point(aes(y = perSorb2_50per*100, x = Yn), shape = 21, fill = "grey30", size = 5) +
  ylab("% DRP sorbed") +
  xlab("Year") +
  theme_bw() +
  annotate("text", x = 1975, y = 450, label = "a", size = 8, fontface = "bold") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black"))

  # trends in DRP, Q, and SS loads  
  ggplot(Mall %>% 
           select(Yn, gDRPmarJun_50per, Qm3marJun_50per, gSSmarJun_50per) %>% 
           pivot_longer(cols = c(gDRPmarJun_50per:gSSmarJun_50per), names_to = "vars", values_to = "values"),
         aes(y = values, x = Yn)) +
    geom_point() +
    geom_line(alpha = 0.6) +
    stat_smooth() +
    facet_grid(vars ~., scales = "free_y") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12))





################
# Cyano fig
################

Fig4D <- ggplot(Mhab %>% 
                  filter(var == "cyanoIndex" | var == "cyanoIndexWoS") %>% 
                  mutate( Yn2 = str_sub(Yn,3),
                          var = fct_recode(var, "Observed" ="cyanoIndex", "Without P sorption" = "cyanoIndexWoS")),
                aes(y = values, x = var,label = Yn2))+
    geom_boxplot(fill = "grey80", alpha = 0.7) + #"#7fbf7b"
    geom_text(position = "jitter", size = 6) +
    scale_y_log10() +
    # annotate("text", x = 0.5, y = 90, label = "e", size = 8, fontface = "bold") +
    theme_bw() +
    ylab("Max cyanobacterial index") +
    xlab(NULL) +
    theme(axis.title = element_text(size = 24),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 20),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(size = 1, color = "black"),
          axis.line.y = element_line(size = 1, color = "black"))

# effect size
  MhabW <- Mhab %>% 
    pivot_wider(names_from = "var", values_from = "values")
  
  MhabWs <- MhabW %>% 
    summarize(across(cyanoIndex:cyanoIndexWoS, mean, na.rm = TRUE)) %>% 
    mutate(effectSize = 1-cyanoIndex/cyanoIndexWoS)
  
# just a look at the relationship between cyano Index and sorption
ggplot(Mhab %>% 
         pivot_wider(names_from = "var", values_from = "values") ,
       aes(y = cyanoIndex,x = mtPsorbMarJunP_50per)) +
  geom_point()



################
# SS is drive sorption
################

# THis is primarily driven by sediment concentrations
ggplot(Mall, aes(y = gPsorbMarJunP_50per/1e6, x = gSSmarJun_50per/1e6, label = Y)) +
  geom_text()   

# so what drive SS?
Mall2 <- Mall %>% 
  mutate(Yg2 = ifelse(Yn < 2003, "<2003", ">2003"),
         mtSSmarJun_50per = gSSmarJun_50per/1e6, # Kilotonne or 1 billion g
         Qm3e6_marJun_50per = Qm3marJun_50per/1e6,
         L10.mtSSmarJun_50per = log10(mtSSmarJun_50per),
         L10.Qm3e6_marJun_50per = log10(Qm3e6_marJun_50per),
         Yn2 = str_sub(Yn,3))

################
# Differences before and after 2003
################
Mall2 %>% 
  group_by(Yg2) %>% 
  summarize_at(vars(perSorb2_50per), list(mean = mean), na.rm = TRUE)

  

################
# SS is drive sorption _more
################
# First let's see if SS changes over time
#best model - log10
  SSvQintYg.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per * Yg, Mall2); 
  SSvQintYg2.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per * Yg2, Mall2); 
  SSvQintYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per * Yn, Mall2); 
  SSvQpYg.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yg, Mall2)
  SSvQpYg2.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yg2, Mall2)
  SSvQpYn.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per + Yn, Mall2)
  SSvQ.l <- lm(L10.mtSSmarJun_50per ~ L10.Qm3e6_marJun_50per, Mall2)

# what's the best model - Q+Yn, Q+Yg2
model.sel(SSvQintYg.l, SSvQintYn.l, SSvQintYg2.l, SSvQpYg.l, SSvQpYn.l, SSvQpYg2.l, SSvQ.l)
summary(SSvQpYn.l)
# summary(SSvQpYg2.l)

# Decline in SS loads over time
ggplot(Mall2, aes(y = mtSSmarJun_50per, x = Yn)) +
  geom_point() +
  # scale_x_log10() +
  stat_smooth(method = "lm")

ggplot(Mall2, aes(y = gDRPmarJun_50per/1e6, x = Yn)) +
  geom_point() +
  geom_line()

Mall2 %>% 
  ungroup() %>% 
  # group_by(Yg2) %>% 
  summarize(Qm3e6_marJun_50per50per = quantile(Qm3e6_marJun_50per, probs = 0.5))
SSload1999 <- 10^predict(SSvQpYn.l, as.data.frame(as.matrix(cbind(L10.Qm3e6_marJun_50per = log10(2420), Yn = 1999))))
SSload2019 <- 10^predict(SSvQpYn.l, as.data.frame(as.matrix(cbind(L10.Qm3e6_marJun_50per = log10(2420), Yn = 2019))))
PerDeclineSSDecade <- (SSload2019 -SSload1999)/mean(c(SSload2019,SSload1999))/2
AbsDeclineSSDecade <- (SSload2019 -SSload1999)/2

# decline in P sorption from 2003/2019
mtPsorb1999 <- SSload1999*quantile(Mall2$gPsorbedgDM_50per, probs = 0.5)
mtPsorb2019 <- SSload2019*quantile(Mall2$gPsorbedgDM_50per, probs = 0.5)
PerDeclinePsorbDecade <- (mtPsorb1999 - mtPsorb2019)/mean(c(mtPsorb1999, mtPsorb2019))/2
PerDeclinePsorbDecade2 <- (mtPsorb1999 - mtPsorb2019)/mtPsorb1999/2
AbsDeclinePsorbDecade <- (mtPsorb1999 - mtPsorb2019)/2
AbsDeclinePsorbYr <- AbsDeclinePsorbDecade/10
ChangeInHabs <- 0.48 * 10^(AbsDeclinePsorbYr * 0.00387) # says a * 10^-3 = 0.001?
ChangeInHabs*10

SSpredFake <- as.data.frame(cbind(L10.Qm3e6_marJun_50per = rep(seq(min(Mall2$L10.Qm3e6_marJun_50per), max(Mall2$L10.Qm3e6_marJun_50per), length = 1000), times = 2),
                            Yg2 = c(rep("<2003", times = 1000), rep(">2003", times = 1000)))) %>% 
        mutate(L10.Qm3e6_marJun_50per = as.numeric(L10.Qm3e6_marJun_50per))

SSpredFake <- SSpredFake %>% 
                mutate(L10.Qm3e6_marJun_50per = as.numeric(L10.Qm3e6_marJun_50per),
                       L10.mtSSmarJun_50per = predict(SSvQpYg2.l, SSpredFake),
                       Qm3e6_marJun_50per = 10^(L10.Qm3e6_marJun_50per),
                       mtSSmarJun_50per = 10^(L10.mtSSmarJun_50per))


# check to make sure relationship between res and year is linear
# Man, that's legit!
Mall3 <- Mall2 %>% 
  mutate(SSvQ.lres = residuals(SSvQ.l),
         SSflowWeightedMean_mgL = log10(mtSSmarJun_50per/Qm3e6_marJun_50per*1000))



Fig4Cinset <- ggplot(Mall3, aes(y = SSflowWeightedMean_mgL, x = Yn)) +
  stat_smooth(method = "lm", color = "black", fill = "grey80") +
  geom_point(shape = 21, size = 2, fill = "grey30") +
    # ylab(expression(paste(log[10]," Flow-weighted mean concentration (mg DM ",L^-1,")"))) +
  ylab(expression(paste(log[10]," FWMC (mg DM ",L^-1,")"))) +
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
Fig4B <- ggplot(Mall2, aes(y = gPsorbMarJunP_50per/1e6, x = gDRPmarJun_50per/1e6, color = Yg2)) +
  geom_point(size = 4) +
  # geom_text() +
  geom_errorbar(sorpRange0) +
  geom_segment(aes(x = 0, xend = 300, y = 0, yend = 300), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 500, y = 0, yend = 500*0.5), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 500, y = 0, yend = 500*0.25), color = "grey", linetype = "dashed") +
  annotate("text", x = 280, y = 300, label = "100%", color = "grey40", size = 6) +
  annotate("text", x = 500, y = 260, label = "50%", color = "grey40", size = 6) +
  annotate("text", x = 500, y = 135, label = "25%", color = "grey40", size = 6) +
  # annotate("text", x = 0, y = 300, label = "d", size = 8, fontface = "bold") +
  scale_color_manual(values = c("#af8dc3", "#7fbf7b"), name = "Year", labels = c("≤ 2002","≥ 2003")) +
  # geom_abline(intercept = 0, slope =1, color = "grey") +
  # geom_abline(intercept = 0, slope =0.5, color = "grey") +
  # geom_abline(intercept = 0, slope = 0.25, color = "grey") +
  theme_bw() +
  xlab("DRP load (tons P)") +
  ylab("P sorption (tons P)") +
  ylim(0,300) +
  xlim(0,500) +
  theme(
    legend.position = c(0.2,0.9),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 22),
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
# SS is drive sorption
################


data.tb <- tibble(x = 5100, y = 50,
                  plot = list(Fig4Cinset))

Fig4C <- ggplot() + 
  geom_point(data = Mall2, 
             aes(y = mtSSmarJun_50per/1000,x = Qm3e6_marJun_50per, color = Yg2), size = 4) +
  geom_line(data = SSpredFake, 
            aes(y = mtSSmarJun_50per/1000, x = Qm3e6_marJun_50per, color = Yg2)) +
  geom_plot(data = data.tb, aes(x,y, label = plot)) +
  scale_x_log10(expression(paste("Discharge (",m^3," x ", 10^6,")"))) +
  scale_y_log10(expression(paste("SS load (", 10^3," tons DM)")),
  sec.axis = sec_axis(~ . * quantile(Mall$gPsorbedgDM_50per, probs = 0.5), #conver to kT DM then metric ton P
  "P sorption (tons P)"), limits =c(50,1300)) +
  scale_color_manual(values = c("#af8dc3", "#7fbf7b"), name = "Year", labels = c("≤ 2002","≥ 2003")) +
  # annotate("text", x = 1000, y = 1250, label = "f", size = 8, fontface = "bold") +
  theme_bw() +
  theme(
    # legend.position = c(0.22,0.9),
    legend.position = "none",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
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
# Assemble PLOT
################



Fig4A.g <- ggplotGrob(Fig4A)
Fig4B.g <- ggplotGrob(Fig4B)
Fig4C.g <- ggplotGrob(Fig4D)
Fig4D.g <- ggplotGrob(Fig4C)


Fig4A.gtf <- gtable_frame(Fig4A.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4B.gtf <- gtable_frame(Fig4B.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4C.gtf <- gtable_frame(Fig4C.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4D.gtf <- gtable_frame(Fig4D.g,  width = unit(0.9, "null"), height = unit(0.9, "null"))

# Fig4ABC.gtf <- gtable_frame(gtable_cbind(Fig4AB.gtf, Fig4C.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4AB.gtf <- gtable_frame(gtable_cbind(Fig4A.gtf, Fig4B.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))
Fig4CD.gtf <- gtable_frame(gtable_cbind(Fig4C.gtf, Fig4D.gtf),  width = unit(0.9, "null"), height = unit(0.9, "null"))

Fig4.gtf <- gtable_rbind(Fig4AB.gtf, Fig4CD.gtf)




cairo_pdf("05_Figures/08_Fig4.pdf", height = 12, width = 15)
grid.newpage()
grid.draw(Fig4.gtf)
grid.text("a", x = unit(0.01,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.01,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.48,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("d", x = unit(0.49,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()
  

  ################
  # Effects sizes
  ################
  Mass3Sum <- Mall3 %>% 
    group_by(Yg2) %>% 
    summarize_at(vars(gPsorbedgDM_50per:DRPMarJunWvWOsorp_50per, mtSSmarJun_50per:L10.Qm3e6_marJun_50per), list(mean = mean), na.rm = TRUE)
  
  # DRP loads without sorption
  Mass3Sum$DRPMarJunWvWOsorp_50per_mean
  
  # HABs w/o sorption
  
  ################
  # SAVE/LOAD
  ################
  
  save.image("03_Rdata/08_MaumeeFigs_Rdat")
  # load("08_MaumeeFigs_Rdat")
  
  