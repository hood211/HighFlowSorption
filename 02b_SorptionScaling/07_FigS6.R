# JMH Feb 2022
# Fig. S6 was Fig.3, hence the names
# a version of Fig. 3 assuming colloid-P:DRP = 0


# libraries ----
library(tidyverse)
library(grid)
library(egg)

# Data ----
maum <- read.csv(file.path(here::here("04a_generatedDataOnGit"),"06_Scale2WVall_MjHf.csv"), row.names = 1) %>% 
          mutate(Yg = ifelse(Y >=2003, "G2003", "L2003"),
                 Yg = fct_relevel(Yg,"L2003", "G2003"))

maumCyanos <- read.csv(file.path(here::here("04a_generatedDataOnGit"),"06_Scale2WVall_MjHf_cyanos.csv"), row.names = 1) %>% 
          filter(var %in% c("cyano_modelMean_obs", "cyano_modelMean_woSorp")) %>% 
          mutate(var = fct_recode(var, "Observed" ="cyano_modelMean_obs", "Without P sorption" = "cyano_modelMean_woSorp"),
                 Yn2 = str_sub(Yn,3))

# Calcs 4 results ----
  ## How much higher was cyanos? ----
  CyanosPerHigher <- read.csv(file.path(here::here("04a_generatedDataOnGit"),"06_Scale2WVall_MjHf_cyanos.csv"), row.names = 1) %>% 
    filter(var %in% c("cyano_modelMean_obs", "cyano_modelMean_woSorp")) %>% 
    pivot_wider(id_cols = Yn, names_from = var, values_from = values) %>% 
    mutate(CyanosPerHigher = cyano_modelMean_woSorp/cyano_modelMean_obs*100)
  
summary(CyanosPerHigher$CyanosPerHigher)

  ## Dif in P sorb before/after 2003 ----
  maum %>% 
    group_by(Yg) %>% 
    summarize_at(vars(perSorbP_50per), list(mean), na.rm = T) %>% 
    pivot_wider(names_from = Yg, values_from = perSorbP_50per) %>% 
    mutate(perHigher = L2003/G2003-1)
  
  ## Per higher with SD
  maum %>% 
    group_by(Yg) %>% 
    mutate(PerHigherWoPsorp = 1-gDRP_MjHf_50per/gDRPwoSorp_MjHf_50per) %>% 
    summarize_at(vars(perSorbP_50per, PerHigherWoPsorp), 
                 list(mean = mean,sd = sd), na.rm = T)

# Fig 3a ----
# Time series of DRP exports with and without sorption
# These are for Mar-June and only high flows
# Note these all use what I call in earlier scripts "potential" sorption
# Which allows > 100% of the DRP to be sorbed. This is reasonable because we are making this estimate based on what is at the bottom of the reach.


Fig3a <- ggplot() +
  geom_ribbon(data = maum, aes(ymin = gDRP_MjHf_50per/1e6, ymax = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per/1e6, x = Y), fill = "gold2", alpha = 40/100) +
  geom_segment(data = maum, aes(y = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per/1e6, yend =  gDRP_MjHf_50per/1e6+5, x = Y, xend=Y), 
               arrow = arrow(length = unit(0.1, "inches")), color = "black", size = 1.15)+
  geom_point(data = maum, aes(y = gDRP_MjHf_50per/1e6, x = Y), fill = "black", shape = 21, size = 3) +
  geom_point(data = maum, aes(y = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per/1e6, x = Y), fill = "grey", shape = 21, size = 3) +
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

# Fig 3b ----
# DRP sorp v. DRP load 
sorpRange0 <- aes(ymin = gPsorbP_MjHf_2.5per/1e6, ymax = gPsorbP_MjHf_97.5per/1e6, x = gDRP_MjHf_50per/1e6)
Fig3b <- ggplot(maum, aes(y = gPsorbP_MjHf_50per/1e6, x = gDRP_MjHf_50per/1e6,  fill = Yg)) +
  geom_errorbar(sorpRange0) +
  geom_point(size = 4, shape = 21) +
  geom_segment(aes(x = 0, xend = 290, y = 0, yend = 290), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 560, y = 0, yend = 560*0.5), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 560, y = 0, yend = 560*0.25), color = "grey", linetype = "dashed") +
  annotate("text", x = 280, y = 300, label = "100%", color = "grey40", size = 6) +
  annotate("text", x = 600, y = 290, label = "50%", color = "grey40", size = 6) +
  annotate("text", x = 600, y = 145, label = "25%", color = "grey40", size = 6) +
  scale_fill_manual(values = c("dodgerblue", "firebrick"), name = NULL, labels = c("< 2003","≥ 2003")) +
  theme_bw() +
  xlab("DRP load (tons P)") +
  ylab("P sorption (tons P)") +
  ylim(0,300) +
  xlim(0,600) +
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
    axis.line.y = element_line(size = 1, color = "black"))


# Fig 3c ----
# estimate of cyanos with and without sorption

Fig3c <- ggplot(maumCyanos,
                aes(y = values, x = var, label = Yn2))+
  geom_boxplot(fill = "grey90", alpha = 0.7, outlier.color = "transparent") + #"#7fbf7b"
  geom_text(position = "jitter", size = 6) +
  # ylim(0,30) +
  theme_bw() +
  ylab(expression(atop("Cyanobacteria bloom", paste("(kiloton dry mass)")))) +
  xlab(NULL) +
  theme(axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 25, hjust = 1),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) 

# Print ----
# NOW THIS IS FIG S6
Fig3a.g <- ggplotGrob(Fig3a)
Fig3b.g <- ggplotGrob(Fig3b)
Fig3c.g <- ggplotGrob(Fig3c)



Fig3a.gtf <- gtable_frame(Fig3a.g,  width = unit(1.8, "null"), height = unit(0.9, "null"))
Fig3b.gtf <- gtable_frame(Fig3b.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))
Fig3c.gtf <- gtable_frame(Fig3c.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))

Fig3bC.gtf <- gtable_frame(gtable_cbind(Fig3b.gtf, Fig3c.gtf),  width = unit(1.8, "null"), height = unit(0.9, "null"))

Fig3.gtf <- gtable_rbind(Fig3a.gtf, Fig3bC.gtf)




png("05_Figures/07_FigS6.png", units = "in", height = 12, width = 14, res = 300)
grid.newpage()
grid.draw(Fig3.gtf)
grid.text("a", x = unit(0.01,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
grid.text("b", x = unit(0.01,"npc"), y = unit(0.52,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
grid.text("c", x = unit(0.51,"npc"), y = unit(0.52,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
dev.off()


# save
# write.csv(maum, file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maum_0col.csv"))
# write.csv(maumCyanos, file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maumCyano_0col.csv"))

# save.image(file.path(here::here("03_Rdata"), "07_FigS6_Rdata"))
# load(file.path(here::here("03_Rdata"), "07_FigS6_Rdata"))
